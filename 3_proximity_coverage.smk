import csv
import os
import re
from collections import defaultdict
from functools import lru_cache

configfile: "0_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]
MASH_DISTANCES = config.get("mash_distances", f"{WORKDIR}/3_proximity_coverage/mash/distances.tsv")

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz")
ASSEMBLIES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2")
READS = SAMPLES

# Iterate mapping sizes R2..R9 (skip single read-set and skip all-reads).
MAX_ITER = min(len(READS) - 1, 9)
ITERATION_SIZES = [i for i in range(2, MAX_ITER + 1)]
ITERATIONS = [f"R{i}" for i in ITERATION_SIZES]


def ordered_reads_by_proximity(distances, assembly):
    """Return read sets ordered by increasing Mash distance to the focal assembly's read set."""
    others = [
        (distances[assembly].get(read, float("inf")), read)
        for read in READS
        if read != assembly
    ]
    others.sort(key=lambda x: (x[0], x[1]))
    return [assembly] + [read for _, read in others]


@lru_cache()
def load_mash_distances(path):
    dist = defaultdict(dict)
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            q = row["query"]
            r = row["reference"]
            d = float(row["distance"])
            dist[q][r] = d
            dist[r][q] = d
    for sample in READS:
        dist[sample].setdefault(sample, 0.0)
    return dist


def select_reads(iter_label, assembly):
    """Return ordered read list for iteration label (R2..R9) for the given assembly."""
    size = int(re.sub(r"[^0-9]", "", iter_label))
    ckpt = checkpoints.mash_distances.get()
    dist_path = ckpt.output[0]
    distances = load_mash_distances(dist_path)
    ordered = ordered_reads_by_proximity(distances, assembly)
    return ordered[:size]


rule all:
    input:
        lambda wildcards: checkpoints.mash_distances.get().output[0],
        expand(f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth", iteration=ITERATIONS, assembly=ASSEMBLIES),
        expand(f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_maxbin.depth", iteration=ITERATIONS, assembly=ASSEMBLIES)


rule mash_sketch_reads:
    input:
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/3_proximity_coverage/mash/{{sample}}.msh"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(4 * 1024, int(input.size_mb * 2) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 50) * 2 ** (attempt - 1))
    message: "Creating Mash sketch for {wildcards.sample}..."
    shell:
        """
        module load mash/2.3
        mash sketch -r -k 21 -s 10000 -o {output} {input.r1} {input.r2}
        """


checkpoint mash_distances:
    input:
        expand(f"{WORKDIR}/3_proximity_coverage/mash/{{sample}}.msh", sample=READS)
    output:
        MASH_DISTANCES
    threads: 2
    resources:
        mem_mb=lambda wildcards, input, attempt: max(2 * 1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(10, int(input.size_mb / 200) * 2 ** (attempt - 1))
    message: "Calculating Mash distances across read sets..."
    shell:
        """
        module load mash/2.3
        mash dist {input} | awk 'BEGIN{{OFS="\\t";print "query","reference","distance","pvalue","shared"}}{{print $1,$2,$3,$4,$5}}' > {output}
        """

rule assembly_map:
    input:
        index=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2",
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_2.fq.gz",
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{reads}}_vs_{{assembly}}.bam"
    benchmark:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/benchmark/assembly_map/{{reads}}_vs_{{assembly}}.txt"
    params:
        basename=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8 * 1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.iteration} (proximity-ranked): reads {wildcards.reads} to assembly {wildcards.assembly}..."
    run:
        selected = select_reads(wildcards.iteration, wildcards.assembly)
        if wildcards.reads not in selected:
            raise ValueError(f"{wildcards.reads} not selected for {wildcards.iteration} vs {wildcards.assembly}")
        shell("""
            module load bowtie2/2.4.2 samtools/1.21
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """)


rule assembly_map_depth:
    input:
        bams=lambda wildcards: expand(
            f"{WORKDIR}/3_proximity_coverage/{wildcards.iteration}/bowtie2/{{reads}}_vs_{wildcards.assembly}.bam",
            reads=select_reads(wildcards.iteration, wildcards.assembly),
        )
    output:
        metabat2=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth",
        maxbin2=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_maxbin.depth"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8 * 1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Summarising proximity coverage {wildcards.iteration} for assembly {wildcards.assembly}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input.bams}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """
