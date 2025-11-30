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
ITERATIONS = [f"P{i}" for i in ITERATION_SIZES]


def ordered_reads_by_proximity(distances, assembly):
    """Return read sets ordered by increasing Mash distance to the focal assembly's read set.

    If the assembly has no matching read set (e.g. co-assemblies), fall back to the closest
    available component samples; if none match, return all reads alphabetically.
    """
    if assembly in READS:
        others = [
            (distances[assembly].get(read, float("inf")), read)
            for read in READS
            if read != assembly
        ]
        others.sort(key=lambda x: (x[0], x[1]))
        return [assembly] + [read for _, read in others]

    # For co-assemblies like A_B, use the min distance to any component that exists in READS.
    parts = [p for p in re.split(r"[^\w]+", assembly) if p in READS]
    if parts:
        scored = []
        for read in READS:
            dists = [distances[p].get(read, float("inf")) for p in parts]
            scored.append((min(dists) if dists else float("inf"), read))
        scored.sort(key=lambda x: (x[0], x[1]))
        return [read for _, read in scored]

    # No matching components; default to alphabetical order
    return sorted(READS)


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
        expand(f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_drep/{{assembly}}/data_tables/genomeInformation.csv", iteration=ITERATIONS, assembly=ASSEMBLIES),
        #expand(f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_drep/{{assembly}}/data_tables/genomeInformation.csv", iteration=ITERATIONS, assembly=ASSEMBLIES)


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

rule assembly_map_depth:
    input:
        bams=lambda wildcards: expand(
            f"{WORKDIR}/2_all_coverage/bowtie2/{{reads}}_vs_{{assembly}}.bam",
            reads=select_reads(wildcards.iteration, wildcards.assembly),
            assembly=wildcards.assembly,
        )
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8 * 1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Summarising proximity coverage {wildcards.iteration} for assembly {wildcards.assembly}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bams}
        """

rule maxbin_depth:
    input:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_{{reads}}_maxbin.depth"
    threads: 1
    params:
        idx = lambda wildcards: READS.index(wildcards.reads)
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    shell:
        r"""
        set -euo pipefail

        depthfile="{input}"

        # Column: 1=contig, 2=len, 3=totalAvgDepth, then 4,6,8,... are sample depths
        col=$(( 4 + 2*{params.idx} ))

        awk -v col="$col" 'NR>1 && NF>=col {{ print $1 "\t" $col+0 }}' "$depthfile" > "{output}"
        """

rule metabat2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.fna",
        depth=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2/{{assembly}}.tsv"
    params:
        basename=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} ({wildcards.iteration}) using metabat2..."
    shell:
        """
        module load metabat2/2.17
        metabat2 -i {input.assembly} -a {input.depth} -o {params.basename} -m 1500 --saveCls

        #Standardise bin names
        dir="$(dirname {params.basename})"
        base="$(basename {params.basename})"

        i=0
        for f in "$dir"/"$base".*.fa; do
            [ -e "$f" ] || continue
            i=$((i+1))
            new=$(printf "%s/MET_{wildcards.iteration}_%s_%03d.fna" "$dir" "$base" "$i")
            mv "$f" "$new"
        done

        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "*$(basename {params.basename})_*.fna" | sort > {output}
        """

rule metabat2_checkm:
    input:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_checkm/{{assembly}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/3_proximity_coverage/{wildcards.iteration}/metabat2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_checkm/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of MetaBAT2 bins for {wildcards.assembly} ({wildcards.iteration}) using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fna -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "genome","completeness","contamination"; next}} {{print $1".fna",$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule metabat2_drep:
    input:
        genomes=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2/{{assembly}}.tsv",
        genomeinfo=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_checkm/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_drep/{{assembly}}/data_tables/genomeInformation.csv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/3_proximity_coverage/{wildcards.iteration}/metabat2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/metabat2_drep/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Dereplicating MetaBAT2 bins for {wildcards.assembly} ({wildcards.iteration}) at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -g {input.genomes} -p {threads} -pa 0.95 --genomeInfo {input.genomeinfo}
        """

rule maxbin2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.fna",
        depth=lambda wildcards: expand(
            f"{WORKDIR}/3_proximity_coverage/{{iteration}}/bowtie2/{{assembly}}_{{reads}}_maxbin.depth",
            iteration=wildcards.iteration,
            assembly=wildcards.assembly,
            reads=READS,
        )
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2/{{assembly}}.tsv"
    params:
        basedir=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2/{{assembly}}",
        basename=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2/{{assembly}}/{{assembly}}",
        abund=lambda wildcards, input: " ".join(f"-abund {f}" for f in input.depth)
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 3) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} ({wildcards.iteration}) using maxbin2..."
    shell:
        """
        module load bowtie2/2.4.2 fraggenescan/1.32 maxbin2/2.2.7 hmmer/3.3.2
        rm -rf {params.basedir}
        mkdir -p {params.basedir}
        /opt/shared_software/shared_envmodules/conda/maxbin2-2.2.7/bin/run_MaxBin.pl -contig {input.assembly} {params.abund} -max_iteration 10 -out {params.basename} -min_contig_length 1500
    
        #Standardise bin names
        dir="$(dirname {params.basename})"
        base="$(basename {params.basename})"

        i=0
        for f in "$dir"/"$base".*.fasta; do
            [ -e "$f" ] || continue
            i=$((i+1))
            new=$(printf "%s/MAX_{wildcards.iteration}_%s_%03d.fna" "$dir" "$base" "$i")
            mv "$f" "$new"
        done

        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "*$(basename {params.basename})_*.fna" | sort > {output}

        # Fail if output is empty
        if [ ! -s {output} ]; then
            echo "ERROR: No bins produced for assembly {wildcards.assembly} (iteration {wildcards.iteration}); {output} is empty." >&2
            rm -f {output}
            exit 1
        fi
        """

rule maxbin2_checkm:
    input:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_checkm/{{assembly}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/3_proximity_coverage/{wildcards.iteration}/maxbin2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_checkm/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of Maxbin2 bins for {wildcards.assembly} ({wildcards.iteration}) using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fna -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "genome","completeness","contamination"; next}} {{print $1".fna",$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule maxbin2_drep:
    input:
        genomes=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2/{{assembly}}.tsv",
        genomeinfo=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_checkm/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_drep/{{assembly}}/data_tables/genomeInformation.csv"
    params:
        outdir=f"{WORKDIR}/3_proximity_coverage/{{iteration}}/maxbin2_drep/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb * 1000) * 2 ** (attempt - 1)))
    message: "Dereplicating Maxbin2 bins for {wildcards.assembly} ({wildcards.iteration}) at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -g {input.genomes} -p {threads} -pa 0.95 --genomeInfo {input.genomeinfo}
        """
