import random

configfile: "0_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz")
ASSEMBLIES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2")
READS = SAMPLES

# Iterate mapping sizes R2..R9 (skip single read-set and skip all-reads).
MAX_ITER = min(len(READS) - 1, 9)
ITERATION_SIZES = [i for i in range(2, MAX_ITER + 1)]
ITERATIONS = [f"R{i}" for i in ITERATION_SIZES]
RANDOM_SEED = 42

def choose_reads_for_iteration(size):
    """Return dict mapping assembly -> list of read datasets (self + random others) of given size."""
    selection = {}
    extras = size - 1
    for assembly in ASSEMBLIES:
        other_reads = [r for r in READS if r != assembly]
        rng = random.Random(f"{assembly}-{size}-{RANDOM_SEED}")
        chosen = rng.sample(other_reads, extras)
        selection[assembly] = [assembly] + chosen
    return selection

# Precompute selections per iteration
ASSEMBLY_READS = {f"R{size}": choose_reads_for_iteration(size) for size in ITERATION_SIZES}
READ_ASSEMBLY_ITER_PAIRS = {
    (iter_label, r, a)
    for iter_label, mapping in ASSEMBLY_READS.items()
    for a, reads in mapping.items()
    for r in reads
}

rule all:
    input:
        expand(f"{WORKDIR}/4_random_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth", iteration=ITERATIONS, assembly=ASSEMBLIES),
        expand(f"{WORKDIR}/4_random_coverage/{{iteration}}/bowtie2/{{assembly}}_maxbin.depth", iteration=ITERATIONS, assembly=ASSEMBLIES)

rule assembly_map:
    input:
        index=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2",
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_2.fq.gz"
    output:
        f"{WORKDIR}/4_random_coverage/{{iteration}}/bowtie2/{{reads}}_vs_{{assembly}}.bam"
    benchmark:
        f"{WORKDIR}/4_random_coverage/{{iteration}}/benchmark/assembly_map/{{reads}}_vs_{{assembly}}.txt"
    params:
        basename=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.iteration}: reads {wildcards.reads} to assembly {wildcards.assembly}..."
    run:
        key = (wildcards.iteration, wildcards.reads, wildcards.assembly)
        if key not in READ_ASSEMBLY_ITER_PAIRS:
            raise ValueError(f"Combination {wildcards.iteration} {wildcards.reads} vs {wildcards.assembly} not selected for random coverage.")
        shell("""
            module load bowtie2/2.4.2 samtools/1.21
            bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """)


rule assembly_map_depth:
    input:
        lambda wildcards: expand(
            f"{WORKDIR}/4_random_coverage/{wildcards.iteration}/bowtie2/{{reads}}_vs_{wildcards.assembly}.bam",
            reads=ASSEMBLY_READS[wildcards.iteration][wildcards.assembly],
        )
    output:
        metabat2=f"{WORKDIR}/4_random_coverage/{{iteration}}/bowtie2/{{assembly}}_metabat.depth",
        maxbin2=f"{WORKDIR}/4_random_coverage/{{iteration}}/bowtie2/{{assembly}}_maxbin.depth"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Summarising coverage {wildcards.iteration} for assembly {wildcards.assembly}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """
