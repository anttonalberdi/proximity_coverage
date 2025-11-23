import os

configfile: "0_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz")
ASSEMBLIES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2")
READS = SAMPLES

rule all:
    input:
        expand(f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_metabat.depth", assembly=ASSEMBLIES),
        expand(f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_maxbin.depth", assembly=ASSEMBLIES)

rule assembly_map:
    input:
        index=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.rev.2.bt2",
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{reads}}_2.fq.gz"
    output:
        f"{WORKDIR}/2_all_coverage/bowtie2/{{reads}}_vs_{{assembly}}.bam"
    benchmark:
        f"{WORKDIR}/2_all_coverage/benchmark/assembly_map/{{reads}}_vs_{{assembly}}.txt"
    params:
        basename=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping reads from {wildcards.reads} to assembly {wildcards.assembly}..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule assembly_map_depth:
    input:
        lambda wildcards: expand(f"{WORKDIR}/2_all_coverage/bowtie2/{{reads}}_vs_{wildcards.assembly}.bam", reads=READS)
    output:
        metabat2=f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_metabat.depth",
        maxbin2=f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_maxbin.depth"
    benchmark:
        f"{WORKDIR}/2_all_coverage/benchmark/assembly_map_depth/{{assembly}}.txt"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 100) * 2 ** (attempt - 1)))
    message: "Calculating mapping stats across all reads for assembly {wildcards.assembly}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """