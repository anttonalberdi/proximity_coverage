import os

configfile: "0_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]
READS = config["reads"]
REFERENCE = config["reference_genome"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{READS}/{{sample}}_1.fq.gz")

rule all:
    input:
        expand(f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.rev.2.bt2", sample=SAMPLES)

rule fastp:
    input:
        r1=f"{READS}/{{sample}}_1.fq.gz",
        r2=f"{READS}/{{sample}}_2.fq.gz"
    output:
        r1=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}_2.fq.gz",
        html=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}.html",
        json=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}.json"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 1024 * 3) * 2 ** (attempt - 1))
    message: "Quality-filtering sample {wildcards.sample}..."
    shell:
        """
        module load fastp/0.23.4
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --trim_poly_g \
            --trim_poly_x \
            --low_complexity_filter \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        """

# Build a Bowtie2 index from the given reference genome FASTA file.  
# Produces the index files required for read mapping and saves a copy of the reference sequence.

rule reference_index:
    input:
        REFERENCE
    output:
        index = f"{WORKDIR}/0_preprocessing/reference/{REF_BASENAME}.rev.1.bt2"
    params:
        refdir = f"{WORKDIR}/0_preprocessing/reference/",
        basename = f"{WORKDIR}/0_preprocessing/reference/{REF_BASENAME}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Indexing reference genome..."
    shell:
        """
        module load bowtie2/2.4.2
        cp {input} {params.refdir}
        bowtie2-build {input} {params.basename}
        cat {input} > {params.basename}.fna
        """

# Align quality-filtered paired-end reads to the corresponding reference genome using Bowtie2.  
# The alignments are converted to BAM format and sorted with samtools for downstream analyses.

rule reference_map:
    input:
        index = rules.reference_index.output.index,
        r1=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/fastp/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/0_preprocessing/bowtie2/{{sample}}.bam"
    params:
        basename = f"{WORKDIR}/0_preprocessing/reference/{REF_BASENAME}"
    threads: 16
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} against reference genome..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} -p {threads} | samtools view -bS - | samtools sort -o {output}
        """

# Split mapped BAM files into metagenomic (unmapped) and host (mapped) read sets.  
# Outputs paired FASTQ files for unmapped reads, a host-only BAM, and text files counting reads/bases in each category.

rule extract_metagenomic_reads:
    input:
        rules.reference_map.output
    output:
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_2.fq.gz"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 5) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 100) * 2 ** (attempt - 1))
    message: "Extracting metagenomic reads of {wildcards.sample}..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        samtools view -b -f12 -@ {threads} {input} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -
        """

rule assembly:
    input:
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.fna"
    params:
        outputdir=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1020*1024,max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 20) * 2 ** (attempt - 1))
    message: "Assembling {wildcards.sample}..."
    shell:
        """
        module load megahit/1.2.9
        rm -rf {params.outputdir}

        megahit \
            -t {threads} \
            --verbose \
            --min-contig-len 1500 \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outputdir}
        mv {params.outputdir}/final.contigs.fa {output}
        """

rule assembly_index:
    input:
        f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.fna"
    output:
        index=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.rev.2.bt2"
    params:
        basename=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Indexing assembly {wildcards.sample}..."
    shell:
        """
        module load bowtie2/2.4.2
        bowtie2-build {input} {params.basename}
        """