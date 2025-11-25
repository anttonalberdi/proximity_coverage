import os

configfile: "0_preprocessing.yaml"

# Config variables
WORKDIR = config["workdir"]

# List genome and target wildcards
SAMPLES, = glob_wildcards(f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz")

rule all:
    input:
        expand(f"{WORKDIR}/1_single_coverage/metabat2_checkm/{{sample}}/quality_report.tsv", sample=SAMPLES),
        expand(f"{WORKDIR}/1_single_coverage/maxbin2_checkm/{{sample}}/quality_report.tsv", sample=SAMPLES)

rule assembly_map:
    input:
        index=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.rev.2.bt2",
        r1=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_1.fq.gz",
        r2=f"{WORKDIR}/0_preprocessing/bowtie/{{sample}}_2.fq.gz"
    output:
        f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}.bam"
    benchmark:
        f"{WORKDIR}/1_single_coverage/benchmark/{{sample}}.txt"
    params:
        basename=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 3) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Mapping {wildcards.sample} reads to assembly..."
    shell:
        """
        module load bowtie2/2.4.2 samtools/1.21
        bowtie2 -x {params.basename} -1 {input.r1} -2 {input.r2} | samtools view -bS - | samtools sort -o {output}
        """

rule assembly_map_depth:
    input:
        f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}.bam"
    output:
        metabat2=f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}_metabat.depth",
        maxbin2=f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}_maxbin.depth"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 1000) * 2 ** (attempt - 1)))
    message: "Calculating mapping states of assembly {wildcards.sample}..."
    shell:
        """
        module load metabat2/2.17
        jgi_summarize_bam_contig_depths --outputDepth {output.metabat2} {input}
        cut -f1,3 {output.metabat2} | tail -n+2 > {output.maxbin2}
        """

rule metabat2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.fna",
        depth=f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}_metabat.depth"
    output:
        f"{WORKDIR}/1_single_coverage/metabat2/{{sample}}.tsv"
    params:
        basename=f"{WORKDIR}/1_single_coverage/metabat2/{{sample}}/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.sample} using metabat2..."
    shell:
        """
        module load metabat2/2.17
        metabat2 -i {input.assembly} -a {input.depth} -o {params.basename} -m 1500 --saveCls

        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "$(basename {params.basename}).*.fa" | sort > {output}
        """

rule metabat2_checkm:
    input:
        f"{WORKDIR}/1_single_coverage/metabat2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/1_single_coverage/metabat2_checkm/{{sample}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/1_single_coverage/metabat2/{wildcards.sample}",
        outdir=f"{WORKDIR}/1_single_coverage/metabat2_checkm/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of MetaBAT2 bins for {wildcards.sample} using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fa -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "Name","Completeness","Contamination"; next}} {{print $1,$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule metabat2_drep:
    input:
        f"{WORKDIR}/1_single_coverage/metabat2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/1_single_coverage/metabat2_drep/{{sample}}/dereplicated_genomes.csv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/1_single_coverage/metabat2/{wildcards.sample}",
        outdir=f"{WORKDIR}/1_single_coverage/metabat2_drep/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Dereplicating MetaBAT2 bins for {wildcards.sample} at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3 checkm-genome/1.2.3
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        dRep dereplicate {params.outdir} -g {input} -p {threads} -pa 0.95
        """

rule maxbin2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.fna",
        depth=f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}_maxbin.depth"
    output:
        f"{WORKDIR}/1_single_coverage/maxbin2/{{sample}}.tsv"
    params:
        basedir=f"{WORKDIR}/1_single_coverage/maxbin2/{{sample}}",
        basename=f"{WORKDIR}/1_single_coverage/maxbin2/{{sample}}/{{sample}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 3) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.sample} using maxbin2..."
    shell:
        """
        module load maxbin2/2.2.7 hmmer/3.3.2
        rm -rf {params.basedir}
        mkdir -p {params.basedir}
        /opt/shared_software/shared_envmodules/conda/maxbin2-2.2.7/bin/run_MaxBin.pl -contig {input.assembly} -abund {input.depth} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        
        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "$(basename {params.basename}).*.fasta" | sort > {output}
        """

rule maxbin2_checkm:
    input:
        f"{WORKDIR}/1_single_coverage/maxbin2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/1_single_coverage/maxbin2_checkm/{{sample}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/1_single_coverage/maxbin2/{wildcards.sample}",
        outdir=f"{WORKDIR}/1_single_coverage/maxbin2_checkm/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of Maxbin2 bins for {wildcards.sample} using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fasta -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "Name","Completeness","Contamination"; next}} {{print $1,$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule maxbin2_drep:
    input:
        f"{WORKDIR}/1_single_coverage/maxbin2/{{sample}}.tsv"
    output:
        f"{WORKDIR}/1_single_coverage/maxbin2_drep/{{sample}}/dereplicated_genomes.csv"
    params:
        outdir=f"{WORKDIR}/1_single_coverage/metabat2_drep/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 1000) * 2 ** (attempt - 1)))
    message: "Dereplicating MetaBAT2 bins for {wildcards.sample} at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3 checkm-genome/1.2.3
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        dRep dereplicate {params.outdir} -g {input} -p {threads} -pa 0.95
        """














rule semibin2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{sample}}.fna",
        bam=f"{WORKDIR}/1_single_coverage/bowtie2/{{sample}}.bam"
    output:
        f"{WORKDIR}/1_single_coverage/semibin2/{{sample}}/contig_bins.tsv"
    params:
        outdir=f"{WORKDIR}/1_single_coverage/semibin2/{{sample}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: min(1000*1024,max(8*1024, int(input.size_mb * 30) * 2 ** (attempt - 1))),
        runtime=lambda wildcards, input, attempt: min(20000,max(15, int(input.size_mb / 2) * 2 ** (attempt - 1)))
    message: "Binning contigs from assembly {wildcards.sample} using semibin2..."
    shell:
        """
        module load {params.semibin2_module} {params.bedtools_module} {params.hmmer_module}
        SemiBin2 single_easy_bin -i {input.assembly} -b {input.bam} -o {params.outdir} -m 1500 -t {threads} --compression none
        """

rule semibin2_table:
    input:
         f"{WORKDIR}/1_single_coverage/semibin2/{{sample}}/contig_bins.tsv"
    output:
        f"{WORKDIR}/1_single_coverage/semibin2/{{sample}}.tsv"
    params:
        fastadir=f"{WORKDIR}/1_single_coverage/semibin2/{{sample}}/output_bins"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(1*1024, int(input.size_mb * 10) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(5, int(input.size_mb / 5) * 2 ** (attempt - 1))
    shell:
        """
        tail -n +2 {input} > {output}
        """