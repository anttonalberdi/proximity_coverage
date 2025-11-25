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
        expand(f"{WORKDIR}/2_all_coverage/metabat2_drep/{{assembly}}/data_tables/genomeInformation.csv", assembly=ASSEMBLIES),
        expand(f"{WORKDIR}/2_all_coverage/maxbin2_drep/{{assembly}}/data_tables/genomeInformation.csv", assembly=ASSEMBLIES)

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
        f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_metabat.depth"
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
        jgi_summarize_bam_contig_depths --outputDepth {output} {input}
        """

rule maxbin_depth:
    input:
        f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_{{reads}}_maxbin.depth"
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
        depth=f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_metabat.depth"
    output:
        f"{WORKDIR}/2_all_coverage/metabat2/{{assembly}}.tsv"
    params:
        basename=f"{WORKDIR}/2_all_coverage/metabat2/{{assembly}}/{{assembly}}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 5) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using metabat2..."
    shell:
        """
        module load metabat2/2.17
        metabat2 -i {input.assembly} -a {input.depth} -o {params.basename} -m 1500 --saveCls

        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "$(basename {params.basename}).*.fa" | sort > {output}
        """

rule metabat2_checkm:
    input:
        f"{WORKDIR}/2_all_coverage/metabat2/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/2_all_coverage/metabat2_checkm/{{assembly}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/2_all_coverage/metabat2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/2_all_coverage/metabat2_checkm/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of MetaBAT2 bins for {wildcards.assembly} using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fa -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "genome","completeness","contamination"; next}} {{print $1".fa",$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule metabat2_drep:
    input:
        genomes=f"{WORKDIR}/2_all_coverage/metabat2/{{assembly}}.tsv",
        genomeinfo=f"{WORKDIR}/2_all_coverage/metabat2_checkm/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/2_all_coverage/metabat2_drep/{{assembly}}/data_tables/genomeInformation.csv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/2_all_coverage/metabat2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/2_all_coverage/metabat2_drep/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Dereplicating MetaBAT2 bins for {wildcards.assembly} at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -g {input.genomes} -p {threads} -pa 0.95 --genomeInfo {input.genomeinfo}
        """


rule maxbin2:
    input:
        assembly=f"{WORKDIR}/0_preprocessing/megahit/{{assembly}}.fna",
        depth = lambda wildcards: expand(
            f"{WORKDIR}/2_all_coverage/bowtie2/{{assembly}}_{{reads}}_maxbin.depth",
            assembly=wildcards.assembly,
            reads=READS
        )
    output:
        f"{WORKDIR}/2_all_coverage/maxbin2/{{assembly}}.tsv"
    params:
        basedir=f"{WORKDIR}/2_all_coverage/maxbin2/{{assembly}}",
        basename=f"{WORKDIR}/2_all_coverage/maxbin2/{{assembly}}/{{assembly}}",
        abund=lambda wildcards, input: " ".join(f"-abund {f}" for f in input.depth)
    threads: 1
    resources:
        mem_mb=lambda wildcards, input, attempt: max(8*1024, int(input.size_mb * 50) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: max(15, int(input.size_mb / 3) * 2 ** (attempt - 1))
    message: "Binning contigs from assembly {wildcards.assembly} using maxbin2..."
    shell:
        """
        module load maxbin2/2.2.7 hmmer/3.3.2
        rm -rf {params.basedir}
        mkdir -p {params.basedir}
        /opt/shared_software/shared_envmodules/conda/maxbin2-2.2.7/bin/run_MaxBin.pl -contig {input.assembly} {params.abund} -max_iteration 10 -out {params.basename} -min_contig_length 1500
        
        # Generate summary file for dRep
        find "$(dirname {params.basename})" -maxdepth 1 -type f -name "$(basename {params.basename}).*.fasta" | sort > {output}
        """

rule maxbin2_checkm:
    input:
        f"{WORKDIR}/2_all_coverage/maxbin2/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/2_all_coverage/maxbin2_checkm/{{assembly}}.tsv"
    params:
        bins_dir=lambda wildcards: f"{WORKDIR}/2_all_coverage/maxbin2/{wildcards.assembly}",
        outdir=f"{WORKDIR}/2_all_coverage/maxbin2_checkm/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(90, int(input.size_mb * 100) * 2 ** (attempt - 1)))
    message: "Assessing quality of Maxbin2 bins for {wildcards.assembly} using CheckM2..."
    shell:
        """
        module load checkm2/1.0.2
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        checkm2 predict -i {params.bins_dir}/*.fasta -o {params.outdir} -t {threads} --database_path /maps/datasets/globe_databases/checkm2/20250215/CheckM2_database/uniref100.KO.1.dmnd

        # Prepare genome info for drep
        awk -F'\t' 'BEGIN{{OFS=","}} NR==1{{print "genome","completeness","contamination"; next}} {{print $1".fasta",$2,$3}}' {params.outdir}/quality_report.tsv > {output}
        """

rule maxbin2_drep:
    input:
        genomes=f"{WORKDIR}/2_all_coverage/maxbin2/{{assembly}}.tsv",
        genomeinfo=f"{WORKDIR}/2_all_coverage/maxbin2_checkm/{{assembly}}.tsv"
    output:
        f"{WORKDIR}/2_all_coverage/maxbin2_drep/{{assembly}}/data_tables/genomeInformation.csv"
    params:
        outdir=f"{WORKDIR}/2_all_coverage/maxbin2_drep/{{assembly}}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, input, attempt: max(64*1024, int(input.size_mb * 1000) * 2 ** (attempt - 1)),
        runtime=lambda wildcards, input, attempt: min(20000, max(15, int(input.size_mb * 1000) * 2 ** (attempt - 1)))
    message: "Dereplicating MetaBAT2 bins for {wildcards.assembly} at 95% ANI..."
    shell:
        """
        module load drep/3.6.2 fastani/1.33 mash/2.3
        rm -rf {params.outdir}
        dRep dereplicate {params.outdir} -g {input.genomes} -p {threads} -pa 0.95 --genomeInfo {input.genomeinfo}
        """
