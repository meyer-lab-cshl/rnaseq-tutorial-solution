import pandas as pd

##### functions #####
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

# parameters
samplesfile = "samples.txt"
STRAND = "reverse"
DESIGN="~ condition"
SPECIES="human"
SAMPLESFILE="samples.txt"

##### load sample sheets #####
samples = pd.read_table(SAMPLESFILE).set_index(["sample", "unit"], drop=False)

##### setup report #####
report: "report/workflow.rst"

##### target rule #####
rule all:
    input:
        expand("results/diffexp/{contrast}.diffexp.txt",
            contrast = "AA_vs_control"),
        "results/diffexp/pca.pdf",
        "qc/multiqc_report.html"

##### rules #####
rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    threads: 1
    conda:
        "envs/align.yaml"
    params:
        length=75,
        Nbases=11
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases {params.Nbases} \
            --sjdbOverhang {params.length}
        """

rule cutadapt:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters="CTGACCTCAAGTCTGCACACGAGAAGGCTAG",
    threads: 1
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt \
            -a {params.adapters} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j {threads} \
            {input} \
        > {output.qc}
        """

rule align:
    input:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        genome="genome/STARINDEX/Genome"
    output:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
        "star/{sample}-{unit}.ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        indexdir="genome/STARINDEX"
    threads: 4
    conda:
        "envs/align.yaml"
    shell:
        """
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.indexdir} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --outFileNamePrefix star/{wildcards.sample}-{wildcards.unit}. \
            --quantMode GeneCounts \
            --outSAMtype BAM SortedByCoordinate
        """

rule index:
    input:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
    conda:
        "envs/index.yaml"
    shell:
        "samtools index {input}"

rule rseqc_infer:
    input:
        bed="genome/human.GRCh38.chr22.bed",
        bam="star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "qc/rseqc/{sample}-{unit}.infer_experiment.txt"
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.log"
    conda:
        "envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

rule rseqc_coverage:
    input:
        bed="genome/human.GRCh38.chr22.bed",
        bam="star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "qc/rseqc/{sample}-{unit}.geneBodyCoverage.txt"
    log:
        "logs/rseqc/rseqc_coverage/{sample}-{unit}.log"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o qc/rseqc/{wildcards.sample}-{wildcards.unit} 2> {log}
        """

rule multiqc:
    input:
        expand("qc/rseqc/{samples.sample}-{samples.unit}.geneBodyCoverage.txt",
            samples=samples.itertuples()),
        expand("qc/rseqc/{samples.sample}-{samples.unit}.infer_experiment.txt",
            samples=samples.itertuples())
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir qc \
            --filename multiqc_report.html \
            trimmed star qc/rseqc > {log}
        """

rule count_matrix:
    input:
        expand("star/{samples.sample}-{samples.unit}.ReadsPerGene.out.tab",
            samples=samples.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=samples['sample'].tolist(),
        strand="reverse"
    log:
        "logs/counts/count_matrix.log"
    conda:
       "envs/pandas.yaml"
    script:
        "scripts/count-matrix.py"

rule setup_de:
    input:
        counts="counts/all.tsv",
        samples=SAMPLESFILE
    output:
        dds="deseq2/all.rds"
    params:
        species=SPECIES,
        design=DESIGN
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/setup.log"
    script:
        "scripts/setup_deseq2.R"

rule deseq2:
    input:
        dds="deseq2/all.rds",
    output:
        table=report("results/diffexp/{contrast}.diffexp.txt",
            "report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.pdf",
            "report/ma.rst"),
        up="results/diffexp/deg-sig-up_{contrast}.csv",
        down="results/diffexp/deg-sig-down_{contrast}.csv"
    params:
        design=DESIGN,
        samples=SAMPLESFILE
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/deseq2.R"

rule pca:
    input:
        dds="deseq2/all.rds"
    output:
        pca_plot=report("results/diffexp/pca.pdf", "report/pca.rst"),
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/pca.log"
    script:
        "scripts/pca.R"
