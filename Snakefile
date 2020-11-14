import pandas as pd

##### parameters #####
SAMPLES = ['Id1_AA', 'Id2_AA', 'Id3_control', 'Id4_control']
UNITS = ['rep1', 'rep2', 'rep1', 'rep2']
DESIGN="~ condition"
ORGANISM="human"

samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)

##### target rule #####
rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.csv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast="AA_vs_control"),
        #"results/pca.svg",
        #"qc/multiqc_report.html"
        "counts/all.tsv"

##### rules #####
rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    threads: 4
    conda:
        "envs/align.yaml"
    params:
        length=49,
    shell:
        """
        STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases 11 \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.length}
        """
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

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
        gtf="genome/human.GRCh38.chr22.gtf",
        genome="genome/STARINDEX/Genome"
    output:
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
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
            --runThreadN {threads} \
            --runMode alignReads \
            --quantMode GeneCounts \
            --sjdbGTFfile {input.gtf} \
            --genomeDir {params.indexdir} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --outReadsUnmapped Fastq \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix star/{wildcards.sample}-{wildcards.unit}/
        """

rule count_matrix:
    input:
        expand("star/{sample}-{unit}/ReadsPerGene.out.tab", zip,
            sample=SAMPLES,
            unit=UNITS)
    output:
        "counts/all.tsv"
    params:
        samples=SAMPLES,
        strand="reverse",
        column=4
    log:
        "logs/counts/count_matrix.log"
    conda:
       "envs/pandas.yaml"
    script:
        "scripts/count-matrix.py"

rule deseq2:
    input:
        counts="counts/all.tsv",
        samples="samples.txt"
    output:
        table=report("results/diffexp/{contrast}.diffexp.csv",
                     "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg",
                       "../report/ma.rst"),
        up="results/diffexp/deg-sig-up.{contrast}.csv",
        down="results/diffexp/deg-sig-down.{contrast}.csv"
    params:
        contrast=['AA', 'control'],
        annotationhub=ORGANISM,
        design="~condition",
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/deseq2.R"

rule multiqc:
    input:
        expand("star/{samples}-{units}/Aligned.sortedByCoord.out.bam",
            samples=SAMPLES,
            units=UNITS),
        expand("qc/rseqc/{samples}-{units}.infer_experiment.txt",
            samples=SAMPLES,
            units=UNITS),
        expand("qc/rseqc/{samples}-{units}.stats.txt",
            samples=SAMPLES,
            units=UNITS),
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
            qc/rseqc star > {log}
        """
