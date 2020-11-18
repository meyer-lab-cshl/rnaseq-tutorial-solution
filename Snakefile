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
        expand(["results/diffexp/{contrast}.diffexp.txt",
                "results/diffexp/{contrast}.ma-plot.pdf"],
               contrast="AA_vs_control"),
        "qc/multiqc_report.html",

##### setup report #####
report: "report/workflow.rst"


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
        length=75,
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases 11 \
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


rule count_matrix:
    input:
        expand("star/{sample}-{unit}.ReadsPerGene.out.tab", zip,
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
        table=report("results/diffexp/{contrast}.diffexp.txt", "report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.pdf", "report/ma.rst"),
        pca_plot=report("results/diffexp/{contrast}.pca-plot.pdf", "report/pca.rst"),
        up="results/diffexp/deg-sig-up_{contrast}.csv",
        down="results/diffexp/deg-sig-down_{contrast}.csv"
    params:
        contrast=['AA', 'control'],
        organism=ORGANISM,
        design="~condition",
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/deseq2.R"

rule multiqc:
    input:
        expand("star/{samples}-{units}.Aligned.sortedByCoord.out.bam", zip,
            samples=SAMPLES,
            units=UNITS),
        expand("qc/rseqc/{samples}-{units}.infer_experiment.txt", zip,
            samples=SAMPLES,
            units=UNITS),
        expand("qc/rseqc/{samples}-{units}.geneBodyCoverage.txt", zip,
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
            trimmed star qc/rseqc > {log}
        """
