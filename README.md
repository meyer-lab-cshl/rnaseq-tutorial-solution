# RNA-seq tutorial: solutions

This repository contains the fully developed, [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
RNAseq analysis pipeline built in this [course](https://github.com/meyer-lab-cshl/rnaseq-tutorial). If you haven't done so, please visit the [course repository](https://github.com/meyer-lab-cshl/rnaseq-tutorial) and follow
the instructions there, including the setup.

During the lectures, we built an analysis pipeline for RNAseq including alignment, quality control and differential expression
analysis. The [homework](HOMEWORK.md) assignment asked you to extend the pipeline, evaluate your results and
learn how to generate an analysis report. 

The files you see in this directory are:
* part of the pipeline we built: [Snakefile](Snakefile),[scripts](scripts) and [envs](envs)
* contain the example data to run the analysis: [genome](genome), [reads](reads) and [samples.txt](samples.txt))
* reflect the results of the analysis, including the Homework assignment: [results](results), [qc](qc), and [report](report.html)

