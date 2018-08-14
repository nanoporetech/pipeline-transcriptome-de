Pipeline for differential gene expression analysis using long reads
==================================================================

This pipeline uses [minimap2](https://github.com/lh3/minimap2) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to automate simple [differential transcript expression](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene) workflows on long read data.

Dependencies 
------------
- [miniconda](https://conda.io/miniconda.html) - install it according to the [instructions](https://conda.io/docs/user-guide/install/index.html).
- [snakemake](https://anaconda.org/bioconda/snakemake) install using `conda`.
- [pandas](https://anaconda.org/conda-forge/pandas) - install using `conda`.
- The rest of the dependencies are automatically installed using the `conda` feature of `snakemake`.

Installation
------------

Clone the repository:

```bash
git clone XXX
```

Input
-----

The input files and parameters are specified in `config.yml`:

- `transcriptome` - the input transcriptome.
- `control_samples` - a dictionary with control sample names and paths to the fastq files.
- `treated_samples` - a dictionary with treated sample names and paths to the fastq files.

Usage
-----

Edit `config.yml` to set the input datasets and parameters then issue:

```bash
snakemake --use-conda -j <num_cores> all
```

and issue snakemake <target> to invoke the target of your choice. 

Output
-----

- `alignments/*.bam` - sorted and indexed transcriptome alignments.
- `merged/all_counts.tsv` - the transcript count table including all samples.
- `de_analysis/coldata.tsv` - the condition table supplied to DESeq2.
- `de_analysis/deseq2_results.tsv` - the results of the DESeq2 analysis. This table contains information on which genes are up-regulated in the experimental group (treated) with reference to the control group (untreated). The table contains:
    * `row names` = transcript ID
    * `baseMean` = the base mean value for each gene
    * `log2FoldChange` = The log2 fold change of each transcript in the experiment group compared with the control (up- or down-regulated)
    * `lfcSE` = log 2-fold change standard error
    * `stat` = The Wald test statistic (Log2foldchange/lfcSE)
    * `pvalue` = non-corrected pvalue (stat value compared against a normal distribution)
    * `padj` = [Bonferroni corrected](https://en.wikipedia.org/wiki/Bonferroni_correction) pvalue (acounting for multiple comparisons)
- `de_analysis/deseq2_plots.pdf` - visualisations of the DE analysis results and diagnostic plots.

Layout
------

* `README.md`
* `Snakefile`         - master snakefile
* `config.yml`        - YAML configuration file
* `snakelib/`         - snakefiles collection included by the master snakefile
* `lib/`              - python files included by analysis scripts and snakefiles
* `scripts/`          - analysis scripts
* `data/`             - input data needed by pipeline - use with caution to avoid bloated repo
* `results/`          - pipeline results to be commited - use with caution to avoid bloated repo
* `requirements.txt`  - list of python package dependencies

Useful snakemake targets
------------------------

```
help                    list all targets and descriptions
info                    print pipeline information
clean_workdir           delete working directory. WARNING: all data will be lost!
clean_resdir            delete results directory. WARNING: all data will be lost!
```
