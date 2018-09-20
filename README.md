Pipeline for differential gene expression analysis using long reads
==================================================================

This pipeline uses [snakemake](https://snakemake.readthedocs.io/en/stable/), [minimap2](https://github.com/lh3/minimap2) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to automate simple [differential transcript expression](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene) workflows on long read data.

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
-  TODO

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

Reference
--------

This pipeline is largely based on the approach desscribed in the following paper:

- Love MI, Soneson C and Patro R. *Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification.* F1000Research 2018, 7:952
(doi: [10.12688/f1000research.15398.2](http://dx.doi.org/10.12688/f1000research.15398.2))

