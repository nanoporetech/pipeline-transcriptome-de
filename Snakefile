
import os
from os import path
import pandas as pd
from collections import OrderedDict

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
RESDIR =  config["resdir"]
SNAKEDIR = path.dirname(workflow.snakefile)
PY2_EXEC = "python2 {}/scripts".format(SNAKEDIR)

include: "snakelib/utils.snake"

control_samples = config["control_samples"]
treated_samples = config["treated_samples"]

all_samples = config["control_samples"].copy()
all_samples.update(config["treated_samples"])
datasets = [path.basename(x).rsplit(".", 1)[0] for x in all_samples.values()]

rule build_minimap_index: ## build minimap2 index
    input:
        genome = config["transcriptome"]
    output:
        index = "index/transcriptome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """

rule map_reads: ## map reads using minimap2
    input:
       index = rules.build_minimap_index.output.index,
       fastq = lambda wildcards: all_samples[wildcards.sample]
    output:
       bam = "alignments/{sample}.bam"
    params:
        opts = config["minimap2_opts"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

rule count_reads:
    input:
        bam = rules.map_reads.output.bam
    output:
        tsv = "counts/{sample}.tsv"
    params:
        min_mq = config["minimum_mapping_quality"],
    conda: "env.yml"
    shell: """
        {SNAKEDIR}/scripts/bam_count_reads.py -a {params.min_mq} -t {output.tsv} {input.bam}
    """

rule merge_counts:
    input:
        count_tsvs = expand("counts/{sample}.tsv", sample=all_samples.keys()),
    output:
        tsv = "merged/all_counts.tsv"
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/scripts/merge_count_tsvs.py -z -o {output.tsv} {input.count_tsvs}
    """

rule write_coldata:
    input:
    output:
        coldata = "de_analysis/coldata.tsv"
    run:
        samples, conditions, types = [], [], []
        for sample in control_samples.keys():
            samples.append(sample)
            conditions.append("untreated")
            types.append("single-read")
        for sample in treated_samples.keys():
            samples.append(sample)
            conditions.append("treated")
            types.append("single-read")

        df = pd.DataFrame(OrderedDict([('sample', samples),('condition', conditions),('type', types)]))
        df.to_csv(output.coldata, sep="\t", index=False)

rule deseq_analysis:
    input:
        coldata = rules.write_coldata.output.coldata,
        tsv = rules.merge_counts.output.tsv,
    output:
        res = "de_analysis/deseq2_results.tsv",
        pdf = "de_analysis/deseq2_plots.pdf",
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/scripts/deseq2_analysis.R
    """


rule all:
    input:
        mapped_samples = expand("alignments/{sample}.bam", sample=all_samples.keys()),
        count_tsvs = expand("counts/{sample}.tsv", sample=all_samples.keys()),
        merged_tsv = "merged/all_counts.tsv",
        coldata = "de_analysis/coldata.tsv",
        res = "de_analysis/deseq2_results.tsv",
        
