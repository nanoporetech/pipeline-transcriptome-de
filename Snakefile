
import os
from os import path

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
        min_mq = config["minimum_mapping_quality"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

rule all:
    input:
        mapped_samples = expand("alignments/{sample}.bam", sample=all_samples.keys())
