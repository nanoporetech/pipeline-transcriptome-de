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
       bam = "alignments/{sample}.bam",
       sbam = "sorted_alignments/{sample}.bam",
    params:
        opts = config["minimap2_opts"],
        msec = config["maximum_secondary"],
        psec = config["secondary_score_ratio"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax map-ont -p {params.psec} -N {params.msec} {params.opts} {input.index} {input.fastq}\
    | samtools view -Sb > {output.bam};
    samtools sort -@ {threads} {output.bam} -o {output.sbam};
    samtools index {output.sbam};
    """

rule count_reads:
    input:
        bam = rules.map_reads.output.bam,
        trs = config["transcriptome"],
    output:
        tsv_dir = "counts/{sample}_salmon",
        tsv = "counts/{sample}_salmon/quant.sf",
    params:
    conda: "env.yml"
    threads: config["threads"]
    shell: """
        salmon quant --noErrorModel -p {threads} -t {input.trs} -l SF -a {input.bam} -o {output.tsv_dir}
    """

rule merge_counts:
    input:
        count_tsvs = expand("counts/{sample}_salmon/quant.sf", sample=all_samples.keys()),
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

rule write_de_params:
    input:
    output:
        de_params = "de_analysis/de_params.tsv"
    run:
        d = OrderedDict()
        d["Annotation"] = [config["annotation"]]
        d["min_samps_gene_expr"] = [config["min_samps_gene_expr"]]
        d["min_samps_feature_expr"] = [config["min_samps_feature_expr"]]
        d["min_gene_expr"] = [config["min_gene_expr"]]
        d["min_feature_expr"] = [config["min_feature_expr"]]
        df = pd.DataFrame(d)
        df.to_csv(output.de_params, sep="\t", index=False)


rule de_analysis:
    input:
        de_params = rules.write_de_params.output.de_params,
        coldata = rules.write_coldata.output.coldata,
        tsv = rules.merge_counts.output.tsv,
    output:
        res_dge = "de_analysis/results_dge.tsv",
        pdf_dge = "de_analysis/results_dge.pdf",
        res_dtu_gene = "de_analysis/results_dtu_gene.tsv",
        res_dtu_trs = "de_analysis/results_dtu_transcript.tsv",
        res_dtu_stager = "de_analysis/results_dtu_stageR.tsv",
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/scripts/de_analysis.R
    """


rule all:
    input:
        count_tsvs = expand("counts/{sample}_salmon", sample=all_samples.keys()),
        merged_tsv = "merged/all_counts.tsv",
        coldata = "de_analysis/coldata.tsv",
        de_params = "de_analysis/de_params.tsv",
        res_dge = "de_analysis/results_dge.pdf",
