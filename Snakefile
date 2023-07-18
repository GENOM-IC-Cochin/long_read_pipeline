configfile: "config.yaml"

import os
import sys

import pandas as pd

SAMPLES = [fastq_name.split(".")[0] for fastq_name in os.listdir(config["fastq_dir"]) if "fastq" in fastq_name or "fq" in fastq_name]
ORGANISMS = ["pNLGV", "human"]

samplesheet = pd.read_csv(config["samplesheet"], sep = "\t")

if not set(samplesheet["sampleID"].to_list()) == set(SAMPLES):
    sys.exit("""Mismatch between the samplesheet and the fastq names :
                the sampleID column must consist of the names of the fastq files.""")

with open("bam_list.txt", "w") as fileout:
    for _, row in samplesheet.iterrows():
        fileout.write(f"{config['bam_dir']}{row['sampleID']}_human_s.bam:{row['sampleID']}\n")


qc_rep = "qc/rnaseqc_report/multiqc_report.html"
counts = config["quant_dir"] + "ALL_SAMPLES/ALL_SAMPLES.transcript_model_grouped_counts.tsv"
diff_res = "analysis/adjusted_gene_transcript_pval_005.csv"
pca_plot = "analysis/pca_plot.png"
diff_plot = "analysis/signif_genes_and_isoforms.png"
nofilter_files = expand(config["bam_dir"] + "{sample}.bam", sample=SAMPLES)
filter_files = expand(config["bam_dir"] + "{sample}_{organism}_s.bam", sample=SAMPLES, organism=ORGANISMS)

rule_all_input_list = [qc_rep, counts, diff_res, pca_plot, diff_plot]


if config["filter_path"] == "yes":
    rule_all_input_list.extend(filter_files)
    include: "filter_rules.smk"
else:
    rule_all_input_list.extend(nofilter_files)
    include: "no_filter_rules.smk"

rule all:
    input:
        rule_all_input_list


# does not work : Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.
# I dont know what those parameters are for -x splice
# rule minimap_index:
#     input:
#         config["ref_fa"]
#     output:
#         config["ref_mmi"]
#     shell:
#         "minimap2 -d -x splice {output} {input}"

# Maybe add annotated junctions in bed format, to prioritize annotated splice junction.
# paftools.js gff2bed anno.gff > anno.bed # to create the bed file.

rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq=config["fastq_dir"] + "{sample}.fastq.gz"
    output:
        bam="aligned_tosep/{sample}.bam",
        bai="aligned_tosep/{sample}.bam.bai"
    threads:
        10
    params:
        thread_index=9
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} | "
        "samtools sort -@ {threads} -o {output.bam} && "
        "samtools index -@ {params.thread_index} {output.bam}"


rule rnaseqc:
    input:
        gtf=config["collaps_gtf"],
        bam=lambda wildcards : config["bam_dir"] + wildcards.sample + "_human_s.bam" if config["filter_path"] == "yes" else config["bam_dir"] + wildcards.sample + ".bam"
    output:
        "qc/rnaseqc/{sample}.metrics.tsv"
    params:
        dir="qc/rnaseqc",
        sample_name=lambda wildcards: wildcards.sample
    shell:
        "rnaseqc -s {params.sample_name} {input.gtf} {input.bam} {params.dir}"


rule rnaseqc_multiqc:
    input:
        expand("qc/rnaseqc/{sample}.metrics.tsv", sample=SAMPLES)
    output:
        "qc/rnaseqc_report/multiqc_report.html"
    params:
        in_dir="qc/rnaseqc/",
        out_dir="qc/rnaseqc_report"
    shell:
        "multiqc -f {params.in_dir} -o {params.out_dir}"


rule collapse_annotations:
    # python script from https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
    input:
        config["gtf"]
    output:
        config["collaps_gtf"]
    shell:
        "python3 collapse_annotation.py {input} {output}"


rule isoquant:
    input:
        bam=expand(config["bam_dir"] + "{sample}_human_s.bam" if config["filter_path"] == "yes" else config["bam_dir"] + "{sample}.bam", sample=SAMPLES),
        bam_list="bam_list.txt",
        gtf=config["gtf"],
        fa=config["ref_fa"]
    threads:
        16
    params:
        output_dir=config["quant_dir"]
    output:
        quants=config["quant_dir"] + "ALL_SAMPLES/ALL_SAMPLES.transcript_model_grouped_counts.tsv",
        gtf=config["quant_dir"] + "ALL_SAMPLES/ALL_SAMPLES.transcript_models.gtf"
    shell:
        "isoquant.py --data_type ont --reference {input.fa} "
        "--genedb {input.gtf} --bam_list {input.bam_list} --force "
        "--output {params.output_dir} -t {threads} --complete_genedb "
        "--prefix ALL_SAMPLES"


rule analysis_script:
    input:
        gtf=config["quant_dir"] + "ALL_SAMPLES/ALL_SAMPLES.transcript_models.gtf",
        quants=config["quant_dir"] + "ALL_SAMPLES/ALL_SAMPLES.transcript_model_grouped_counts.tsv"
    output:
        "analysis/adjusted_gene_transcript_pval_005.csv",
        "analysis/signif_genes_and_isoforms.png",
        "analysis/pca_plot.png"
    params:
        design=config["samplesheet"],
        output_dir="analysis"
    script:
        "analysis_script.R"
