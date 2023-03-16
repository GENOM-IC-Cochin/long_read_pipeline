SAMPLES = [
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode08",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode09",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode10",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode11",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode12"
]

ORGANISMS = ["pNLGV", "human"]

configfile: "config.yaml"

rule all:
    input:
        expand(config["storage_dir"] + "aligned/{sample}_{organism}_sorted_indexed.bam", sample=SAMPLES, organism=ORGANISMS),
        "qc/rnaseqc_report/multiqc_report.html",
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_sorted_indexed00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_sorted_indexed.transcript_counts.tsv"


# Maybe add annotated junctions in bed format, to prioritize annotated splice junction.
# paftools.js gff2bed anno.gff > anno.bed # to create the bed file.

rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq=config["storage_dir"] + "{sample}.fastq"
    output:
        human=temp("aligned/{sample}_human.bam"),
        pNLGV=temp("aligned/{sample}_pNLGV.bam")

    threads:
        10
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} | "
        "samtools view -hb -U {output.human} pNLGV:1-15489 > {output.pNLGV}"


# rule sam2bam:
#     input:
#         "aligned/{sample}.sam"
#     output:
#         temp("aligned/{sample}.bam")
#     shell:
#         "samtools view -hbo {output} {input}"


rule sort_index:
    input:
        "aligned/{sample}_{organism}.bam"
    output:
        config["storage_dir"] + "aligned/{sample}_{organism}_sorted_indexed.bam"
    threads:
        4
    shell:
        "samtools sort -@ {threads} {input} -o {output} && samtools index {output}"


rule rnaseqc:
    input:
        gtf=config["collaps_gtf"],
        bam= config["storage_dir"] + "aligned/{sample}_human_sorted_indexed.bam"
    output:
        "qc/rnaseqc/{sample}.metrics.tsv"
    params:
        dir="qc/rnaseqc"
    shell:
        "rnaseqc -s {wildcards.sample} {input.gtf} {input.bam} {params.dir}"


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
        bam=expand(config["storage_dir"] + "aligned/{sample}_human_sorted_indexed.bam", sample=SAMPLES),
        bam_list="bam_list.txt",
        gtf=config["gtf"],
        fa=config["ref_fa"]
    threads: 16
    params:
        output_dir=config["quant_dir"]
    output:
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_sorted_indexed00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_sorted_indexed.transcript_counts.tsv"
    shell:
        "isoquant.py --data_type ont --reference {input.fa} "
        "--genedb {input.gtf} --bam_list {input.bam_list} --force "
        "--output {params.output_dir} -t {threads} --complete_genedb "
        "--read_group file_name"



# Maybe add cramino at some point, but would be redundant with ToulligQC?
