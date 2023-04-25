SAMPLES = [
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode08",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode09",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode10",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode11",
    "20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode12"
]


configfile: "config.yaml"

rule all:
    input:
        expand(config["storage_dir"] + "aligned/{sample}.bam", sample=SAMPLES),
        "qc/rnaseqc_report/multiqc_report.html",
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07/00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07.transcript_model_grouped_counts.tsv"


# Maybe add annotated junctions in bed format, to prioritize annotated splice junction.
# paftools.js gff2bed anno.gff > anno.bed # to create the bed file.

rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq=config["storage_dir"] + "{sample}.fastq"
    output:
        config["storage_dir"] + "aligned/{sample}.bam"
    threads:
        10
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} | samtools sort -@ {threads} -o {output} -"


rule index:
    input:
        config["storage_dir"] + "aligned/{sample}.bam"
    output:
        config["storage_dir"] + "aligned/{sample}.bam.bai"
    threads:
        4
    params:
        add_threads=3
    shell:
        "samtools index -@ {params.add_threads} {input}"


rule rnaseqc:
    input:
        gtf=config["collaps_gtf"],
        bam= config["storage_dir"] + "aligned/{sample}.bam"
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
        bam=expand(config["storage_dir"] + "aligned/{sample}.bam", sample=SAMPLES),
        bam_list="bam_list.txt",
        bam_index=expand(config["storage_dir"] + "aligned/{sample}.bam.bai", sample=SAMPLES),
        gtf=config["gtf"],
        fa=config["ref_fa"]
    threads: 16
    params:
        output_dir=config["quant_dir"]
    output:
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07/00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07.transcript_model_grouped_counts.tsv"
    shell:
        "isoquant.py --data_type ont --reference {input.fa} "
        "--genedb {input.gtf} --bam_list {input.bam_list} --force "
        "--output {params.output_dir} -t {threads} --complete_genedb "
        "--read_group file_name"



# Maybe add cramino at some point, but would be redundant with ToulligQC?
