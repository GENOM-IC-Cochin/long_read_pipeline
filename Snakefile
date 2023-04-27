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
        expand(config["storage_dir"] + "aligned/{sample}_{organism}_s.bam", sample=SAMPLES, organism=ORGANISMS),
        "qc/rnaseqc_report/multiqc_report.html",
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_human_s/00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_human_s.transcript_model_grouped_counts.tsv"


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
        fq=config["storage_dir"] + "{sample}.fastq"
    output:
        bam="aligned_sep/{sample}.bam",
        bai="aligned_sep/{sample}.bam.bai"
    threads:
        10
    params:
        thread_index=9
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} | "
        "samtools sort -@ {threads} -o {output.bam} && "
        "samtools index -@ {params.thread_index} {output.bam}"


rule samtools_filter:
    input:
        bam="aligned_sep/{sample}.bam",
        bed_filter="viral_chr.bed"
    output:
        human=temp("aligned_sep/{sample}_human.bam"),
        pNLGV=temp("aligned_sep/{sample}_pNLGV.bam")
    threads:
        4
    params:
        thread_supp=3
    shell:
        "samtools view -hb -@ {params.thread_supp} -U {output.human} -L {input.bed_filter} {input.bam} > {output.pNLGV}"


rule sort_index:
    input:
        "aligned_sep/{sample}_{organism}.bam"
    output:
        index=config["storage_dir"] + "aligned/{sample}_{organism}_s.bam.bai",
        sortd=config["storage_dir"] + "aligned/{sample}_{organism}_s.bam"
    threads:
        4
    params:
        thread_index=3
    shell:
        "samtools sort -@ {threads} {input} -o {output.sortd} && "
        "samtools index -@ {params.thread_index} {output.sortd}"


rule rnaseqc:
    input:
        gtf=config["collaps_gtf"],
        bam= config["storage_dir"] + "aligned/{sample}_human_s.bam"
    output:
        "qc/rnaseqc/{sample}_human.metrics.tsv"
    params:
        dir="qc/rnaseqc",
        sample_name=lambda wildcards: wildcards.sample + "_human"
    shell:
        "rnaseqc -s {params.sample_name} {input.gtf} {input.bam} {params.dir}"


rule rnaseqc_multiqc:
    input:
        expand("qc/rnaseqc/{sample}_human.metrics.tsv", sample=SAMPLES)
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
        bam=expand(config["storage_dir"] + "aligned/{sample}_human_s.bam", sample=SAMPLES),
        bam_list="bam_list.txt",
        gtf=config["gtf"],
        fa=config["ref_fa"]
    threads:
        16
    params:
        output_dir=config["quant_dir"]
    output:
        config["quant_dir"] + "00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_human_s/00_20200818_NLGV_4GSTm_deltaTat_A2020_pass_barcode07_human_s.transcript_model_grouped_counts.tsv"
    shell:
        "isoquant.py --data_type ont --reference {input.fa} "
        "--genedb {input.gtf} --bam_list {input.bam_list} --force "
        "--output {params.output_dir} -t {threads} --complete_genedb "
        "--read_group file_name"



# Maybe add cramino at some point, but would be redundant with ToulligQC?
