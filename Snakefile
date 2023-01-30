SAMPLES = [
        "SRR8360560",
        "SRR8360562",
        "SRR8360563",
        "SRR8360564",
        "SRR8360565",
        "SRR8360566"
        ]

configfile: "config.yaml"

rule all:
    input:
        expand("aligned/{sample}_sorted_indexed.bam", sample=SAMPLES),
        "qc/rnaseqc_report/multiqc_report.html",
        expand("quants/00_{sample}_sorted_indexed/00_{sample}_sorted_indexed.transcript_tpm.tsv", sample=SAMPLES)


# Maybe add annotated junctions in bed format, to prioritize annotated splice junction.
# paftools.js gff2bed anno.gff > anno.bed # to create the bed file.

rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq="{sample}.fastq"
    output:
        temp("aligned/{sample}.sam")
    threads:
        10
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} > {output}"


rule sam2bam:
    input:
        "aligned/{sample}.sam"
    output:
        temp("aligned/{sample}.bam")
    shell:
        "samtools view -hbo {output} {input}"


rule sort_index:
    input:
        "aligned/{sample}.bam"
    output:
        "aligned/{sample}_sorted_indexed.bam"
    shell:
        "samtools sort {input} -o {output} && samtools index {output}"


rule rnaseqc:
    input:
        gtf=config["collaps_gtf"],
        bam="aligned/{sample}_sorted_indexed.bam"
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
        bam="aligned/{sample}_sorted_indexed.bam",
        gtf=config["gtf"],
        fa=config["ref_fa"]
    threads: 16
    params:
        output_dir="quants"
    output:
        "quants/00_{sample}_sorted_indexed/00_{sample}_sorted_indexed.transcript_tpm.tsv"
    shell:
        "isoquant.py --data_type ont --reference {input.fa} "
        "--genedb {input.gtf} --bam {input.bam} --force "
        "--output {params.output_dir} -t {threads} --complete_genedb "
        "--transcript_quantification unique_only"



# Maybe add cramino at some point, but would be redundant with ToulligQC?
