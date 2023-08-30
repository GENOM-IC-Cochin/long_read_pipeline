# Alignment priorizing the annotated junctions of the gtf file

rule create_bed_align:
    input:
        gtf=config["gtf"]
    output:
        bed="junc.bed"
    shell:
        "paftools.js gff2bed {input.gtf} > {output.bed}"

rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq=config["fastq_dir"] + "{sample}.fastq",
        bed="junc.bed"
    output:
        bam="aligned_tosep/{sample}.bam",
        bai="aligned_tosep/{sample}.bam.bai"
    threads:
        10
    params:
        thread_index= lambda w, threads : threads - 1
    shell:
        "minimap2 -ax splice --junc-bed {input.bed} -t {threads} {input.fa} {input.fq} | "
        "samtools sort -@ {threads} -o {output.bam} && "
        "samtools index -@ {params.thread_index} {output.bam}"
