# Alignment without priorizing the annotated junctions of the gff file
rule minimap_align:
    input:
        fa=config["ref_fa"],
        fq=config["fastq_dir"] + "{sample}.fastq"
    output:
        bam="aligned_tosep/{sample}.bam",
        bai="aligned_tosep/{sample}.bam.bai"
    threads:
        10
    params:
        thread_index= lambda w, threads : threads - 1
    shell:
        "minimap2 -ax splice -t {threads} {input.fa} {input.fq} | "
        "samtools sort -@ {threads} -o {output.bam} && "
        "samtools index -@ {params.thread_index} {output.bam}"
