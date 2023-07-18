rule samtools_filter:
    input:
        bam="aligned_tosep/{sample}.bam",
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
        index=config["bam_dir"] + "{sample}_{organism}_s.bam.bai",
        sortd=config["bam_dir"] + "{sample}_{organism}_s.bam"
    threads:
        4
    params:
        thread_index=3
    shell:
        "samtools sort -@ {threads} {input} -o {output.sortd} && "
        "samtools index -@ {params.thread_index} {output.sortd}"
