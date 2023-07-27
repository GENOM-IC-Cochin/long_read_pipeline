rule samtools_filter:
    input:
        bam="aligned_tosep/{sample}.bam",
        bed_filter=config["bed_file"]
    output:
        org1=temp(f'aligned_sep/{{sample}}_{ORGANISMS[0]}.bam'),
        org2=temp(f'aligned_sep/{{sample}}_{ORGANISMS[1]}.bam')
    threads:
        4
    params:
        thread_supp=3
    shell:
        "samtools view -hb -@ {params.thread_supp} -U {output.org1} -L {input.bed_filter} {input.bam} > {output.org2}"


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
