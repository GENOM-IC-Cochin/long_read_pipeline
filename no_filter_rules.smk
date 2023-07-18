rule mv_bams:
    input:
        bam="aligned_tosep/{sample}.bam",
        bai="aligned_tosep/{sample}.bam.bai"
    output:
        config["bam_dir"] + "{sample}.bam",
        config["bam_dir"] + "{sample}.bam.bai"
    params:
        out_dir=config["bam_dir"]
    shell:
        "mv {input.bam} {input.bai} {params.out_dir}"
