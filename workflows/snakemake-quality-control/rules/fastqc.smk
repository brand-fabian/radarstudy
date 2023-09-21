rule fastqc:
    input: get_fastq
    output: 
        html=[
            temp("output/{sample}/{sample}_1.html"),
            temp("output/{sample}/{sample}_2.html")
        ],
        zip=[
            temp("output/{sample}/{sample}_1_fastqc.zip"),
            temp("output/{sample}/{sample}_2_fastqc.zip")
        ]
    params:
        extra=""
    log: "logs/{sample}/{sample}.fastqc.log"
    resources:
        runtime=720,
        mem_mb=32000
    wrapper: "file:snakemake-wrapper/fastqc"