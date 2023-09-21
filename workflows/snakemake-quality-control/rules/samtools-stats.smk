rule samtools_stats:
    input: unpack(get_bams)
    output:
        stats="output/{sample}/{sample}.samtools-stats.txt",
        flagstat="output/{sample}/{sample}.flagstat.txt",
        idxstats="output/{sample}/{sample}.idxstats.txt"
    log: "logs/{sample}/{sample}.samtools-stats.log"
    resources:
        runtime=1440,
        mem_mb=8000
    wrapper: "file:snakemake-wrapper/samtools/stats"