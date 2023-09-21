rule bcftools_stats:
    input: get_variant
    output: "output/{sample}/{sample}.bcftools-stats.txt"
    log: "logs/bcftools/{sample}.txt"
    threads: 2
    resources:
        runtime=60,
        mem_mb=16000
    wrapper: "file:snakemake-wrapper/bcftools/stats"

rule bcftools_stats_cohort:
    input: get_cohort_vcf()[1]
    output: "output/{project}/{project}.bcftools-stats.txt"
    log: "logs/bcftools/{project}.log"
    threads: 2
    resources:
        runtime=60,
        mem_mb=16000
    wrapper: "file:snakemake-wrapper/bcftools/stats"