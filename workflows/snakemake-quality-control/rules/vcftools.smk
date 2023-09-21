rule relatedness2:
    input:
        vcf=get_cohort_vcf()[1]
    output:
        vcf="output/{project}/{project}.relatedness2"
    params:
        extra="--relatedness2"
    threads: 2
    resources:
        runtime=120,
        mem_mb=8192
    log: "logs/relatedness2/{project}.log"
    wrapper: "file:snakemake-wrapper/vcftools"