rule peddy:
    input:
        vcf=get_cohort_vcf()[1],
        ped=config.get("pedigree", "")
    output:
        ped=temp("output/{project}/{project}.peddy.ped"),
        het_check=temp("output/{project}/{project}.het_check.csv"),
        ped_check=temp("output/{project}/{project}.ped_check.csv"),
        sex_check=temp("output/{project}/{project}.sex_check.csv"),
        background_pca=temp("output/{project}/{project}.background_pca.json")
    params:
        prefix="{project}"
    log: "logs/peddy/{project}.log"
    threads: 8
    resources:
        runtime=480,
        mem_mb=32000
    wrapper: "file:snakemake-wrapper/peddy"
