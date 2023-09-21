do_merge, cohort_vcf = get_cohort_vcf()
if do_merge:
    rule merge_vcfs:
        input: get_variants
        output: 
            vcf=cohort_vcf
        log: "logs/merge_vcfs/{project}.log"
        threads: 2
        resources:
            runtime=720,
            mem_mb=64000
        wrapper: "file:snakemake-wrapper/bcftools/merge"
else:
    rule convert_cohort_vcf:
        input: config.get("cohort_vcf", None)
        output:
            vcf=cohort_vcf
        log: "logs/convert_cohort_vcf/{project}.log"
        threads: 2
        resources:
            runtime=120,
            mem_mb=8192
        wrapper: "file:snakemake-wrapper/bcftools/view"
