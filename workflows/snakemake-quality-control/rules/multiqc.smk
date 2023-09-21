rule multiqc:
    input: get_multiqc_data
    output: report("output/{project}/{project}.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        extra="--interactive"
    log: "logs/{project}/{project}.log"
    resources:
        runtime=30,
        mem_mb=8000
    wrapper: "file:snakemake-wrapper/multiqc"