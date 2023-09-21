rule verifybamid:
    input:
        unpack(lambda wildcards: {
            **get_bam(wildcards),
            "vcf": config["ref"]["1kg-variants"]
        })
    output:
        self_sm="output/{sample}/{sample}.selfSM"
    params:
        extra="--ignoreRG --verbose",
        sample_name="{sample}"
    log: "logs/{sample}/{sample}.verifybamid.log"
    resources:
        runtime=48 * 60,
        mem_mb=16000
    wrapper: "file:snakemake-wrapper/verifybamid"