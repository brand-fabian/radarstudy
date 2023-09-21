rule ngs_chew_fingerprint:
    input:
        unpack(lambda wildcards: {
            **get_bam(wildcards),
            "ref": config["ref"]["genome"]
        })
    output:
        vcf=temp("output/{sample}/{sample}.fingerprint.vcf.gz"),
        fingerprint=temp("output/{sample}/{sample}.npz")
    log: "logs/{sample}/{sample}.ngs-chew-fingerprint.log"
    resources:
        runtime=300,
        mem_mb=32000
    wrapper: "file:snakemake-wrapper/ngs-chew/fingerprint"

rule ngs_chew_compare:
    input: unpack(get_ngs_chew_fingerprints)
    output:
        aab=report("output/{project}/{project}.aab.html", caption="../report/ngs-chew-aab.rst", category="ngs-chew"),
        var_het=report("output/{project}/{project}.var_het.html", caption="../report/ngs-chew-var_het.rst", category="ngs-chew"),
        relatedness=report("output/{project}/{project}.relatedness.html", caption="../report/ngs-chew-relatedness.rst", category="ngs-chew")
    log: "logs/{project}/{project}.ngs-chew-compare.log"
    resources:
        runtime=60,
        mem_mb=32000
    wrapper: "file:snakemake-wrapper/ngs-chew/compare"
