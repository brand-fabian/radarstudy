rule sv_calling:
    input: unpack(get_sv_input)
    output:
        vcf="output/{family}/{family}.sv.vcf.gz",
        tbi="output/{family}/{family}.sv.vcf.gz.tbi",
    params:
        license=license_text,
        reference_dir=config['reference_dir'],
        sample="{family}"
    log: "logs/{family}/{family}.{family}.sv_calling.log"
    group: "sv_calling"
    wrapper: "file:wrapper/dragen/svCalling"

rule sv_to_remote:
    input:
        vcf=rules.sv_calling.output.vcf,
        tbi=rules.sv_calling.output.tbi
    output:
        vcf=S3.remote(config['s3']['output'] + '/{family}/{family}.sv.vcf.gz')
    params:
        sample="{family}"
    log: "logs/{family}/{family}.sv_to_remote.log"
    group: "sv_calling"
    wrapper: "file:wrapper/s3/upload"