rule de_novo_cnv:
    input: unpack(get_dnm_cnv_input)
    output: 
        vcf="output/{family}/{family}.cnv.vcf.gz",
        tbi="output/{family}/{family}.cnv.vcf.gz.tbi",
        md5="output/{family}/{family}.cnv.vcf.gz.md5sum",
    params:
        license=license_text,
        reference_dir=config['reference_dir'],
        sample="{family}"
    log: "logs/{family}/{family}.{family}.de_novo_cnv.log"
    group: "dnm_cnv"
    wrapper: "file:wrapper/dragen/deNovoCnvCalling"

rule de_novo_to_remote:
    input:
        vcf=rules.de_novo_cnv.output.vcf,
        tbi=rules.de_novo_cnv.output.tbi,
        md5=rules.de_novo_cnv.output.md5,
    output:
        vcf=S3.remote(config['s3']['output'] + '/{family}/{family}.cnv.vcf.gz'),
    params:
        sample="{family}"
    log: "logs/{family}/{family}.de_novo_to_remote.log"
    group: "dnm_cnv"
    wrapper: "file:wrapper/s3/upload"
