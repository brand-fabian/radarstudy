rule joint_genotyping:
    input: unpack(gather_samples)
    output:
        hard_filtered="output/{family}/{family}.hard-filtered.vcf.gz",
        hard_filtered_idx="output/{family}/{family}.hard-filtered.vcf.gz.tbi",
        vcf="output/{family}/{family}.vcf.gz",
        vcf_idx="output/{family}/{family}.vcf.gz.tbi"
    group: "dragen_joint"
    log: "logs/{family}/joint_genotyping.log"
    params:
        reference_dir=config['reference_dir'],
        family="{family}",
        license=license_text
    wrapper: "file:wrapper/dragen/jointGenotyping"

rule copy_vcf_to_remote:
    input: 
        hard_filtered=rules.joint_genotyping.output.hard_filtered
    output: S3.remote(config['s3']['output'] + "/{family}/{family}.hard-filtered.vcf.gz")
    group: "dragen_joint"
    log: "logs/{family}/copy_to_remote.log"
    params:
        sample="{family}"
    wrapper: "file:wrapper/s3/upload"
