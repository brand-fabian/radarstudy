rule cnv_calling:
    input: unpack(cnv_input)
    output:
        vcf="output/{sample}/{sample}.cnv.vcf.gz",
        tbi="output/{sample}/{sample}.cnv.vcf.gz.tbi",
        md5="output/{sample}/{sample}.cnv.vcf.gz.md5sum",
        tn="output/{sample}/{sample}.tn.tsv",
        calls="output/{sample}/{sample}.seg.called.merged",
        gff="output/{sample}/{sample}.cnv.gff3",
        bw="output/{sample}/{sample}.tn.bw",
        calls_png="output/{sample}/{sample}.cnv.calls.png",
        gc_png="output/{sample}/{sample}.cnv.gc.curve.png",
        cnv_png="output/{sample}/{sample}.cnv.png"
    params:
        extra="--cnv-enable-plots true --cnv-enable-tracks true",
        self_normalize=config['self_normalize'],
        license=license_text,
        reference_dir=config['reference_dir'],
        sample="{sample}"
    log: "logs/{sample}/{sample}.cnv_calling.log"
    group: "cnv_calling"
    wrapper: "file:wrapper/dragen/cnvCalling"

rule cnv_to_remote:
    input:
        tn=rules.cnv_calling.output.tn,
        cnv_png=rules.cnv_calling.output.cnv_png
    output:
        tn=S3.remote(config['s3']['output'] + '/{sample}/{sample}.tn.tsv'),
        cnv_png=S3.remote(config['s3']['output'] + '/{sample}/{sample}.cnv.png')
    log: "logs/{sample}/{sample}.cnv_to_remote.log"
    params:
        sample="{sample}"
    group: "cnv_calling"
    wrapper: "file:wrapper/s3/upload"
