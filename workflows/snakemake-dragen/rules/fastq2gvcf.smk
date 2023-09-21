rule fastq2gvcf:
    input: unpack(get_vc_input)
    output:
        hard_filtered="output/{sample}/{sample}.hard-filtered.gvcf.gz",
        hard_filtered_idx="output/{sample}/{sample}.hard-filtered.gvcf.gz.tbi",
        hard_filtered_md5="output/{sample}/{sample}.hard-filtered.gvcf.gz.md5sum",
        bam="output/{sample}/{sample}.bam",
        bai="output/{sample}/{sample}.bam.bai",
        bam_md5="output/{sample}/{sample}.bam.md5sum"
    params:
        reference_dir=config['reference_dir'],
        sample="{sample}",
        license=license_text
    log: "logs/{sample}/fastq2gvcf.log"
    group: "dragen_gvcf"
    wrapper: "file:wrapper/dragen/fastq2gvcf"

rule gvcf2remote:
    input: 
        hard_filtered=rules.fastq2gvcf.output.hard_filtered,
        hard_filtered_idx=rules.fastq2gvcf.output.hard_filtered,
        hard_filtered_md5=rules.fastq2gvcf.output.hard_filtered_md5,
        bam=rules.fastq2gvcf.output.bam,
        bai=rules.fastq2gvcf.output.bai,
        bam_md5=rules.fastq2gvcf.output.bam_md5
    output:
        hard_filtered=S3.remote(config['s3']['output'] + '/{sample}/{sample}.hard-filtered.gvcf.gz'),
        hard_filtered_idx=S3.remote(config['s3']['output'] + '/{sample}/{sample}.hard-filtered.gvcf.gz.tbi'),
        bam=S3.remote(config['s3']['output'] + '/{sample}/{sample}.bam'),
        bai=S3.remote(config['s3']['output'] + '/{sample}/{sample}.bam.bai'),
    log: "logs/{sample}/copy_to_remote.log"
    params:
        sample="{sample}"
    group: "dragen_gvcf"
    wrapper: "file:wrapper/s3/upload"
