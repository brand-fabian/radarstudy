rule subset_bam:
  input: unpack(get_bam)
  output:
    bam="output/{family_id}/{family_id}.bam"
  params:
    sites=config['sites'],
    pedigree=config['pedigree'],
    grch37=config['ref']['grch37']['genome'],
    grch38=config['ref']['grch38']['genome'],
    family_id="{family_id}",
    flanking_region=config["flanking-region"]
  threads: 12
  resources:
    mem_mb=32000,
    runtime=720
  log: "logs/subset_bam/{family_id}.log"
  wrapper: "file:scripts/subset_bam"

rule subset_sample_bam:
  input: unpack(get_sample_bam)
  output:
    bam="output/{sample}/{sample}.bam"
  params:
    sites=config['sites'],
    pedigree=config['pedigree'],
    grch37=config['ref']['grch37']['genome'],
    grch38=config['ref']['grch38']['genome'],
    sample="{sample}",
    flanking_region=config["flanking-region"]
  threads: 12
  resources:
    mem_mb=32000,
    runtime=720
  log: "logs/subset_bam/{sample}.log"
  wrapper: "file:scripts/subset_bam"

rule subset_vcf:
  input: unpack(get_cohort_vcf)
  output:
    vcf="output/{family_id}/{family_id}.sites.vcf.gz"
  params:
    sites=config['sites'],
    pedigree=config['pedigree'],
    family_id="{family_id}",
    flanking_bases=config['flanking-region'],
  threads: 2
  resources:
    mem_mb=8192,
    runtime=60
  log: "logs/subset_vcf/{family_id}.log"
  wrapper: "file:scripts/subset_vcf"

rule dnm_vcf:
  input: unpack(get_cohort_vcf)
  output:
    vcf="output/{family_id}/{family_id}.dnm_sites.vcf.gz"
  params:
    sites=config['sites'],
    pedigree=config['pedigree'],
    family_id="{family_id}",
    flanking_bases=3,
  threads: 2
  resources:
    mem_mb=8192,
    runtime=60
  log: "logs/subset_vcf/{family_id}.log"
  wrapper: "file:scripts/subset_vcf"