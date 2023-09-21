rule whatshap:
  input: unpack(get_whatshap_input)
  output:
    vcf="output/{family_id}/{family_id}.phased.vcf.gz"
  params:
    family_id="{family_id}",
    reference=config["ref"]["grch37"]["genome"]
  threads: 2
  resources:
    mem_mb=32000,
    runtime=720
  log: "logs/whatshap/{family_id}.log"
  wrapper: "file:wrapper/whatshap"
