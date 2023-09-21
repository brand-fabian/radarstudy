rule unfazed:
  input: unpack(get_unfazed_input)
  output:
    vcf="output/{family_id}/{family_id}.unfazed.vcf.gz"
  params:
    family_id="{family_id}",
    reference=config["ref"]["grch37"]["genome"]
  threads: 12
  resources:
    mem_mb=64000,
    runtime=720
  log: "logs/unfazed/{family_id}.log"
  wrapper: "file:wrapper/unfazed"
