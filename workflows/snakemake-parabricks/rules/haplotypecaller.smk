rule haplotypecaller:
  input: unpack(get_alignment)
  output:
    gvcf="output/{sample}/{sample}.gvcf.gz"
  params:
    sample="{sample}",
    reference=config["ref"]["genome"],
    recal_in="output/{sample}/{sample}.recal_out",
  threads: 48
  resources:
    runtime=150,
    mem_mb=384000,
    gres="localtmp:200G,gpu:1",
    licenses="parabricks@nvlm.meb.uni-bonn.de:1",
  log: "logs/haplotypecaller/{sample}.log"
  wrapper: "file:wrapper/parabricks/haplotypecaller"