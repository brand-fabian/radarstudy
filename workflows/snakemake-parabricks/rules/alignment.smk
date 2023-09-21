rule fq2bam:
  input: unpack(get_fq2bam_input)
  output:
    bam="output/{sample}/{sample}.bam",
    bai="output/{sample}/{sample}.bam.bai",
    recal="output/{sample}/{sample}.recal_out",
    dups_metrics="output/{sample}/{sample}.duplicate_metrics.txt",
  params:
    sample="{sample}",
    reference=config["ref"]["genome"],
    extra=" ".join([
      f"--knownSites {f}" for f in config["ref"]["known-variants"]
    ])
  threads: 24
  resources:
    mem_mb=384000,
    runtime=150,
    gres="localtmp:200G,gpu:1",
    licenses="parabricks@nvlm.meb.uni-bonn.de:1",
  log: "logs/fq2bam/{sample}.log"
  wrapper: "file:wrapper/parabricks/fq2bam"
