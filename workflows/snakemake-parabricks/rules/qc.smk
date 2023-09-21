# def multiple_metrics_output(metrics=[
#     # CollectAlignmentSummaryMetrics
#     "alignment.txt",
#     "read_length_histogram.pdf"
#     # CollectInsertSizeMetrics
#     "insert_size.txt",
#     "insert_size.pdf",
#     "insert_size.png",
#     # QualityScoreDistribution
#     "qualityscore.pdf",
#     "qualityscore.png",
#     "qualityscore.txt",
#     # MeanQualityByCycle
#     "mean_quality_by_cycle.pdf",
#     "mean_quality_by_cycle.png",
#     "mean_quality_by_cycle.txt",
#     # CollectBaseDistributionByCycle
#     "base_distribution_by_cycle.txt",
#     "base_distribution_by_cycle.png",
#     "base_distribution_by_cycle.pdf",
#     # CollectGcBiasMetrics
#     "gcbias_0.png",
#     "gcbias_summary.txt",
#     "gcbias.pdf",
#     "gcbias_detail.txt",
#     # CollectSequencingArtifactMetrics
#     "sequencingArtifact.bait_bias_detail_metrics.txt",
#     "sequencingArtifact.bait_bias_summary_metrics.txt",
#     "sequencingArtifact.pre_adapter_detail_metrics.txt",
#     "sequencingArtifact.pre_adapter_summary_metrics.txt",
#     "sequencingArtifact.error_summary_metrics.txt",
# ]):
#     return {
#         metric: "output/{{sample}}/{{sample}}.{metric}.txt".format(metric=metric)
#         for metric in filter(lambda x: x.endswith('txt'), metrics)
#     }


rule wgs_metrics:
  input: unpack(get_alignment)
  output: "output/{sample}/{sample}.wgs_metrics.txt"
  params:
    reference=config["ref"]["genome"]
  threads: 32
  resources:
    runtime=120,
    mem_mb=64000,
    gres="localtmp:100G",
    licenses="parabricks@nvlm.meb.uni-bonn.de:1,gpu:1",
  log: "logs/wgs_metrics/{sample}.log"
  wrapper: "file:wrapper/parabricks/bammetrics"


rule multiple_metrics:
  input: unpack(get_alignment)
  output: **get_multiple_metrics_output(METRICS_SOURCE.PARABRICKS)
  params:
    reference=config["ref"]["genome"]
  threads: 32
  resources:
    runtime=120,
    mem_mb=64000,
    gres="localtmp:100G,gpu:1",
    licenses="parabricks@nvlm.meb.uni-bonn.de:1",
  log: "logs/multiple_metrics/{sample}.log"
  wrapper: "file:wrapper/parabricks/collectmultiplemetrics"