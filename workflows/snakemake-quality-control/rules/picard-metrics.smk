if metrics_source == METRICS_SOURCE.PICARD:
    rule multiple_metrics:
        input: unpack(get_bam)
        output: **get_multiple_metrics_output(METRICS_SOURCE.PICARD)
        log: "logs/multiple_metrics/{sample}.log"
        params:
            reference=config['ref']['genome']
        threads: 2
        resources:
            runtime=600,
            mem_mb=24000,
        wrapper: "file:snakemake-wrapper/picard/CollectMultipleMetrics"

    rule wgs_metrics:
        input: unpack(get_bam)
        output: "output/{sample}/{sample}.wgs_metrics.txt"
        log: "logs/wgs_metrics/{sample}.log"
        params:
            reference=config['ref']['genome']
        threads: 2
        resources:
            runtime=600,
            mem_mb=24000,
        wrapper: "file:snakemake-wrapper/picard/CollectWgsMetrics"
elif metrics_source == METRICS_SOURCE.PARABRICKS:
    rule multiple_metrics:
        input: unpack(get_bam)
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
        wrapper: "file:snakemake-wrapper/parabricks/collectmultiplemetrics"

    rule wgs_metrics:
        input: unpack(get_bam)
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
        wrapper: "file:snakemake-wrapper/parabricks/bammetrics"
