rule mosdepth:
    input: unpack(get_bams)
    output:
        dist="output/{sample}/{sample}.mosdepth.global.dist.txt",
        summary="output/{sample}/{sample}.mosdepth.summary.txt",
        regions="output/{sample}/{sample}.regions.bed.gz",
        regions_dist="output/{sample}/{sample}.mosdepth.region.dist.txt",
        thresholds="output/{sample}/{sample}.thresholds.bed.gz",
    log: "logs/{sample}/{sample}.mosdepth.log"
    params:
        extra=("--no-per-base --fast-mode --by {refseq}".format(refseq=config['ref']['refseq'])
            + " --thresholds {}".format(config['mosdepth']['thresholds']) if 'thresholds' in config['mosdepth'] else ""
            + " --quantize {}".format(config['mosdepth']['quantize']) if 'quantize' in config['mosdepth'] else "")
    threads: 4
    resources:
        runtime=300,
        mem_mb=32000
    wrapper: "file:snakemake-wrapper/mosdepth"