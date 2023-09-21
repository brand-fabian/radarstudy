__author__ = 'Fabian Brand'
__email__ = 'brand@imbie.uni-bonn.de'

from snakemake.shell import shell
import os

shell.executable('/bin/bash')

fqs = snakemake.input.get('fq', None)
bam = snakemake.input.get('bam', None)
sample = snakemake.params.get('sample', 'SM0')
output_dir = os.path.dirname(snakemake.output[0])

if fqs is None and bam is not None:
    input_cmd = '-b {} --pair-by-name true'.format(bam)

    # Bam input does not allow read groups settings
    read_groups = ''
else:
    # Fq is set
    if len(fqs) == 2:
        # Default paired end processing
        input_cmd = '-1 {fq1} -2 {fq2}'.format(
            fq1=fqs[0],
            fq2=fqs[1]
        )
        read_groups = "--RGID {sample} --RGLB {sample} --RGPL Illumina --RGSM {sample} --vc-sample-name {sample}".format(
            sample=snakemake.params.get('sample', 'SM0')
        )
    elif len(fqs) > 2:
        # Use input list for multiple fastqs
        if not len(fqs) % 2 == 0:
            raise Exception("Please provide paired end fastq sets.")
        fastq_pairs = set()
        for fq in fqs:
            if fq.endswith('_1.fq.gz'):
                for other in fqs:
                    if other == fq.replace('_1.fq.gz', '_2.fq.gz'):
                        fastq_pairs.add((fq, other))
            elif fq.endswith('_R1_001.fastq.gz'):
                for other in fqs:
                    if other == fq.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'):
                        fastq_pairs.add((fq, other))
        fastq_list = os.path.join(output_dir, f'{sample}.input_list.csv')
        # Make output directory
        if not os.path.exists(os.path.dirname(fastq_list)):
            os.makedirs(os.path.dirname(fastq_list))
        with open(fastq_list, 'w') as fd:
            fd.write('RGID,RGSM,RGLB,RGPL,Lane,Read1File,Read2File\n')
            for pair in fastq_pairs:
                fd.write(f'{sample},{sample},WGS,Illuminna,1,{pair[0]},{pair[1]}\n')
        
        input_cmd = f'--fastq-list {fastq_list} --fastq-list-sample-id {sample}'
        read_groups = ''

shell("""
    /opt/edico/bin/dragen -r {snakemake.params.reference_dir} \
        {input_cmd} \
        --output-file-prefix {snakemake.params.sample} \
        --output-directory {output_dir} \
        --enable-variant-caller true \
        --vc-emit-ref-confidence GVCF \
        --enable-duplicate-marking true \
        --enable-map-align-output true \
        {read_groups} \
        --lic-server {snakemake.params.license} 2>&1 > {snakemake.log}
""")
