def ped_input(wildcards):
    if 'fastq1' in samples.iloc[0]:
        return S3.remote(samples.iloc[0].fastq1, stay_on_remote=True)
    else:
        return S3.remote(input_files.iloc[0].path, stay_on_remote=True)

rule pedigree:
    input: ped_input # Placeholder input
    output: "output/{family}/{family}.ped"
    params:
        family="{family}"
    run:
        family = samples[samples['family_id'] == params.family]
        if not os.path.exists(os.path.dirname(output[0])):
            os.path.makedirs(os.path.dirname(output[0]))
        with open(output[0], 'w+') as out_f:
            for idx, row in family.iterrows():
                out_f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    row.family_id,
                    row.sample_id,
                    row.father_id if is_valid(row.father_id) else '0',
                    row.mother_id if is_valid(row.mother_id) else '0',
                    transform_sex(row.sex) if is_valid(row.sex) else '0',
                    row.is_affected if is_valid(row.is_affected) else '0', 
                ))
        shell("aws s3 cp {output} s3://{bucket}/{output}".format(
            output=output[0],
            bucket=config['s3']['output'].strip('s3://')
        ))