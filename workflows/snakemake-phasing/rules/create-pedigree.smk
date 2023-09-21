rule create_pedigree:
  input:
    ped=config["pedigree"]
  output:
    ped="output/{family_id}/{family_id}.ped"
  params:
    family_id="{family_id}"
  log: "logs/create_pedigree/{family_id}.log"
  threads: 1
  resources:
    runtime=5
  run:
    import os
    import pandas
    if not os.path.exists(os.path.dirname(output[0])):
      os.path.makedirs(os.path.dirname(output[0]), exist_ok=True)
    ped_df = pandas.read_csv(
      input[0],
      sep="\t",
      names=["family_id", "s", "father", "mother", "sex", "is_affected"]
    )
    with open(output[0], 'w+') as out_f:
      for _, row in ped_df[ped_df.family_id == params.family_id].iterrows():
        out_f.write("\t".join([
          row.family_id,
          row.s,
          str(row.father),
          str(row.mother),
          str(row.sex),
          str(row.is_affected),
        ]) + "\n")
