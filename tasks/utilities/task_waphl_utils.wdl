version 1.0

task freyja_epi_metadata {
  input {
    String samplename
    File ww_metadata_csv
    Int memory = 4
    String docker = "staphb/freyja:1.3.10"
  }
  command <<<
  # capture version

  python3 <<CODE
  import pandas as pd

  accession="~{samplename}"
  accession=accession[0:9]
  file="~{ww_metadata_csv}"
  df = pd.read_csv(file)
  df = df[["PHLAccessionNumber", "SubmitterSampleNumber", "WWTPName",  "SampleCollectDate"]]
  print(df)

  if accession in df['PHLAccessionNumber'].unique():
    submitter_num=df.loc[df['PHLAccessionNumber'] == accession, 'SubmitterSampleNumber'].iloc[0]
    wwtpname=df.loc[df['PHLAccessionNumber'] == accession, 'WWTPName'].iloc[0]
    collection_date=df.loc[df['PHLAccessionNumber'] == accession, 'SampleCollectDate'].iloc[0]

    submitter_file = open("SUBMITTER_NUM", "w")
    n = submitter_file.write(str(submitter_num))
    submitter_file.close()

    wwtpname_file = open("WWTPNAME", "w")
    n = wwtpname_file.write(str(wwtpname))
    wwtpname_file.close()

    collection_date_file = open("COLLECTION_DATE", "w")
    n = collection_date_file.write(str(collection_date))
    collection_date_file.close()

  else:
    submitter_file = open("SUBMITTER_NUM", "w")
    n = submitter_file.write("")
    submitter_file.close()

    wwtpname_file = open("WWTPNAME", "w")
    n = wwtpname_file.write("")
    wwtpname_file.close()

    collection_date_file = open("COLLECTION_DATE", "w")
    n = collection_date_file.write("")
    collection_date_file.close()


  CODE

  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: 2
    docker: "~{docker}"
    disks: "local-disk 100 HDD"
    maxRetries: 3
  }
  output {
    String? SubmitterSampleNumber = read_string("SUBMITTER_NUM")
    String? WWTPName = read_string("WWTPNAME")
    String? SampleCollectDate = read_string("COLLECTION_DATE")
    String epi_metadata_docker = "~{docker}"
  }
}

task freyja_epi_output {
  input {
    String samplename
    String? samplecollectdate
    String? wwtpname
    String? submittersamplenumber
    String? wwtp_county
    String? wwtp_zip
    File freyja_demixed
    File freyja_depths
    Int memory = 4
    String docker = "staphb/freyja:1.3.10"
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION
  # update freyja reference files if specified

  awk '{ total += $4; count++ } END { print total/count }' ~{freyja_depths} | tee AVG_COVERAGE


  python3 <<CODE

  import csv
  import fileinput
  import pandas as pd
  import numpy as np
  from datetime import date


  today = date.today()
  print("today", today)

  id = "~{samplename}"
  sc_date = "~{samplecollectdate}"
  print("sc_date", sc_date)
  location = "~{wwtpname}"
  submitter = "~{submittersamplenumber}"
  county = "~{wwtp_county}"
  zip = "~{wwtp_zip}"

  with open("~{freyja_demixed}") as f:
    lines = f.readlines()
    for line in lines:
      if "lineages" in line:
        lineages = line.split("\t")[1].strip().split(" ")
        lineages.append("unreportable")
        text_file = open("LINEAGES", "w")
        n = text_file.write(str(lineages))
        text_file.close()

      if "abundances" in line:
        abundances = line.split("\t")[1].strip().split(" ")
        abundances_float = [float(item) for item in abundances]
        unreportable = 1-(sum(abundances_float))
        abundances.append(unreportable)

        float_file = open("UNREPORTABLE", "w")
        n = float_file.write(str(unreportable))
        float_file.close()

        text_file = open("ABUNDANCES", "w")
        n = text_file.write(str(abundances))
        text_file.close()

      if "resid" in line:
        line=line.split("\t")[1]
        #print("resid", line)
        text_file = open("RESID", "w")
        n = text_file.write(line)
        text_file.close()
      if "coverage" in line:
        line=line.split("\t")[1]
        uncov = 100 - float(line)

        float_file = open("UNCOVERED", "w")
        n = float_file.write(str(uncov))
        float_file.close()
        #print("coverage", line)
        text_file = open("COVERAGE", "w")
        n = text_file.write(line)
        text_file.close()

  print(abundances, lineages)
  assert len(abundances) == len(lineages), "error: There should be one relative abundance for every lineage"

  print(len(abundances))
  id_list=[id]*len(abundances)
  date_list=[sc_date]*len(abundances)
  missing_data = "not_missing"
  print("date_list", date_list)
  if sc_date == "":
    print("test1")
    missing_data = "Missing"
    date_list=[None] * len(abundances)
  freyja_date_list=[today]*len(abundances)

  submitter_list=[submitter]*len(abundances)
  if submitter == "":
    missing_data = "Missing"
    submitter_list=[None] * len(abundances)
    
  location_list=[location]*len(abundances)
  if location == "":
    missing_data = "Missing"
    location_list=[None] * len(abundances)

  county_list=[county]*len(abundances)
  if county == "":
    missing_data = "Missing"
    county_list=[None] * len(abundances)

  zip_list=[zip]*len(abundances)
  if zip == "":
    missing_data = "Missing"
    zip_list=[None] * len(abundances)

  epi = open("MISSING_EPI", "w")
  print(missing_data)
  if missing_data == "Missing":
    a = epi.write('Missing Epi Data')
  epi.close()


  df = pd.DataFrame({'PHL_ID':id_list, 'Sample_ID':submitter_list, 'Sample_Collection_date':date_list,
  'Sample_Site':location_list, 'Sample_county':county_list, 'Sample_zip':zip_list, "lineages":lineages, "abundances":abundances, "freyja_date":freyja_date_list})
  df.replace(np.nan, 'NULL', regex=True)
  df.to_csv('~{samplename}_for_epi.tsv', sep="\t", header=False, index=False)
  CODE

  >>>
  runtime {
    memory: "~{memory} GB"
    cpu: 2
    docker: "~{docker}"
    disks: "local-disk 100 HDD"
    maxRetries: 3
  }
  output {
    String freyja_version = read_string("FREYJA_VERSION")
    String? freyja_lineages = read_string("LINEAGES")
    String? freyja_abundances = read_string("ABUNDANCES")
    String? freyja_resid = read_string("RESID")
    String? freyja_coverage = read_string("COVERAGE")
    Float? freyja_avg_coverage = read_float("AVG_COVERAGE")
    Float? freyja_uncovered = read_float("UNCOVERED")
    File? freyja_epi_file = "~{samplename}_for_epi.tsv"
    Float? freyja_unreportable = read_float("UNREPORTABLE")
    String? missing_epi = read_string("MISSING_EPI")
  }
}
