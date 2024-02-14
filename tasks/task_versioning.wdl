version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
  }
  meta {
    volatile: true
  }
  command {
    PHB_Version="PHB v1.3.0"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$PHB_Version" > PHB_VERSION
  }
  output {
    String date = read_string("TODAY")
    String phb_version = read_string("PHB_VERSION")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: docker
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}

task waphl_version_capture {
  input {
    String? input_1
    String? input_2
    String? input_3
    String? input_4
    String? input_5
    String? input_6
    String? input_7
    String? input_8
    String? input_9
    String? input_10
    String? input_11
    String? input_12
    String? input_13
    String? input_14
    String version = "PHBG-WAPHL v1.0.0-beta"
    String? timezone
  }
  meta {
    volatile: true
  }
  command <<<
    #touch input.tsv
    echo "~{input_1}" "~{input_2}" "~{input_3}" "~{input_4}" "~{input_5}" \
    "~{input_6}" "~{input_7}" "~{input_8}" "~{input_9}" "~{input_10}" \
    "~{input_11}" "~{input_12}" "~{input_13}" "~{input_14}">input.tsv

    #echo $docker_array

    #for item in "${!docker_array[@]}"; do
    #  echo $item>>input.tsv
    #done

    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "~{version}" > PHBG_WAPHL_VERSION

    python3<<CODE

    import os

    directory = '.'

    with open("input.tsv", "r") as file1, \
     open('versions.tsv', mode='w') as out_file:
        tool_list = ["pipeline_version", "utilities"]
        version_list = ["~{version}", "1.1"]
        for line in file1:
          for l in line.split(" "):
            l=l.split(":")
            if len(l)>1:
                tool = ''.join(l[:-1]).split("/")[-1]
                version = l[-1]
                if tool.upper().strip("-").strip("_") not in tool_list:
                    tool_list.append(tool.upper().strip("-").strip("_"))
                    version_list.append(version)
        for i in range(len(tool_list)):
          out_file.write(tool_list[i]+"\t"+version_list[i]+"\n")

    file1.close()
    out_file.close()

    CODE
  >>>
  output {
    String date = read_string("TODAY")
    String phbg_waphl_version = read_string("PHBG_WAPHL_VERSION")
    File tool_versions = "versions.tsv"
    File input_file = "input.tsv"
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "quay.io/theiagen/utility:1.1"
    disks: "local-disk 10 HDD"
    continueOnReturnCode: "True"
  }
}
