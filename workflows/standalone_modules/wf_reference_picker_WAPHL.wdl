version 1.0

workflow basespace_fetch {

  input {
    String    sample_name
    String    dataset_name
    String    api_server
    String    access_token
  }

  call fetch_bs {
    input:
      sample=sample_name,
      dataset=dataset_name,
      api=api_server,
      token=access_token
  }

  output {
    File    read1   =fetch_bs.read1
    File    read2   =fetch_bs.read2
  }
}

task fetch_bs {

  input {
    String    sample
    String    dataset
    String    api
    String    token
  }

  command <<<

    bs --api-server=~{api} --access-token=~{token} download biosample -n ~{dataset} -o .

    mv ~{dataset}*ds*/*_R1_* ~{sample}_R1.fastq.gz
    mv ~{dataset}*ds*/*_R2_* ~{sample}_R2.fastq.gz

  >>>

  output {
    File    read1="${sample}_R1.fastq.gz"
    File    read2="${sample}_R2.fastq.gz"
  }

  runtime {
    docker:       "theiagen/basespace_cli:1.2.1"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
