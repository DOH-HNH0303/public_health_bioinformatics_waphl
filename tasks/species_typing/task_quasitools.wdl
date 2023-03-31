version 1.0

task quasitools {
  input {
    File read1
    File? read2
    File mutation_db = "gs://theiagen-public-files/terra/hivgc-files/mutation_db.tsv"
    String samplename
    String docker = "quay.io/biocontainers/quasitools:0.7.0--pyh864c0ab_1"
  }
  command <<<
    # date and version capture
    date | tee DATE

    # unzip reads
    if [[ "~{read1}" == *.gz ]]; then
      gunzip -c ~{read1} > r1.fastq
    else
      mv ~{read1} r1.fastq
    fi
    read1_unzip=r1.fastq

    # do the same on read2
    # unzip file if necessary
    if ~ [ -z ~{read2} ]; then
      if [[ "~{read2}" == *.gz ]]; then
        gunzip -c ~{read2} > r2.fastq
      else
        mv ~{read2} r2.fastq
      fi
      read2_unzip=r2.fastq
    else
      ont="-rt 5"
    fi


    # Print and save version
    quasitools --version > QUASITOOLS_VERSION && sed -i -e 's/^/quasitools /' QUASITOOLS_VERSION

    # Run hydra
    set -e
    quasitools hydra \
      -mf 0.05 \
      "${ont}" \
      -sc 7 \
      -lc 50 \
      -vq 7 \
      -md 10 \
      -ma 1 \
      -me ${read1_unzip} ${read2_unzip} \
      -o "~{samplename}"
  >>>
  output {
    String quasitools_version = read_string("QUASITOOLS_VERSION")
    String quasitools_date = read_string("DATE")
    File coverage_file = "~{samplename}/coverage_file.csv"
    File dr_report = "~{samplename}/dr_report.csv"
    File hydra_vcf = "~{samplename}/hydra.vcf"
    File mutations_report = "~{samplename}/mutation_report.aavf"
  }
  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 4
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}