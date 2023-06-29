version 1.0

task ska {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    File reference
    String cluster_name
    Int kmer_size = 19
    String docker_image = "staphb/ska:1.0"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
  assembly_array=(~{sep=' ' assembly_fasta})
  echo "${#assembly_array[@]}">assembly.txt
  assembly_array_len=$(echo "${#assembly_array[@]}")
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")

  # Ensure assembly, and samplename arrays are of equal length
  if [ "$assembly_array_len" -ne "$samplename_array_len" ]; then
    echo "Assembly array (length: $assembly_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
    exit 1
  fi

  ls /data
  echo "ls /data"

  # create file of filenames for kSNP3 input
  # for f in *fasta; do echo '${f%.fasta}    ${f}'; done > kingc_isolates.list

  touch ~{cluster_name}_isolates.txt
  for f in ${!assembly_array[@]}; do
    assembly=${assembly_array[$f]}
    samplename=${samplename_array[$f]}
    echo -e "${samplename}    ${assembly}" >> ~{cluster_name}_isolates.txt
  done

  #for f in ${!assembly_array[@]};do echo "${f%.fna}    ${f}"; done >  ~{cluster_name}_isolates.txt && \
  ls
  python3 /data/generate_ska_alignment.py --reference ~{reference} --fasta ~{cluster_name}_isolates.txt --out ~{cluster_name}.aln


  >>>
  output {
    File ska_aln = "~{cluster_name}.aln"
    File ska_input = "~{cluster_name}_isolates.txt"
    String ska_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 3
    continueOnReturnCode: "True"
  }
}
