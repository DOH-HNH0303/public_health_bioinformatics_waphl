version 1.0

task ksnp4 {
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String cluster_name
    Int kmer_size = 19
    String ksnp4_args = "" # add -ML to calculate a maximum likelihood tree or -NJ to calculate a neighbor-joining tree
    String docker_image = "staphb/ksnp4:4.0"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
  assembly_array=(~{sep=' ' assembly_fasta})
  assembly_array_len=$(echo "${#assembly_array[@]}")
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")
  
  # Ensure assembly, and samplename arrays are of equal length
  if [ "$assembly_array_len" -ne "$samplename_array_len" ]; then
    echo "Assembly array (length: $assembly_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
    exit 1
  fi
  # ensure kSNP file naming convention is met
  assembly_renamed_array=()
  for index in ${!assembly_array[@]}; do
      assembly=${assembly_array[$index]}
      # ensure kSNP file naming convention is met by removing non-id, dot-separated info, 
      # e.g. sample01.ivar.consensus.fasta will be renamed to sample01.fasta
      assembly_renamed=$(echo $assembly | sed 's/\.\(.*\)\././')
      echo "ASSEMBLY: $assembly"
      echo "ASSEMBLY_RENAMED: $assembly_renamed"
      mv $assembly $assembly_renamed
      assembly_renamed_array+=($assembly_renamed)
  done

  # create file of filenames for kSNP4 input
  touch ksnp4_input.tsv
  for index in ${!assembly_renamed_array[@]}; do
    assembly=${assembly_renamed_array[$index]}
    samplename=${samplename_array[$index]}
    echo -e "${assembly}\t${samplename}" >> ksnp4_input.tsv
  done
  
  echo "ksnp4_input.tsv:: "
  cat ksnp4_input.tsv

  # run ksnp4 on input assemblies
  kSNP4 -in ksnp4_input.tsv -outdir ksnp4 -k ~{kmer_size} -core -vcf ~{ksnp4_args} -debug

  echo "ls ksnp4"
  ls ksnp4
  echo ""
  
  # rename ksnp4 outputs with cluster name 
  # sometimes the core nwk and fasta outputs do not have content
  mv -v ksnp4/core_SNPs_matrix.fasta ksnp4/~{cluster_name}_core_SNPs_matrix.fasta
  mv -v ksnp4/tree.core.tre ksnp4/~{cluster_name}_core.nwk

  if [ -s ksnp4/~{cluster_name}_core_SNPs_matrix.fasta ]; then # is the file not-empty?
    echo "The core SNP matrix was produced" | tee SKIP_SNP_DIST # then do NOT skip
  else
    echo "The core SNP matrix could not be produced" | tee SKIP_SNP_DIST # otherwise, skip
  fi

  mv -v ksnp4/VCF.*.vcf ksnp4/~{cluster_name}_core.vcf
  mv -v ksnp4/SNPs_all_matrix.fasta ksnp4/~{cluster_name}_pan_SNPs_matrix.fasta
  mv -v ksnp4/tree.parsimony.tre ksnp4/~{cluster_name}_pan_parsimony.nwk

  if [ -f ksnp4/tree.ML.tre ]; then  
    mv -v ksnp4/tree.ML.tre ksnp4/~{cluster_name}_ML.nwk
  fi 
  if [ -f ksnp4/tree.NJ.tre ]; then  
    mv -v ksnp4/tree.NJ.tre ksnp4/~{cluster_name}_NJ.nwk
  fi 

  >>>
  output {
    File ksnp4_core_matrix = "ksnp4/${cluster_name}_core_SNPs_matrix.fasta"
    File ksnp4_core_tree = "ksnp4/${cluster_name}_core.nwk"
    File ksnp4_core_vcf = "ksnp4/${cluster_name}_core.vcf"
    File ksnp4_pan_matrix = "ksnp4/~{cluster_name}_pan_SNPs_matrix.fasta"
    File ksnp4_pan_parsimony_tree = "ksnp4/~{cluster_name}_pan_parsimony.nwk"
    File? ksnp4_ml_tree = "ksnp4/~{cluster_name}_ML.nwk"
    File? ksnp4_nj_tree = "ksnp4/~{cluster_name}_NJ.nwk"
    File number_snps = "ksnp4/COUNT_SNPs"
    File ksnp4_input = "ksnp4_input.tsv"
    String skip_core_snp_dists = read_string("SKIP_SNP_DIST")
    Array[File] ksnp_outs = glob("ksnp4/*")
    String ksnp4_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 0
  }
}
