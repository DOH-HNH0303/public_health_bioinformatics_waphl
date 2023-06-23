version 1.0

task python {
  input {
    File kraken2_report
    String docker_image = "amancevice/pandas:1.4.2-alpine"
    Int mem_size_gb = 8
    Int CPUs = 4
  }
  command <<<
  ls>list.txt
  head ~{kraken2_report}


  >>>
  output {
    File list = "list.txt"
    String python_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{mem_size_gb} GB"
    cpu: CPUs
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task get_dt_results {
  input {
    String docker_image = "broadinstitute/terra-tools:tqdm"
    String samplename
    Int mem_size_gb = 8
    Int CPUs = 4
    File?    tblastn_dt_omega_report
    File?    tblastn_dt_beta_report
    File?    tblastn_dt_beta_homologue_report
  }
  command <<<
  python <<CODE
  print("test 1")
  dt_array=["P00587", "P00588", "P00589"]
  file_array=["~{tblastn_dt_omega_report}", "~{tblastn_dt_beta_report}", "~{tblastn_dt_beta_homologue_report}"]
  for i in range(len(file_array)):

    with open(file_array[i], 'r') as file:
      count = 0
      print("test 2")
      for line in file:
        count += 1
        print("test 3")
        if count==1:
          print("test 4")
          data_array = line.split()
          eval=data_array[-2]
          bitscore=data_array[-1]
          print("test 5")

          if float(eval)<=1e-50:
             text="positive"
             print("test 6")
          elif float(data) <=0.01:
             text="possible homolog"
             print("test 7")
          else:
             text="negative"
             print("test 8")
      print("test 9")
      if count==0:
        text="negative"
        eval="NULL"
        bitscore="NULL"
        print("test 10")

    print("test 11")
    eval_name=dt_array[i]+"_EVALUE"
    print("test 11.2")
    f = open(eval_name, "w")
    print("test 11.3")
    f.write(eval)
    print("test 11.4")
    f.close()

    print("test 12")
    bitscore_name=dt_array[i]+"_BITSCORE"
    f = open(bitscore_name, "w")
    f.write(bitscore)
    f.close()

    print("test 13")
    result_name=dt_array[i]+"_RESULT"
    f = open(result_name, "w")
    f.write(text)
    f.close()

  CODE
  ls


  >>>
  output {
    String dt_omega =read_string("P00587_RESULT")
    String dt_beta =read_string("P00588_RESULT")
    String dt_beta_homologue =read_string("P00589_RESULT")
    String? dt_omega_evalue =read_string("P00587_EVALUE")
    String? dt_beta_evalue =read_string("P00588_EVALUE")
    String? dt_beta_homologue_evalue =read_string("P00589_EVALUE")
    String? dt_omega_bitscore =read_string("P00587_BITSCORE")
    String? dt_beta_bitscore =read_string("P00588_BITSCORE")
    String? dt_beta_homologue_bitscore =read_string("P00589_BITSCORE")

    String wget_dt_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{mem_size_gb} GB"
    cpu: CPUs
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task join_genus_species {
  input {
    String? genus
    String? species
    Int mem_size_gb = 8
  }
  command <<<
  echo "~{genus} ~{species}" >GENUS_SPECIES
  >>>
  output {
    String genus_species = read_string("GENUS_SPECIES")
  }
  runtime {
    docker: "quay.io/theiagen/utility:1.1"
    memory: "~{mem_size_gb} GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task split_by_clade {
  input {
    File snp_matrix
    String cluster_name
    String docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    Int threads = 6
    Int snp_clade
  }
  command <<<
    # date and version control
    date | tee DATE
    python3<<CODE

    import pandas as pd

    '''Takes snp_distance_matrix.tsv as input and returns txt file with list
    of '''
    input = "~{snp_matrix}"
    output = "~{cluster_name}_output.txt"
    #cluster_dist = 150 #SNP ingdistance determines cluster
    cluster_dist = ~{snp_clade}

    if not cluster_dist:
      cluster_dist = 150

    df =pd.read_csv(input, sep='\t', header=0)

    seqs=df[df.columns[0]].tolist()
    df = df.replace("-",0)

    done_list=[]
    seq_list =[]
    for i in seqs:
      snp_dist=df[i].tolist()
      res = [idx for idx, val in enumerate(snp_dist) if int(val) <= cluster_dist]

      ids=[]
      for j in res:
        val=df.iloc[j].loc[df.columns[0]]
        ids.append(val)

      if res not in done_list:
        done_list.append(res)
        seq_list.append(ids)

    with open(output, 'w') as fp:
        for li in seq_list:
            for item in li:
                fp.write("%s\t" % item)
            fp.write("\n")
        print('Done')

    CODE
  >>>
  output {
    String date = read_string("DATE")
    File clade_list_file = "~{cluster_name}_output.txt"
    Array[Array[String]] clade_list = read_tsv("~{cluster_name}_output.txt")
    String split_clade_docker_image = docker
  }
  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task scatter_by_clade {
  input {
    Array[File] assembly_files
    String filetype = "gff"
    String cluster_name
    String docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    Int threads = 6
    Array[String] clade_list
  }
  command <<<
    date | tee DATE
    mkdir files_dir
    for x in ~{sep=' ' assembly_files}
    do
        mv "${x}" ./$(basename "${x}")
    done;
    ls
    echo ""
    for x in ~{sep=' ' clade_list}
    do
        if [ ~{filetype} == "fasta" ]; then
            mv "${x}_contigs.fasta" files_dir/"${x}_contigs.fasta"
        elif [ ~{filetype} == "gff" ]; then
            mv "${x}.gff" files_dir/"${x}.gff" 
        else
            echo "Please add filetype to task"
            ls "${x}"*
        fi
    done;
    echo "ls files_dir"
    ls files_dir
    cd files_dir
    python3<<CODE

    import os
    # assign directory
    directory = '.'
    filetype = "gff"
    # iterate over files in
    # that directory
    with open("file_list.txt", "w") as file1:
        # Writing data to a file
      for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
          if filetype == "gff":
            if f.endswith("gff"):
              if f.startswith("./"):
                f = f.split("./")[-1]
              f = f.split('.')[0]
              print(f)
              file1.write(f+"\n")
          elif filetype == "fasta":
            if f.endswith("fasta") or f.endswith("fna") or f.endswith("aln"):
              if f.startswith("./"):
                f = f.split("./")[-1]
              f = f.split('.')[0]
              f = f.split('_contigs')[0]
              print(f)
              file1.write(f+"\n")
      file1.close()

    CODE
    mv file_list.txt ../file_list.txt
  >>>
  output {
    String date = read_string("DATE")
    Array[File] clade_files = glob("files_dir/*")
    Array[String] samplename = read_lines("file_list.txt")
    String scatter_clade_docker_image = docker
  }
  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task generate_none {
  input {
    String docker = "quay.io/broadinstitute/py3-bio:0.1.2"

  }
  command <<<
    date | tee DATE
    touch none.txt
  >>>
  output {
    String date = read_string("DATE")
    File none_file = "none.txt"
    String none_string = ""
    String none_docker_image = "~docker"
  }
  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task concat_fastq {
  input {
    String docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    File read1_cleaned
    File read2_cleaned
    String samplename

  }
  command <<<
    date | tee DATE
    cat ~{read1_cleaned} ~{read2_cleaned}>~{samplename}_total_reads.fastq
  >>>
  output {
    String date = read_string("DATE")
    File comb_fastq = "~{samplename}_total_reads.fastq"
    String none_docker_image = "~docker"
  }
  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
