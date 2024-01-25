version 1.0

task iqtree_old {
  input {
    File alignment
    String cluster_name
    String iqtree_model = "GTR+I+G" # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    Int iqtree_bootstraps = 1000 #  Ultrafast bootstrap replicates
    Int alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String iqtree_opts = ""
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/iqtree:1.6.7"
    Int disk_size = 100
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE
    iqtree --version | grep version | sed 's/.*version/version/;s/ for Linux.*//' | tee VERSION

    numGenomes=`grep -o '>' ~{alignment} | wc -l`
    if [ "$numGenomes" -gt 3 ]
    then
      cp ~{alignment} ./msa.fasta
      iqtree \
      -nt AUTO \
      -s msa.fasta \
      -m ~{iqtree_model} \
      -bb ~{iqtree_bootstraps} \
      -alrt ~{alrt} \
      ~{iqtree_opts}

      cp msa.fasta.contree ~{cluster_name}_iqtree.tree
    fi
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File ml_tree = "~{cluster_name}_iqtree.tree"
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}


task iqtree {
  input {
    File? alignment
    String cluster_name
    String iqtree_model = "GTR+I+G" # For comparison to other tools use HKY for bactopia, GTR+F+I for grandeur, GTR+G4 for nullarbor, GTR+G for dryad
    String iqtree_bootstraps = 1000 #  Ultrafast bootstrap replicates
    String alrt = 1000 # SH-like approximate likelihood ratio test (SH-aLRT) replicates
    String? iqtree_opts = ""
    String docker = "staphb/iqtree:1.6.7"
    String? details
  }
  command <<<
    # date and version control
    date | tee DATE
    iqtree --version | grep version | sed 's/.*version/version/;s/ for Linux.*//' | tee VERSION

    numGenomes=`grep -o '>' ~{alignment} | wc -l`
    if [ $numGenomes -gt 3 ]
    then
      cp ~{alignment} ./msa.fasta
      iqtree \
      -nt AUTO \
      -s msa.fasta \
      -m ~{iqtree_model} \
      -bb ~{iqtree_bootstraps} \
      -alrt ~{alrt} \
      ~{iqtree_opts}>> terminal_output.txt 2>&1|| \
      iqtree \
      -nt AUTO \
      -s msa.fasta \
      -m "GTR+I+G" \
      -bb ~{iqtree_bootstraps} \
      -alrt ~{alrt} \
      ~{iqtree_opts} >> terminal_output.txt 2>&1 || \
      iqtree \
      -s msa.fasta \
      -m "GTR+I+G" \
      ~{iqtree_opts} >> terminal_output.txt 2>&1

      echo "test1"
      cp msa.fasta.contree ~{details}~{cluster_name}_msa.tree || touch none.tree
      cp msa.fasta.iqtree ~{details}~{cluster_name}_msa.iqtree || touch none.iqtree
    fi
    ls

    echo "test2"
    if [ -f "~{details}~{cluster_name}_msa.iqtree" ]; then
    echo "test3"
    if grep -q "Model of substitution:" "~{details}~{cluster_name}_msa.iqtree"; then
      echo "test4"
      cat ~{details}~{cluster_name}_msa.iqtree | grep "Model of substitution" | sed -r 's/Model of substitution: //' | tee IQTREE_MODEL # SomeString was found
      echo "test5"
    elif grep -q "Best-fit model according to BIC" ~{details}~{cluster_name}_msa.iqtree; then
      echo "test6"
      cat ~{details}~{cluster_name}_msa.iqtree | grep "Best-fit model according to BIC" | sed -r 's/Best-fit model according to BIC//' | tee IQTREE_MODEL
      echo "test7"
    else
      echo "test8"
      echo ~{iqtree_model} | tee IQTREE_MODEL
      echo "test9"
    fi
    fi

    if [ ! -f IQTREE_MODEL ]; then
    touch IQTREE_MODEL
    fi
    if [ ! -f *.tree ]; then
    echo "test9.2"
    touch none.tree
    ls none.tree
    fi
    echo "test10"
    if [ ! -f *.iqtree ]; then
    echo "test11"
    touch none.iqtree
    echo "test12"
    fi
    echo "test13"
    touch IQTREE_COMMENT
    if [[ -f "terminal_output.txt" ]]; then
        echo "test14"
        if grep -q "WARNING: Your alignment contains too many identical sequences!" terminal_output.txt; then
            echo "test15"
            echo "Too few unique sequences to generate tree">>IQTREE_COMMENT
        elif grep -q "ERROR: It makes no sense to perform bootstrap with less than 4 sequences" terminal_output.txt; then
            echo "test15"
            echo "Too few unique sequences to perform bootstrapping">>IQTREE_COMMENT
        else
            echo "test16"
            touch IQTREE_COMMENT
        fi
    echo "test17"
    elif [ $numGenomes -le 3 ]; then
        echo "test18"
        echo "Too few unique sequences to generate tree">>IQTREE_COMMENT
    fi
    echo "test end"

    ls
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File? iqtree_terminal = "terminal_output.txt"#select_first(["terminal_output3.txt", "terminal_output2.txt", "terminal_output1.txt"])
    File? ml_tree = select_first(glob("*.tree")) # ["~{cluster_name}_msa.tree", "none.tree"])
    File? iqtree_report = select_first(["~{details}~{cluster_name}_msa.iqtree", "none.iqtree"]) #glob("*.iqtree"))
    String? iqtree_model_used = read_string("IQTREE_MODEL")
    String? iqtree_comment = read_string("IQTREE_COMMENT")
  }
  runtime {
    docker: "~{docker}"
    memory: "32 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
