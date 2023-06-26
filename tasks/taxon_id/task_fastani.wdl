version 1.0

task fastANI {

  input {
    String    samplename
    File    assembly
    String    genus

  }

  command <<<

    #ls ../

    if [[ ~{genus} == "Tatlockia" || ~{genus} == "Fluoribacter" ]]; then
        echo "Number is Even"
        fastANI -q ~{assembly} --rl /data/Legionella_list.txt -o fastani_~{samplename}.txt
    else
        fastANI -q ~{assembly} --rl /data/~{genus}_list.txt -o fastani_~{samplename}.txt
    fi

    python3 /data/pull_organism.py fastani_~{samplename}.txt

  >>>

  output {
    File?    fastani_report="fastani_~{samplename}.txt"
    String?    fastani_genus=read_string("TOPGENUS")
    String?    fastani_species=read_string("TOPSPECIES")
    String?    fastani_strain=read_string("TOPSTRAIN")
    Float?    fastani_aniestimate=read_float("ANIESTIMATE")

  }

  runtime {
    docker:       "hnh0303/fastani:1.33-legionella_corynebacterium"
    memory:       "16 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
