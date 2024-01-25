version 1.0
workflow ww_for_epi {

  input {
    Array[String]    samplename
    Array[File]    freyja_epi_file
    String?    batch_name

  }

  call batch_data {
    input:
      samplename=samplename,
      freyja_epi_file=freyja_epi_file,
      batch_name=batch_name
  }

  output {
    String    ww_batch_upload=batch_data.batch_tsv
    #String    ww_batch_for_lab=batch_data.lab_batch_tsv

  }
}


task batch_data {
  input {
    Array[String]    samplename
    Array[File]    freyja_epi_file
    String?    batch_name
    Int    total = length(freyja_epi_file)

  }

  command {
    file_array=~{sep='\t' freyja_epi_file}
    echo $file_array
    name_array=(~{sep=' ' samplename})
    touch ~{batch_name}batch_output.tsv
    #touch ~{batch_name}lab_batch_output.tsv

    # echo "Sample_ID	Sample_Collection_Date	Sample_Site	Variant_name	Variant_proportion	Pipeline_date" >~{batch_name}batch_output.tsv
    # for c in '~{sep=" " freyja_epi_file}';do
    #         echo $c
    #         cat $c >>~{batch_name}batch_output.tsv
    #     done

    echo "Lab_ID	Sample_ID	Sample_Collection_Date	Sample_Site	Sample_County	Sample_Zip	Variant_name	Variant_proportion	Pipeline_date" >~{batch_name}batch_output.tsv
    for c in '~{sep=" " freyja_epi_file}';do
            echo $c
            cat $c >>~{batch_name}batch_output.tsv
        done

  }

  output {
    File    batch_tsv="~{batch_name}batch_output.tsv"
    #File    lab_batch_tsv="~{batch_name}lab_batch_output.tsv"

  }

  runtime {
    docker:       "ubuntu:bionic-20220902"
    memory:       "16 GB"
    cpu:          4
    disks:        "local-disk 100 SSD"
    preemptible:  1
    continueOnReturnCode: "True"
  }
}
