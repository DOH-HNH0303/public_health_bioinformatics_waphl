version 1.0

import "../tasks/phylogenetic_inference/task_lyveset.wdl" as lyveset
import "../tasks/task_versioning.wdl" as versioning

workflow lyvset_workflow {
  input {
    Array[File] read1
    Array[File] read2
    String dataset_name
    File reference_genome
  }
  call lyveset.lyveset {
    input:
      read1 = read1,
      read2 = read2,
      dataset_name = dataset_name,
      reference_genome = reference_genome
  }
  call versioning.version_capture{
    input:
  }
  output {
    String lyveset_wf_version = version_capture.phb_version
    String lyvset_wf_analysis_date = version_capture.date

    String lyvset_docker_image = lyveset.lyveset_docker_image
    File lyveset_distance_matrix = lyveset.lyveset_distance_matrix
    File lyvset_raxmpl_tree = lyveset.lyveset_raxml_tree
  }
}