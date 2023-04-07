version 1.0

import "../../../tasks/utilities/task_file_handling.wdl" as file_handling
import "../../../tasks/task_versioning.wdl" as versioning

workflow concatenate_column_content {
  input {
    Array[File] files_to_cat
    String concatenated_file_name
    Array[String]? samplename
    Boolean samplename_column = false
  }
  if (!samplename_column) {
    call file_handling.cat_files {
      input:
        files_to_cat = files_to_cat,
        concatenated_file_name = concatenated_file_name
    }
  }
  if (samplename_column) {
    call file_handling.cat_tables {
      input:
        files_to_cat = files_to_cat,
        samplename = samplename,
        concatenated_file_name = concatenated_file_name
    }   
  }
  call versioning.version_capture {
    input:
  }
  output {
    String concatenate_column_content_version = version_capture.phb_version
    String concatenate_column_content_analysis_date = version_capture.date
    File concatenated_files = cat_files.concatenated_files
  }
}