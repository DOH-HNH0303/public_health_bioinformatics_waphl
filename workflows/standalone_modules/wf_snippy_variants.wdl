version 1.0

import "../../tasks/gene_typing/task_snippy_variants.wdl" as snippy
import "../../tasks/task_versioning.wdl" as versioning

workflow snippy_variants_wf {
  meta {
    description: "Perform SNP analysis using snippy"
  }
  input {
    File reference
    File read1
    File? read2
    String samplename
    String? docker
  }
  call snippy.snippy_variants {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      reference = reference,
      docker = docker
  }
  call versioning.version_capture{
    input:
  }
  output {
    String snippy_variants_wf_version = version_capture.phb_version
    String snippy_version = snippy_variants.snippy_variants_version
    String snippy_variants_query = snippy_variants.snippy_variants_query
    String snippy_variants_hits = snippy_variants.snippy_variants_hits
    File snippy_variants_gene_query_results = snippy_variants.snippy_variants_gene_query_results
    File snippy_results = snippy_variants.snippy_variants_results
    File snippy_bam = snippy_variants.snippy_variants_bam
    File snippy_bai = snippy_variants.snippy_variants_bai
    File snippy_variants_summary = snippy_variants.snippy_variants_summary
    File snippy_variants_outdir_tarball = snippy_variants.snippy_variants_outdir_tarball
  }
}
