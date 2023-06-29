version 1.0

import "wf_clade_analysis_WAPHL.wdl" as clade_analysis
import "../tasks/phylogenetic_inference/task_ska.wdl" as ska
import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins
import "../tasks/task_versioning.wdl" as versioning
import "wf_ksnp3_WAPHL.wdl" as ksnp3
import "../tasks/utilities/task_utilities.wdl" as utilities
import "../tasks/task_versioning.wdl" as versioning
import "../tasks/utilities/task_summarize_table_waphl.wdl" as summarize

workflow recomb_aware_phylo_analysis {
  input {
    Array[File] assembly_fasta
    Array[File] assembly_gff
    Array[String] samplename
    File reference_genome
    String cluster_name
    # String terra_workspace
    #String terra_project
    # String terra_table
    String iqtree_model = "MFP"
    Int snp_clade = 150
    Float filter_perc = 35.0
  }
  call ska.ska as ska {
    input:
      assembly_fasta = assembly_fasta,
      samplename = samplename,
      cluster_name = cluster_name,
      reference = reference_genome
}
call gubbins.gubbins as gubbins_init {
  input:
    alignment = ska.ska_aln,
    filter_perc = filter_perc,
    cluster_name = cluster_name
}
if (gubbins_init.gubbins_mask == true){
call gubbins.maskrc_svg as mask_gubbins_init  {
  input:
    alignment = ska.ska_aln,
    cluster_name = cluster_name,
    recomb = gubbins_init.recomb_gff,
    base_reconstruct = gubbins_init.base_reconstruct,
    recomb_embl = gubbins_init.recomb_embl,
    polymorph_site_fasta = gubbins_init.polymorph_site_fasta,
    polymorph_site_phylip = gubbins_init.polymorph_site_phylip,
    branch_stats = gubbins_init.branch_stats,
    gubbins_snps = gubbins_init.gubbins_snps,
    gubbins_final_tre = gubbins_init.gubbins_final_tre,
    gubbins_log = gubbins_init.gubbins_log,
    gubbins_node_tre = gubbins_init.gubbins_node_tre
}
}
call ksnp3.ksnp3_workflow as ksnp3  {
  input:
    assembly_fasta = select_first([mask_gubbins_init.masked_fasta_list, assembly_fasta]),
    samplename = samplename,
    cluster_name = cluster_name
}
  call iqtree.iqtree as total_iqtree {
    input:
      alignment = ksnp3.ksnp3_core_matrix,
      cluster_name = cluster_name,
      iqtree_model = iqtree_model
  }
call utilities.split_by_clade as split_by_clade  {
  input:
    snp_matrix = ksnp3.ksnp3_core_snp_matrix,
    cluster_name = cluster_name,
    snp_clade = snp_clade
}
scatter (pair in zip(split_by_clade.clade_list, range(length(split_by_clade.clade_list)))) {
call utilities.scatter_by_clade as scatter_by_clade  {
  input:
    clade_list = pair.left,
    cluster_name = cluster_name,
    assembly_files = assembly_gff
}
call clade_analysis.clade_analysis as clade_analysis  {
  input:
    cluster_name = "~{cluster_name + '_' + pair.right + '_'}clade",
    filter_perc = filter_perc,
    prokka_gff = scatter_by_clade.clade_files,
    samplename = scatter_by_clade.samplename
}

}
call summarize.zip_files as zip_files  {
  input:
    clade_trees = select_all(clade_analysis.clade_iqtree_pan_tree),
    recomb_gff = select_all(clade_analysis.gubbins_clade_recomb_gff),
    pirate_aln_gff = clade_analysis.pirate_aln_pan,
    pirate_gene_presence_absence = select_all(clade_analysis.pirate_for_scoary_csv),
    cluster_name = cluster_name,
    cluster_tree = total_iqtree.ml_tree
    #terra_table = terra_table,
    #terra_workspace = terra_workspace,
    #terra_project = terra_project,
    
}
call versioning.waphl_version_capture as version {
  input:
    input_1 = ska.ska_docker_image,
    input_2 = gubbins_init.gubbins_docker_image,
    input_3 = mask_gubbins_init.maskrc_docker_image,
    input_4 = ksnp3.ksnp3_docker,
    input_5 = total_iqtree.version,
    input_6 = ksnp3.ksnp3_snp_dists_version,
    input_7 = split_by_clade.split_clade_docker_image,
    input_8 = select_first(scatter_by_clade.scatter_clade_docker_image),
    input_9 = select_first(clade_analysis.pirate_docker_image),
    input_11 = select_first(clade_analysis.maskrc_docker_image),
    input_13 = select_first(clade_analysis.iqtree_version),
    input_14 = select_first(clade_analysis.snp_dist_version)
  }
  output {
    File ska_aln = ska.ska_aln
    String gubbins_date = gubbins_init.date
    String ska_docker = ska.ska_docker_image

    File? gubbins_polymorph_site_fasta = gubbins_init.polymorph_site_fasta
    File? gubbins_polymorph_site_phylip = gubbins_init.polymorph_site_phylip
    File? gubbins_branch_stats = gubbins_init.branch_stats
    File? gubbins_recomb_gff = gubbins_init.recomb_gff
    File? gubbins_snps= gubbins_init.gubbins_snps

    File? masked_fastas = mask_gubbins_init.masked_fastas
    #Array[File?] masked_fasta_list = mask_gubbins_init.masked_fasta_list

    File? tree = total_iqtree.ml_tree

    File clade_list_file = split_by_clade.clade_list_file

    Array[File?] gubbins_clade_polymorph_fasta = select_all(clade_analysis.gubbins_clade_polymorph_fasta)
    Array[File?] gubbins_clade_branch_stats = select_all(clade_analysis.gubbins_clade_branch_stats)
    Array[File?] gubbins_clade_recomb_gff = select_all(clade_analysis.gubbins_clade_recomb_gff)

    Array[File] pirate_pangenome_summary = select_all(clade_analysis.pirate_pangenome_summary)
    Array[File] pirate_gene_families_ordered = select_all(clade_analysis.pirate_gene_families_ordered)
    Array[String] pirate_docker_image = select_all(clade_analysis.pirate_docker_image)
    Array[File] pirate_gene_presence_absence = select_all(clade_analysis.pirate_for_scoary_csv)
    Array[File] pirate_aln_pan = clade_analysis.pirate_aln_pan
    Array[File] pirate_aln_core = clade_analysis.pirate_aln_core
    # snp_dists outputs
    Array[String?] clade_snps_dists_version = select_all(clade_analysis.clade_snps_dists_version)
    Array[File?] clade_core_snp_matrix = select_all(clade_analysis.clade_core_snp_matrix)
    Array[File?] clade_pan_snp_matrix = select_all(clade_analysis.clade_pan_snp_matrix)
    # iqtree outputs
    Array[String?] clade_iqtree_version = select_all(clade_analysis.clade_iqtree_version)
    Array[File?] clade_iqtree_core_tree = select_all(clade_analysis.clade_iqtree_core_tree)
    Array[File?] clade_iqtree_pan_tree = select_all(clade_analysis.clade_iqtree_pan_tree)
    Array[String?] clade_iqtree_pan_model = select_all(clade_analysis.clade_iqtree_pan_model)
    Array[String?] clade_iqtree_core_model = select_all(clade_analysis.clade_iqtree_core_model)
    Array[File?] plot_roary = clade_analysis.plot_roary
    File tool_versions = version.input_file
    File zipped_output = zip_files.zipped_output
  }
}
