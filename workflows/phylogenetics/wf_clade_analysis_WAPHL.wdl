version 1.0

import "../tasks/phylogenetic_inference/task_ska.wdl" as ska
import "../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../tasks/phylogenetic_inference/task_gubbins.wdl" as gubbins
import "../tasks/phylogenetic_inference/task_pirate.wdl" as pirate
import "../tasks/phylogenetic_inference/task_snp_dists.wdl" as snp_dists
import "../tasks/task_versioning.wdl" as versioning
import "../tasks/utilities/task_utilities.wdl" as utilities
import "../tasks/gene_typing/task_prokka.wdl" as prokka
import "../tasks/phylogenetic_inference/task_ksnp3.wdl" as ksnp3
import "../tasks/utilities/task_summarize_table_waphl.wdl" as summarize
import "../tasks/utilities/task_report_waphl.wdl" as report


workflow clade_analysis {
  input {

    Array[File] prokka_gff
    Array[String] samplename
    String iqtree_model = "MFP"
    Boolean? core = true
    Boolean? pan = false
    String cluster_name
    Float filter_perc = 25.0
    Boolean summarize = true

  }

call pirate.pirate as pirate {
  input:
    prokka_gff = prokka_gff,
    cluster_name = cluster_name
  }

call gubbins.gubbins as gubbins_clade {
  input:
    alignment = pirate.pirate_pangenome_alignment_fasta,
    filter_perc = filter_perc,
    cluster_name = cluster_name
}

if (pan == true) {
  if (gubbins_clade.gubbins_mask == true) {
    call gubbins.maskrc_svg as pan_mask_gubbins_clade  {
      input:
        alignment = pirate.pirate_pangenome_alignment_fasta,
        cluster_name = cluster_name,
        recomb = gubbins_clade.recomb_gff,
        base_reconstruct = gubbins_clade.base_reconstruct,
        recomb_embl = gubbins_clade.recomb_embl,
        polymorph_site_fasta = gubbins_clade.polymorph_site_fasta,
        polymorph_site_phylip = gubbins_clade.polymorph_site_phylip,
        branch_stats = gubbins_clade.branch_stats,
        gubbins_snps = gubbins_clade.gubbins_snps,
        gubbins_final_tre = gubbins_clade.gubbins_final_tre,
        gubbins_log = gubbins_clade.gubbins_log,
        gubbins_node_tre = gubbins_clade.gubbins_node_tre
    }
    call iqtree.iqtree as masked_pan_iqtree {
      input:
        alignment = pan_mask_gubbins_clade.masked_aln,
        cluster_name = cluster_name,
        iqtree_model = iqtree_model,
        details = "masked_pan_"
    }
}
  if (gubbins_clade.gubbins_mask == false) {
    call iqtree.iqtree as unmasked_pan_iqtree {
      input:
        alignment = pirate.pirate_pangenome_alignment_fasta,
        cluster_name = cluster_name,
        iqtree_model = iqtree_model,
        details = "unmasked_pan_"
    }
  }
  call snp_dists.snp_dists as pan_snp_dists {
    input:
      alignment = select_first([pan_mask_gubbins_clade.masked_aln, pirate.pirate_pangenome_alignment_fasta]),
      cluster_name = cluster_name
  }
}
  if (core == true) {
    if (gubbins_clade.gubbins_mask == true) {
      call gubbins.maskrc_svg as core_mask_gubbins_clade  {
        input:
          alignment = pirate.pirate_pangenome_alignment_fasta,
          cluster_name = cluster_name,
          recomb = gubbins_clade.recomb_gff,
          base_reconstruct = gubbins_clade.base_reconstruct,
          recomb_embl = gubbins_clade.recomb_embl,
          polymorph_site_fasta = gubbins_clade.polymorph_site_fasta,
          polymorph_site_phylip = gubbins_clade.polymorph_site_phylip,
          branch_stats = gubbins_clade.branch_stats,
          gubbins_snps = gubbins_clade.gubbins_snps,
          gubbins_final_tre = gubbins_clade.gubbins_final_tre,
          gubbins_log = gubbins_clade.gubbins_log,
          gubbins_node_tre = gubbins_clade.gubbins_node_tre
      }
      call ksnp3.ksnp3 as ksnp3_clade_core {
        input:
          assembly_fasta = core_mask_gubbins_clade.masked_fasta_list,
          samplename = samplename,
          cluster_name = cluster_name
      }
      call iqtree.iqtree as masked_core_iqtree {
        input:
          alignment =ksnp3_clade_core.ksnp3_core_matrix,
          cluster_name = cluster_name,
          iqtree_model = iqtree_model,
          details = "masked_core_"
      }
    }
    if (gubbins_clade.gubbins_mask == false) {
      call iqtree.iqtree as unmasked_core_iqtree {
        input:
          alignment =pirate.pirate_core_alignment_fasta,
          cluster_name = cluster_name,
          iqtree_model = iqtree_model,
          details = "unmasked_core_"
      }
    }
    call snp_dists.snp_dists as core_snp_dists {
      input:
        alignment = select_first([ksnp3_clade_core.ksnp3_core_matrix,pirate.pirate_core_alignment_fasta]),
        cluster_name = cluster_name
    }
  call utilities.generate_none {
    input:
    }
  if (summarize == true) {
  call summarize.zip_files as zip_files  {
  input:
    recomb_gff = select_all([gubbins_clade.recomb_gff]),
    pirate_aln_gff = select_all([pirate.pirate_pangenome_alignment_gff]),
    pirate_gene_presence_absence = select_all([pirate.pirate_for_scoary_csv]),
    cluster_name = cluster_name,
    cluster_tree = select_first([masked_pan_iqtree.ml_tree, unmasked_pan_iqtree.ml_tree, masked_core_iqtree.ml_tree, unmasked_core_iqtree.ml_tree]),
    #terra_table = terra_table,
    #terra_workspace = terra_workspace,
    #terra_project = terra_project,
    
}
}
  call report.plot_roary_waphl as plot_roary  {
  input:
    recomb_gff = gubbins_clade.recomb_gff,
    pirate_aln_gff = pirate.pirate_pangenome_alignment_gff,
    pirate_for_scoary_csv = pirate.pirate_for_scoary_csv,
    cluster_name = cluster_name,
    treefile = select_first([masked_pan_iqtree.ml_tree, unmasked_pan_iqtree.ml_tree, masked_core_iqtree.ml_tree, unmasked_core_iqtree.ml_tree])
    #terra_table = terra_table,
    #terra_workspace = terra_workspace,
    #terra_project = terra_project,
    
}
  call versioning.waphl_version_capture as version {
    input:
      input_1 = pirate.pirate_docker_image,
      input_2 = gubbins_clade.gubbins_docker_image,
      input_3 = select_first([core_mask_gubbins_clade.maskrc_docker_image, pan_mask_gubbins_clade.maskrc_docker_image, generate_none.none_string]),
      input_4 = ksnp3_clade_core.ksnp3_docker_image,
      input_5 = select_first([masked_pan_iqtree.version, unmasked_pan_iqtree.version, masked_core_iqtree.version, unmasked_core_iqtree.version]),
      input_6 = select_first([pan_snp_dists.version, core_snp_dists.version]),
      input_7 = plot_roary.plot_roary_docker_image
    }
  }

  output {

    String gubbins_date = gubbins_clade.date
    File? gubbins_clade_polymorph_fasta = gubbins_clade.polymorph_site_fasta
    File? gubbins_clade_branch_stats = gubbins_clade.branch_stats
    File? gubbins_clade_recomb_gff = gubbins_clade.recomb_gff
    Boolean clade_recombinants_detected = gubbins_clade.gubbins_mask
    File? gubbins_clade_nonrecomb_vcf = gubbins_clade.gubbins_nonrecomb_vcf
    File? gubbins_clade_snps_vcf = gubbins_clade.gubbins_snps

    String pirate_docker_image = pirate.pirate_docker_image
    String gubbins_docker_image = gubbins_clade.gubbins_docker_image
    String? maskrc_docker_image = select_first([core_mask_gubbins_clade.maskrc_docker_image, pan_mask_gubbins_clade.maskrc_docker_image, generate_none.none_string])
    String? ksnp3_docker_image = ksnp3_clade_core.ksnp3_docker_image
    String? iqtree_version = select_first([masked_pan_iqtree.version, unmasked_pan_iqtree.version, masked_core_iqtree.version, unmasked_core_iqtree.version])
    String? snp_dist_version = select_first([pan_snp_dists.version, core_snp_dists.version])

    File pirate_pangenome_summary = pirate.pirate_pangenome_summary
    File pirate_aln_pan = pirate.pirate_pangenome_alignment_gff
    File pirate_aln_core = pirate.pirate_core_alignment_gff
    File pirate_gene_families_ordered = pirate.pirate_gene_families_ordered
    String pirate_for_scoary_csv = pirate.pirate_for_scoary_csv
    # snp_dists outputs
    String? clade_snps_dists_version = select_first([core_snp_dists.version, pan_snp_dists.version])#core_snp_dists.version
    File? clade_core_snp_matrix = core_snp_dists.snp_matrix
    File? clade_pan_snp_matrix = pan_snp_dists.snp_matrix
    # iqtree outputs
    String? clade_iqtree_version = select_first([unmasked_pan_iqtree.version, masked_pan_iqtree.version, unmasked_core_iqtree.version, masked_core_iqtree.version])#pan_iqtree.version
    File? clade_iqtree_core_tree = select_first([masked_core_iqtree.ml_tree, unmasked_core_iqtree.ml_tree])#core_iqtree.ml_tree
    File? clade_iqtree_pan_tree = select_first([masked_pan_iqtree.ml_tree, unmasked_pan_iqtree.ml_tree])#pan_iqtree.ml_tree
    String? clade_iqtree_pan_model = select_first([masked_pan_iqtree.iqtree_model_used, unmasked_pan_iqtree.iqtree_model_used])#pan_iqtree.iqtree_model
    String? clade_iqtree_core_model = select_first([masked_core_iqtree.iqtree_model_used, unmasked_core_iqtree.iqtree_model_used])#core_iqtree.iqtree_model
    File? software_versions_clade_analysis = version.tool_versions
    File? clade_zipped_output = zip_files.zipped_output
    File? plot_roary = plot_roary.plot_roary_png

  }
}
