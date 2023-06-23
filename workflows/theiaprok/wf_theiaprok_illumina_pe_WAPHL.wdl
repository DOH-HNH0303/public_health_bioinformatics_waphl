version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../utilities/wf_merlin_magic.wdl" as merlin_magic
import "../../tasks/assembly/task_shovill.wdl" as shovill
import "../../tasks/quality_control/task_quast.wdl" as quast
import "../../tasks/quality_control/task_cg_pipeline.wdl" as cg_pipeline
import "../../tasks/taxon_id/task_gambit.wdl" as gambit
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/species_typing/task_serotypefinder.wdl" as serotypefinder
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/task_broad_terra_tools.wdl" as terra_tools
import "../../tasks/quality_control/task_general_qc.wdl" as general_qc

#Theiagen packages added
import "../../tasks/quality_control/task_screen.wdl" as screen
import "../../tasks/quality_control/task_busco.wdl" as busco
import "../../tasks/quality_control/task_mummer_ani.wdl" as ani
import "../../tasks/gene_typing/task_amrfinderplus.wdl" as amrfinderplus
import "../../tasks/gene_typing/task_resfinder.wdl" as resfinder
import "../../tasks/species_typing/task_ts_mlst.wdl" as ts_mlst
import "../../tasks/gene_typing/task_prokka.wdl" as prokka
import "../../tasks/gene_typing/task_plasmidfinder.wdl" as plasmidfinder
import "../../tasks/utilities/task_utilities.wdl" as utilities
import "../../tasks/task_qc_utils.wdl" as qc
import "../../tasks/task_taxon_id.wdl" as taxon_id
import "../../tasks/taxon_id/task_fastani.wdl" as fastani
import "../../tasks/task_denovo_assembly.wdl" as assembly
import "../../tasks/task_read_clean.wdl" as read_clean
import "../../tasks/utilities/task_summarize_table_waphl.wdl" as summarize


workflow theiaprok_illumina_pe {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of paired-end bacterial NGS data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    String? run_id
    File? taxon_table
    String terra_table
    String terra_project
    String terra_workspace

    Int? genome_size
    String? collection_date
    String? originating_lab
    String? city
    String? county
    String? zip
    # by default do not call ANI task, but user has ability to enable this task if working with enteric pathogens or supply their own high-quality reference genome
    Boolean call_ani = false
    Int min_reads = 7472
    Int min_basepairs = 2241820
    Int min_genome_size = 100000
    Int max_genome_size = 18040666
    Int min_coverage = 10
    Int min_proportion = 50
    Boolean call_resfinder = false
    Boolean skip_screen = false
    String wget_dt_docker_image = "inutano/wget:1.20.3-r1"
  }

  call versioning.version_capture{
    input:
  }

  call screen.check_reads as raw_check_reads {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_size = min_genome_size,
      max_genome_size = max_genome_size,
      min_coverage = min_coverage,
      min_proportion = min_proportion,
      skip_screen = skip_screen
  }
  call taxon_id.kraken2 as kraken2_raw {
    input:
    samplename = samplename,
    read1 = read1_raw,
    read2 = read2_raw
  }
  call read_clean.ncbi_scrub_pe {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw
  }
  if (raw_check_reads.read_screen=="PASS") {
    call read_qc.read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = ncbi_scrub_pe.read1_dehosted,
      read2_raw = ncbi_scrub_pe.read2_dehosted
    }

    call taxon_id.kraken2 as kraken2_clean {
      input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean
    }
    call screen.check_reads as clean_check_reads {
      input:
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_size = min_genome_size,
        max_genome_size = max_genome_size,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        skip_screen = skip_screen
    }
    if (clean_check_reads.read_screen=="PASS") {
      call shovill.shovill_pe {
        input:
          samplename = samplename,
          read1_cleaned = read_QC_trim.read1_clean,
          read2_cleaned = read_QC_trim.read2_clean,
          genome_size = select_first([genome_size, clean_check_reads.est_genome_length])
      }
      call quast.quast {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call busco.busco {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call general_qc.general_qc {
        input:
          assembly_fasta = shovill_pe.assembly_fasta
      }
      call cg_pipeline.cg_pipeline {
        input:
          read1 = read1_raw,
          read2 = read2_raw,
          samplename = samplename,
          genome_length = select_first([genome_size, clean_check_reads.est_genome_length])
      }
      call gambit.gambit {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }

      if (call_ani) {
      call ani.animummer as ani {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      }
      call amrfinderplus.amrfinderplus_nuc as amrfinderplus_task {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename,
          organism = gambit.gambit_predicted_taxon
      }
      if (call_resfinder) {
      call resfinder.resfinder as resfinder_task {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename,
          organism = gambit.gambit_predicted_taxon
        }
      }
      call ts_mlst.ts_mlst {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call prokka.prokka {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call plasmidfinder.plasmidfinder {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
  #call shovill.shovill_pe {
  #  input:
  #    samplename = samplename,
  #    read1_cleaned = read_QC_trim.read1_clean,
  #    read2_cleaned = read_QC_trim.read2_clean
  #}
  if (kraken2_clean.kraken2_genus=="Corynebacterium" || kraken2_clean.kraken2_species=="diphtheriae"){
    call taxon_id.ncbi_blast {
      input:
        samplename=samplename,
        assembly=shovill_pe.assembly_fasta
  }
  call utilities.get_dt_results {
    input:
      samplename=samplename,
      tblastn_dt_omega_report=ncbi_blast.tblastn_dt_omega_report,
      tblastn_dt_beta_report=ncbi_blast.tblastn_dt_beta_report,
      tblastn_dt_beta_homologue_report=ncbi_blast.tblastn_dt_beta_homologue_report
}
}

  if (kraken2_clean.kraken2_genus=="Legionella" || kraken2_clean.kraken2_genus=="Tatlockia" ||kraken2_clean.kraken2_genus=="Corynebacterium" || kraken2_clean.kraken2_genus=="Fluoribacter"){
    call fastani.fastANI {
      input:
        samplename = samplename,
        assembly = shovill_pe.assembly_fasta,
        genus=kraken2_clean.kraken2_genus
    }
    call utilities.join_genus_species {
      input:
        genus = fastANI.fastani_genus,
        species = fastANI.fastani_species
    }
    call merlin_magic.merlin_magic as merlin_magic_cdip_legionella {
      input:
        merlin_tag = join_genus_species.genus_species,
        assembly = shovill_pe.assembly_fasta,
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean
    }
}
  if  (kraken2_clean.kraken2_genus!="Legionella" || kraken2_clean.kraken2_genus!="Tatlockia" ||kraken2_clean.kraken2_genus!="Corynebacterium" || kraken2_clean.kraken2_genus!="Fluoribacter"){
    call merlin_magic.merlin_magic {
      input:
        merlin_tag = gambit.merlin_tag,
        assembly = shovill_pe.assembly_fasta,
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean
    }
  }

  call abricate.abricate as abricate_amr {
    input:
      assembly = shovill_pe.assembly_fasta,
      samplename = samplename,
      database = "ncbi",
      gene_type = "AMR"
  }
  call abricate.abricate as abricate_virulence {
    input:
      assembly = shovill_pe.assembly_fasta,
      samplename = samplename,
      database = "vfdb"
  }
  }
  }
  output {
    #Version Captures
    String theiaprok_illumina_pe_version = version_capture.phbg_version
    String theiaprok_illumina_pe_analysis_date = version_capture.date
    #Read Metadata
    String seq_platform = seq_method
    #Sample Screening
    String raw_read_screen = raw_check_reads.read_screen
    String? clean_read_screen = clean_check_reads.read_screen
    #Read QC read_QC_trim.read1_clean
    Int? num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    File? reads_clean1 = read_QC_trim.read1_clean
    File? reads_clean2 = read_QC_trim.read2_clean
    String? num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_software = read_QC_trim.trimmomatic_pe_software
    String? bbduk_docker = read_QC_trim.bbduk_docker
    Float? r1_mean_q = cg_pipeline.r1_mean_q
    Float? r2_mean_q = cg_pipeline.r2_mean_q

    #Assembly and Assembly QC
    File? assembly_fasta = shovill_pe.assembly_fasta
    File? contigs_gfa = shovill_pe.contigs_gfa
    File? contigs_fastg = shovill_pe.contigs_fastg
    File? contigs_lastgraph = shovill_pe.contigs_lastgraph
    String? shovill_pe_version = shovill_pe.shovill_version
    File? quast_report = quast.quast_report
    String? quast_version = quast.version
    Int? genome_length = quast.genome_length
    Int? number_contigs = quast.number_contigs
    Int? n50_value = quast.n50_value
    Float? gc_content = quast.gc_content
    File? cg_pipeline_report = cg_pipeline.cg_pipeline_report
    String? cg_pipeline_docker = cg_pipeline.cg_pipeline_docker
    Float? est_coverage = cg_pipeline.est_coverage
    String? busco_version = busco.busco_version
    String? busco_database = busco.busco_database
    String? busco_results = busco.busco_results
    File? busco_report = busco.busco_report
    Int? number_N = general_qc.number_N
    Int? number_ATCG = general_qc.number_ATCG
    Int? number_Degenerate = general_qc.number_Degenerate
    Int? number_Total = general_qc.number_Total
    #Taxon ID
    File? gambit_report = gambit.gambit_report_file
    File? gabmit_closest_genomes = gambit.gambit_closest_genomes_file
    String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String? gambit_version = gambit.gambit_version
    String? gambit_db_version = gambit.gambit_db_version
    String? gambit_docker = gambit.gambit_docker

    String  kraken2_raw_version              = kraken2_raw.version
    Float   kraken2_raw_human                = kraken2_raw.percent_human
    String  kraken2_raw_report               = kraken2_raw.kraken_report
    String  kraken2_raw_genus              = kraken2_raw.kraken2_genus
    String   kraken2_raw_species                = kraken2_raw.kraken2_species
    String  kraken2_raw_strain               = kraken2_raw.kraken2_strain

    String?  kraken2_clean_version              = kraken2_clean.version
    Float?   kraken2_clean_human                = kraken2_clean.percent_human
    String?  kraken2_clean_report               = kraken2_clean.kraken_report
    String?  kraken2_clean_genus              = kraken2_clean.kraken2_genus
    String?   kraken2_clean_species                = kraken2_clean.kraken2_species
    String?  kraken2_clean_strain               = kraken2_clean.kraken2_strain

    File?    fastani_report   =fastANI.fastani_report
    String?    fastani_genus   =fastANI.fastani_genus
    String?    fastani_species   =fastANI.fastani_species
    String?    fastani_strain   =fastANI.fastani_strain
    Float?    fastani_ani_estimate   =fastANI.fastani_aniestimate
    String?   fastani_genus_species = join_genus_species.genus_species

    File?    dt_omega_tsv=ncbi_blast.tblastn_dt_omega_report
    File?    dt_beta_tsv=ncbi_blast.tblastn_dt_beta_report
    File?    dt_beta_homologue_tsv=ncbi_blast.tblastn_dt_beta_homologue_report
    String?   dt_omega=get_dt_results.dt_omega
    String? dt_beta =get_dt_results.dt_beta
    String? dt_beta_homologue =get_dt_results.dt_beta_homologue
    String? dt_omega_evalue =get_dt_results.dt_omega_evalue
    String? dt_beta_evalue =get_dt_results.dt_beta_evalue
    String? dt_beta_homologue_evalue =get_dt_results.dt_beta_homologue_evalue



    #String gambit_closest_taxon = gambit.gambit_closest_match
    #Midas taxonomy
    #File midas_report = midas.midas_report
    #String midas_predicted_genus = midas.midas_genus
    #String midas_predicted_species = midas.midas_species
    #String midas_predicted_strain = midas.midas_strain
    #String midas_docker = midas.midas_docker

    # ani-mummer
    Float? ani_highest_percent = ani.ani_highest_percent
    Float? ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned
    File? ani_output_tsv = ani.ani_output_tsv
    String? ani_top_species_match = ani.ani_top_species_match
    String? ani_mummer_version = ani.ani_mummer_version
    # NCBI-AMRFinderPlus Outputs
    File? amrfinderplus_all_report = amrfinderplus_task.amrfinderplus_all_report
    File? amrfinderplus_amr_report = amrfinderplus_task.amrfinderplus_amr_report
    File? amrfinderplus_stress_report = amrfinderplus_task.amrfinderplus_stress_report
    File? amrfinderplus_virulence_report = amrfinderplus_task.amrfinderplus_virulence_report
    String? amrfinderplus_amr_genes = amrfinderplus_task.amrfinderplus_amr_genes
    String? amrfinderplus_stress_genes = amrfinderplus_task.amrfinderplus_stress_genes
    String? amrfinderplus_virulence_genes = amrfinderplus_task.amrfinderplus_virulence_genes
    String? amrfinderplus_version = amrfinderplus_task.amrfinderplus_version
    String? amrfinderplus_db_version = amrfinderplus_task.amrfinderplus_db_version

    #AMR Screening
    File? abricate_amr_results = abricate_amr.abricate_results
    String? abricate_amr_database = abricate_amr.abricate_database
    String? abricate_amr_version = abricate_amr.abricate_version
    String? abricate_amr_genes = abricate_amr.abricate_genes


    # Resfinder Outputs
    File? resfinder_pheno_table = resfinder_task.resfinder_pheno_table
    File? resfinder_pheno_table_species = resfinder_task.resfinder_pheno_table_species
    File? resfinder_seqs = resfinder_task.resfinder_hit_in_genome_seq
    File? resfinder_results = resfinder_task.resfinder_results_tab
    File? resfinder_pointfinder_pheno_table = resfinder_task.pointfinder_pheno_table
    File? resfinder_pointfinder_results = resfinder_task.pointfinder_results
    String? resfinder_db_version = resfinder_task.resfinder_db_version
    String? resfinder_docker = resfinder_task.resfinder_docker

    # Virulence Genes
    File? abricate_virulence_results = abricate_virulence.abricate_results
    String? abricate_virulence_database = abricate_virulence.abricate_database
    String? abricate_virulence_version = abricate_virulence.abricate_version
    String? abricate_virulence_genes = abricate_virulence.abricate_genes
    # MLST Typing
    File? ts_mlst_results = ts_mlst.ts_mlst_results
    String? ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st
    String? ts_mlst_version = ts_mlst.ts_mlst_version
    String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
    # Prokka Results
    File? prokka_gff = prokka.prokka_gff
    File? prokka_gbk = prokka.prokka_gbk
    File? prokka_sqn = prokka.prokka_sqn

    # Plasmidfinder Results
    String? plasmidfinder_plasmids = plasmidfinder.plasmidfinder_plasmids
    File? plasmidfinder_results = plasmidfinder.plasmidfinder_results
    File? plasmidfinder_seqs = plasmidfinder.plasmidfinder_seqs
    String? plasmidfinder_docker = plasmidfinder.plasmidfinder_docker
    String? plasmidfinder_db_version = plasmidfinder.plasmidfinder_db_version
    # Ecoli Typing
    File? serotypefinder_report = merlin_magic.serotypefinder_report
    String? serotypefinder_docker = merlin_magic.serotypefinder_docker
    String? serotypefinder_serotype = merlin_magic.serotypefinder_serotype
    File? ectyper_results = merlin_magic.ectyper_results
    String? ectyper_version = merlin_magic.ectyper_version
    String? ectyper_predicted_serotype = merlin_magic.ectyper_predicted_serotype
    #Listeria Typing
    File? lissero_results = merlin_magic.lissero_results
    String? lissero_version = merlin_magic.lissero_version
    String? lissero_serotype = merlin_magic.lissero_serotype
    #Salmonella Typing
    File? sistr_results = merlin_magic.sistr_results
    File? sistr_allele_json = merlin_magic.sistr_allele_json
    File? sister_allele_fasta = merlin_magic.sistr_allele_fasta
    File? sistr_cgmlst = merlin_magic.sistr_cgmlst
    String? sistr_version = merlin_magic.sistr_version
    String? sistr_predicted_serotype = merlin_magic.sistr_predicted_serotype
    File? seqsero2_report = merlin_magic.seqsero2_report
    String? seqsero2_version = merlin_magic.seqsero2_version
    String? seqsero2_predicted_antigenic_profile = merlin_magic.seqsero2_predicted_antigenic_profile
    String? seqsero2_predicted_serotype = merlin_magic.seqsero2_predicted_serotype
    String? seqsero2_predicted_contamination = merlin_magic.seqsero2_predicted_contamination
    # Salmonella serotype Typhi Typing
    File? genotyphi_report_tsv = merlin_magic.genotyphi_report_tsv
    File? genotyphi_mykrobe_json = merlin_magic.genotyphi_mykrobe_json
    String? genotyphi_version = merlin_magic.genotyphi_version
    String? genotyphi_species = merlin_magic.genotyphi_species
    Float? genotyphi_st_probes_percent_coverage = merlin_magic.genotyphi_st_probes_percent_coverage
    String? genotyphi_final_genotype = merlin_magic.genotyphi_final_genotype
    String? genotyphi_genotype_confidence = merlin_magic.genotyphi_genotype_confidence
    #Klebsiella Typing
    File? kleborate_output_file = merlin_magic.kleborate_output_file
    String? kleborate_version = merlin_magic.kleborate_version
    String? kleborate_key_resistance_genes = merlin_magic.kleborate_key_resistance_genes
    String? kleborate_genomic_resistance_mutations = merlin_magic.kleborate_genomic_resistance_mutations
    String? kleborate_mlst_sequence_type = merlin_magic.kleborate_mlst_sequence_type
    # Legionella pneumophila typing
    File? legsta_results = merlin_magic_cdip_legionella.legsta_results
    String? legsta_predicted_sbt = merlin_magic_cdip_legionella.legsta_predicted_sbt
    String? legsta_version = merlin_magic_cdip_legionella.legsta_version
    # Mycobacterium Typing
    File? tbprofiler_output_file = merlin_magic.tbprofiler_output_file
    File? tbprofiler_output_bam = merlin_magic.tbprofiler_output_bam
    File? tbprofiler_output_bai = merlin_magic.tbprofiler_output_bai
    String? tbprofiler_version = merlin_magic.tbprofiler_version
    String? tbprofiler_main_lineage = merlin_magic.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = merlin_magic.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = merlin_magic.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = merlin_magic.tbprofiler_resistance_genes
  }
}
