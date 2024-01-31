version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../../tasks/assembly/task_shovill.wdl" as shovill
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/quality_control/task_cg_pipeline.wdl" as cg_pipeline
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/species_typing/task_serotypefinder.wdl" as serotypefinder
import "../../tasks/quality_control/task_general_qc.wdl" as general_qc

#Theiagen packages added
import "../../tasks/quality_control/task_screen.wdl" as screen
import "../../tasks/quality_control/task_busco.wdl" as busco_task
import "../../tasks/taxon_id/task_gambit.wdl" as gambit_task

import "../../tasks/quality_control/task_mummer_ani.wdl" as ani_task
import "../../tasks/taxon_id/task_kmerfinder.wdl" as kmerfinder_task
#import "../../tasks/quality_control/task_qc_utils.wdl" as qc
import "../../tasks/gene_typing/task_amrfinderplus.wdl" as amrfinderplus
import "../../tasks/gene_typing/task_resfinder.wdl" as resfinder
import "../../tasks/species_typing/task_ts_mlst.wdl" as ts_mlst_task
import "../../tasks/gene_typing/task_bakta.wdl" as bakta_task
import "../../tasks/gene_typing/task_prokka.wdl" as prokka_task
import "../../tasks/gene_typing/task_plasmidfinder.wdl" as plasmidfinder_task
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/task_broad_terra_tools.wdl" as terra_tools
import "../../tasks/utilities/task_utilities.wdl" as utilities

import "../../tasks/taxon_id/task_taxon_id_waphl.wdl" as taxon_id
import "../../tasks/taxon_id/task_fastani.wdl" as fastani
#import "../../tasks/task_denovo_assembly.wdl" as assembly
#import "../../tasks/task_read_clean.wdl" as read_clean
#import "../../tasks/utilities/task_summarize_table_waphl.wdl" as summarize


workflow theiaprok_illumina_pe_waphl {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of paired-end bacterial NGS data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    Int? genome_size
    # export taxon table parameters
    String? run_id
    String? collection_date
    String? originating_lab
    String? city
    String? county
    String? zip
    File? taxon_tables
    #String terra_table
    String terra_project = "NA"
    String terra_workspace = "NA"
    # read screen parameters
    Boolean skip_screen = false
    Int min_reads = 7472
    Int min_basepairs = 2241820
    Int min_genome_size = 100000
    Int max_genome_size = 18040666
    Int min_coverage = 10
    Int min_proportion = 40
        # trimming parameters
    Int trim_minlen = 75
    Int trim_quality_trim_score = 20
    Int trim_window_size = 10
     # module options
    Boolean call_ani = false # by default do not call ANI task, but user has ability to enable this task if working with enteric pathogens or supply their own high-quality reference genome
    Boolean call_kmerfinder = false 
    Boolean call_resfinder = false
    String wget_dt_docker_image = "inutano/wget:1.20.3-r1"
    String genome_annotation = "prokka" # options: "prokka" or "bakta"
    # qc check parameters
    String? expected_taxon  # allow user to provide organism (e.g. "Clostridioides_difficile") string to amrfinder. Useful when gambit does not predict the correct species    # qc check parameters
    File? qc_check_table
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
      skip_screen = skip_screen,
      expected_genome_size = genome_size
  }
  call taxon_id.kraken2 as kraken2_raw {
    input:
    samplename = samplename,
    read1 = read1_raw,
    read2 = read2_raw
  }
  if (raw_check_reads.read_screen=="PASS") {
    call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      samplename = samplename,
        read1_raw = read1_raw,
        read2_raw = read2_raw,
        trim_minlen = trim_minlen,
        trim_quality_trim_score = trim_quality_trim_score,
        trim_window_size = trim_window_size,
        workflow_series = "theiaprok"
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
        skip_screen = skip_screen,
        expected_genome_size = genome_size
    }
    if (clean_check_reads.read_screen == "PASS") {
      call shovill.shovill_pe {
        input:
          samplename = samplename,
          read1_cleaned = read_QC_trim.read1_clean,
          read2_cleaned = read_QC_trim.read2_clean,
          genome_size = select_first([genome_size, clean_check_reads.est_genome_length])
      }
      call quast_task.quast {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call general_qc.general_qc {
        input:
          assembly_fasta = shovill_pe.assembly_fasta
      }
      call cg_pipeline.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = read1_raw,
          read2 = read2_raw,
          samplename = samplename,
          genome_length = select_first([genome_size, quast.genome_length])
      }
      call cg_pipeline.cg_pipeline as cg_pipeline_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          samplename = samplename,
          genome_length = select_first([genome_size, quast.genome_length])
      }
      call gambit_task.gambit {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call busco_task.busco {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      if (call_ani) {
      call ani_task.animummer as ani {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      }
      if (call_kmerfinder) {
        call kmerfinder_task.kmerfinder_bacteria as kmerfinder {
          input:
            assembly = shovill_pe.assembly_fasta,
            samplename = samplename
        }
      }
      call amrfinderplus.amrfinderplus_nuc as amrfinderplus_task {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename,
          organism = select_first([expected_taxon, gambit.gambit_predicted_taxon])
      }
      if (call_resfinder) {
      call resfinder.resfinder as resfinder_task {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename,
          organism = select_first([expected_taxon, gambit.gambit_predicted_taxon])
        }
      }
      call ts_mlst_task.ts_mlst {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      if (genome_annotation == "prokka") {
        call prokka_task.prokka {
          input:
            assembly = shovill_pe.assembly_fasta,
            samplename = samplename
        }
      }
      if (genome_annotation == "bakta") {
        call bakta_task.bakta {
          input:
            assembly = shovill_pe.assembly_fasta,
            samplename = samplename
        }
      }
      call plasmidfinder_task.plasmidfinder {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      # if(defined(qc_check_table)) {
      #   call qc_check.qc_check_phb as qc_check_task {
      #     input:
      #       qc_check_table = qc_check_table,
      #       expected_taxon = expected_taxon,
      #       gambit_predicted_taxon = gambit.gambit_predicted_taxon,
      #       num_reads_raw1 = read_QC_trim.fastq_scan_raw1,
      #       num_reads_raw2 = read_QC_trim.fastq_scan_raw2,
      #       num_reads_clean1 = read_QC_trim.fastq_scan_clean1,
      #       num_reads_clean2 = read_QC_trim.fastq_scan_clean2,
      #       r1_mean_q_raw = cg_pipeline_raw.r1_mean_q,
      #       r2_mean_q_raw = cg_pipeline_raw.r2_mean_q,
      #       combined_mean_q_raw = cg_pipeline_raw.combined_mean_q,
      #       r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength,
      #       r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength,  
      #       combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength,
      #       r1_mean_q_clean = cg_pipeline_clean.r1_mean_q,
      #       r2_mean_q_clean = cg_pipeline_clean.r2_mean_q,
      #       combined_mean_q_clean = cg_pipeline_clean.combined_mean_q,
      #       r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength,
      #       r2_mean_readlength_clean = cg_pipeline_clean.r2_mean_readlength,  
      #       combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength,    
      #       est_coverage_raw = cg_pipeline_raw.est_coverage,
      #       est_coverage_clean = cg_pipeline_clean.est_coverage,
      #       midas_secondary_genus_abundance = read_QC_trim.midas_secondary_genus_abundance,
      #       assembly_length = quast.genome_length,
      #       number_contigs = quast.number_contigs,
      #       n50_value = quast.n50_value,
      #       quast_gc_percent = quast.gc_percent,
      #       busco_results = busco.busco_results,
      #       ani_highest_percent = ani.ani_highest_percent,
      #       ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned
      #   }
      # }
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
    call merlin_magic_workflow.merlin_magic as merlin_magic_cdip_legionella {
      input:
        merlin_tag = join_genus_species.genus_species,
        assembly = shovill_pe.assembly_fasta,
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean
    }
}
  if  (kraken2_clean.kraken2_genus!="Legionella" || kraken2_clean.kraken2_genus!="Tatlockia" ||kraken2_clean.kraken2_genus!="Corynebacterium" || kraken2_clean.kraken2_genus!="Fluoribacter"){
    call merlin_magic_workflow.merlin_magic {
      input:
        merlin_tag = select_first([expected_taxon, gambit.merlin_tag]),
        assembly = shovill_pe.assembly_fasta,
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean
    }
  }
  if(defined(qc_check_table)) {
        call qc_check.qc_check_phb_waphl as qc_check_task_waphl {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = expected_taxon,
            predicted_taxon = select_first([join_genus_species.genus_species, gambit.gambit_predicted_taxon, ""]),
            num_reads_raw1 = read_QC_trim.fastq_scan_raw1,
            num_reads_raw2 = read_QC_trim.fastq_scan_raw2,
            num_reads_clean1 = read_QC_trim.fastq_scan_clean1,
            num_reads_clean2 = read_QC_trim.fastq_scan_clean2,
            r1_mean_q_raw = cg_pipeline_raw.r1_mean_q,
            r2_mean_q_raw = cg_pipeline_raw.r2_mean_q,
            combined_mean_q_raw = cg_pipeline_raw.combined_mean_q,
            r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength,
            r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength,  
            combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength,
            r1_mean_q_clean = cg_pipeline_clean.r1_mean_q,
            r2_mean_q_clean = cg_pipeline_clean.r2_mean_q,
            combined_mean_q_clean = cg_pipeline_clean.combined_mean_q,
            r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength,
            r2_mean_readlength_clean = cg_pipeline_clean.r2_mean_readlength,  
            combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength,    
            est_coverage_raw = cg_pipeline_raw.est_coverage,
            est_coverage_clean = cg_pipeline_clean.est_coverage,
            midas_secondary_genus_abundance = read_QC_trim.midas_secondary_genus_abundance,
            assembly_length = quast.genome_length,
            number_contigs = quast.number_contigs,
            n50_value = quast.n50_value,
            quast_gc_percent = quast.gc_percent,
            busco_results = busco.busco_results,
            ani_highest_percent = ani.ani_highest_percent,
            ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned,
            number_N = general_qc.number_N,
            number_Total = general_qc.number_Total,
            kraken2_clean_human = kraken2_clean.percent_human
        }
  }
  call abricate.abricate as abricate_amr {
    input:
      assembly = shovill_pe.assembly_fasta,
      samplename = samplename,
      database = "ncbi"  #,
      #gene_type = "AMR"
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
    # Version Captures
    String theiaprok_illumina_pe_version = version_capture.phb_version
    String theiaprok_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String raw_read_screen = raw_check_reads.read_screen
    String? clean_read_screen = clean_check_reads.read_screen
    # Read QC - fastq_scan outputs
    Int? num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    String? fastqc_version = read_QC_trim.fastqc_version
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    #String? trimmomatic_software = read_QC_trim.trimmomatic_pe_software
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    # Read QC - bbduk outputs
    File? read1_clean = read_QC_trim.read1_clean
    File? read2_clean = read_QC_trim.read2_clean
    String? bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - cg pipeline outputs
    Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
    Float? r2_mean_q_raw = cg_pipeline_raw.r2_mean_q
    Float? combined_mean_q_raw = cg_pipeline_raw.combined_mean_q
    Float? combined_mean_q_clean = cg_pipeline_clean.combined_mean_q
    Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
    Float? r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength
    Float? combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength
    Float? combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength
    # Read QC - midas outputs
    String? midas_docker = read_QC_trim.midas_docker
    File? midas_report = read_QC_trim.midas_report
    String? midas_primary_genus = read_QC_trim.midas_primary_genus
    String? midas_secondary_genus = read_QC_trim.midas_secondary_genus
    Float? midas_secondary_genus_abundance = read_QC_trim.midas_secondary_genus_abundance
    # Assembly - shovill outputs 
    File? assembly_fasta = shovill_pe.assembly_fasta
    File? contigs_gfa = shovill_pe.contigs_gfa
    File? contigs_fastg = shovill_pe.contigs_fastg
    File? contigs_lastgraph = shovill_pe.contigs_lastgraph
    String? shovill_pe_version = shovill_pe.shovill_version
    # Assembly QC - quast outputs
    File? quast_report = quast.quast_report
    String? quast_version = quast.version
    Int? assembly_length = quast.genome_length
    Int? number_contigs = quast.number_contigs
    Int? n50_value = quast.n50_value
    Float? quast_gc_percent = quast.gc_percent
    # Assembly QC - cg pipeline outputs
    #File? cg_pipeline_report = cg_pipeline.cg_pipeline_report
    #String? cg_pipeline_docker = cg_pipeline.cg_pipeline_docker
    #Float? est_coverage = cg_pipeline.est_coverage
    # Assembly QC - cg pipeline outputs
    File? cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
    String? cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
    Float? est_coverage_raw = cg_pipeline_raw.est_coverage
    File? cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
    Float? est_coverage_clean = cg_pipeline_clean.est_coverage
    # Assembly QC - busco outputs
    String? busco_version = busco.busco_version
    String? busco_database = busco.busco_database
    String? busco_results = busco.busco_results
    File? busco_report = busco.busco_report
    Int? number_N = general_qc.number_N
    Int? number_ATCG = general_qc.number_ATCG
    Int? number_Degenerate = general_qc.number_Degenerate
    Int? number_Total = general_qc.number_Total
    # Taxon ID - gambit outputs
    File? gambit_report = gambit.gambit_report_file
    File? gambit_closest_genomes = gambit.gambit_closest_genomes_file
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

    # ani-mummer outputs
    Float? ani_highest_percent = ani.ani_highest_percent
    Float? ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned
    File? ani_output_tsv = ani.ani_output_tsv
    String? ani_top_species_match = ani.ani_top_species_match
    String? ani_mummer_version = ani.ani_mummer_version
    String? ani_mummer_docker = ani.ani_docker
    # kmerfinder outputs
    String? kmerfinder_docker = kmerfinder.kmerfinder_docker
    File? kmerfinder_results_tsv = kmerfinder.kmerfinder_results_tsv
    String? kmerfinder_top_hit = kmerfinder.kmerfinder_top_hit
    String? kmerfinder_query_coverage = kmerfinder.kmerfinder_query_coverage
    String? kmerfinder_template_coverage = kmerfinder.kmerfinder_template_coverage
    String? kmerfinder_database = kmerfinder.kmerfinder_database
    # NCBI-AMRFinderPlus Outputs
    File? amrfinderplus_all_report = amrfinderplus_task.amrfinderplus_all_report
    File? amrfinderplus_amr_report = amrfinderplus_task.amrfinderplus_amr_report
    File? amrfinderplus_stress_report = amrfinderplus_task.amrfinderplus_stress_report
    File? amrfinderplus_virulence_report = amrfinderplus_task.amrfinderplus_virulence_report
    #String? amrfinderplus_amr_genes = amrfinderplus_task.amrfinderplus_amr_genes
    String? amrfinderplus_amr_core_genes = amrfinderplus_task.amrfinderplus_amr_core_genes
    String? amrfinderplus_amr_plus_genes = amrfinderplus_task.amrfinderplus_amr_plus_genes
    String? amrfinderplus_stress_genes = amrfinderplus_task.amrfinderplus_stress_genes
    String? amrfinderplus_virulence_genes = amrfinderplus_task.amrfinderplus_virulence_genes
    String? amrfinderplus_amr_classes = amrfinderplus_task.amrfinderplus_amr_classes
    String? amrfinderplus_amr_subclasses = amrfinderplus_task.amrfinderplus_amr_subclasses
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
    String? resfinder_predicted_pheno_resistance = resfinder_task.resfinder_predicted_pheno_resistance
    String? resfinder_predicted_xdr_shigella = resfinder_task.resfinder_predicted_xdr_shigella
    String? resfinder_predicted_resistance_Amp = resfinder_task.resfinder_predicted_resistance_Amp
    String? resfinder_predicted_resistance_Azm = resfinder_task.resfinder_predicted_resistance_Azm
    String? resfinder_predicted_resistance_Axo = resfinder_task.resfinder_predicted_resistance_Axo
    String? resfinder_predicted_resistance_Cip = resfinder_task.resfinder_predicted_resistance_Cip
    String? resfinder_predicted_resistance_Smx = resfinder_task.resfinder_predicted_resistance_Smx
    String? resfinder_predicted_resistance_Tmp = resfinder_task.resfinder_predicted_resistance_Tmp
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
    String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
    String? ts_mlst_allelic_profile = ts_mlst.ts_mlst_allelic_profile
    String? ts_mlst_version = ts_mlst.ts_mlst_version
    File? ts_mlst_novel_alleles = ts_mlst.ts_mlst_novel_alleles
    String? ts_mlst_docker = ts_mlst.ts_mlst_docker
    # Prokka Results
    File? prokka_gff = prokka.prokka_gff
    File? prokka_gbk = prokka.prokka_gbk
    File? prokka_sqn = prokka.prokka_sqn
    # Bakta Results
    File? bakta_gbff = bakta.bakta_gbff
    File? bakta_gff3 = bakta.bakta_gff3
    File? bakta_tsv = bakta.bakta_tsv
    File? bakta_summary = bakta.bakta_txt
    String? bakta_version = bakta.bakta_version
    # Plasmidfinder Results
    String? plasmidfinder_plasmids = plasmidfinder.plasmidfinder_plasmids
    File? plasmidfinder_results = plasmidfinder.plasmidfinder_results
    File? plasmidfinder_seqs = plasmidfinder.plasmidfinder_seqs
    String? plasmidfinder_docker = plasmidfinder.plasmidfinder_docker
    String? plasmidfinder_db_version = plasmidfinder.plasmidfinder_db_version
    # QC_Check Results
    # String? qc_check = qc_check_task.qc_check
    # File? qc_standard = qc_check_task.qc_standard
    # QC_Check Results WAPHL
    String? all_qc_check = qc_check_task_waphl.all_qc_check
    String? all_qc_alert = qc_check_task_waphl.all_qc_alert
    # Ecoli Typing
    File? serotypefinder_report = merlin_magic.serotypefinder_report
    String? serotypefinder_docker = merlin_magic.serotypefinder_docker
    String? serotypefinder_serotype = merlin_magic.serotypefinder_serotype
    File? ectyper_results = merlin_magic.ectyper_results
    String? ectyper_version = merlin_magic.ectyper_version
    String? ectyper_predicted_serotype = merlin_magic.ectyper_predicted_serotype
    String? shigatyper_predicted_serotype = merlin_magic.shigatyper_predicted_serotype
    String? shigatyper_ipaB_presence_absence = merlin_magic.shigatyper_ipaB_presence_absence
    String? shigatyper_notes = merlin_magic.shigatyper_notes
    File? shigatyper_hits_tsv = merlin_magic.shigatyper_hits_tsv
    File? shigatyper_summary_tsv = merlin_magic.shigatyper_summary_tsv
    String? shigatyper_version = merlin_magic.shigatyper_version
    String? shigatyper_docker = merlin_magic.shigatyper_docker
    File? shigeifinder_report = merlin_magic.shigeifinder_report
    String? shigeifinder_docker = merlin_magic.shigeifinder_docker
    String? shigeifinder_version = merlin_magic.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence = merlin_magic.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes = merlin_magic.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster = merlin_magic.shigeifinder_cluster
    String? shigeifinder_serotype = merlin_magic.shigeifinder_serotype
    String? shigeifinder_O_antigen = merlin_magic.shigeifinder_O_antigen
    String? shigeifinder_H_antigen = merlin_magic.shigeifinder_H_antigen
    String? shigeifinder_notes = merlin_magic.shigeifinder_notes
    # ShigeiFinder outputs but for task that uses reads instead of assembly as input
    File? shigeifinder_report_reads = merlin_magic.shigeifinder_report
    String? shigeifinder_docker_reads = merlin_magic.shigeifinder_docker
    String? shigeifinder_version_reads = merlin_magic.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence_reads = merlin_magic.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes_reads = merlin_magic.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster_reads = merlin_magic.shigeifinder_cluster
    String? shigeifinder_serotype_reads = merlin_magic.shigeifinder_serotype
    String? shigeifinder_O_antigen_reads = merlin_magic.shigeifinder_O_antigen
    String? shigeifinder_H_antigen_reads = merlin_magic.shigeifinder_H_antigen
    String? shigeifinder_notes_reads = merlin_magic.shigeifinder_notes
    # E coli only typing
    File? virulencefinder_report_tsv = merlin_magic.virulencefinder_report_tsv
    String? virulencefinder_docker = merlin_magic.virulencefinder_docker
    String? virulencefinder_hits = merlin_magic.virulencefinder_hits
    # Shigella sonnei Typing
    File? sonneityping_mykrobe_report_csv = merlin_magic.sonneityping_mykrobe_report_csv
    File? sonneityping_mykrobe_report_json = merlin_magic.sonneityping_mykrobe_report_json
    File? sonneityping_final_report_tsv = merlin_magic.sonneityping_final_report_tsv
    String? sonneityping_mykrobe_version = merlin_magic.sonneityping_mykrobe_version
    String? sonneityping_mykrobe_docker = merlin_magic.sonneityping_mykrobe_docker
    String? sonneityping_species = merlin_magic.sonneityping_species
    String? sonneityping_final_genotype = merlin_magic.sonneityping_final_genotype
    String? sonneityping_genotype_confidence = merlin_magic.sonneityping_genotype_confidence
    String? sonneityping_genotype_name = merlin_magic.sonneityping_genotype_name
    # Listeria Typing
    File? lissero_results = merlin_magic.lissero_results
    String? lissero_version = merlin_magic.lissero_version
    String? lissero_serotype = merlin_magic.lissero_serotype
    # Pseudomonas Aeruginosa Typing
    String? pasty_serogroup = merlin_magic.pasty_serogroup
    Float? pasty_serogroup_coverage = merlin_magic.pasty_serogroup_coverage
    Int? pasty_serogroup_fragments = merlin_magic.pasty_serogroup_fragments
    File? pasty_summary_tsv = merlin_magic.pasty_summary_tsv
    File? pasty_blast_hits = merlin_magic.pasty_blast_hits
    File? pasty_all_serogroups = merlin_magic.pasty_all_serogroups
    String? pasty_version = merlin_magic.pasty_version
    String? pasty_docker = merlin_magic.pasty_docker
    String? pasty_comment = merlin_magic.pasty_comment
    # Salmonella Typing
    File? sistr_results = merlin_magic.sistr_results
    File? sistr_allele_json = merlin_magic.sistr_allele_json
    File? sistr_allele_fasta = merlin_magic.sistr_allele_fasta
    File? sistr_cgmlst = merlin_magic.sistr_cgmlst
    String? sistr_version = merlin_magic.sistr_version
    String? sistr_predicted_serotype = merlin_magic.sistr_predicted_serotype
    String? seqsero2_report = merlin_magic.seqsero2_report
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
    # Klebsiella Typing
    File? kleborate_output_file = merlin_magic.kleborate_output_file
    String? kleborate_version = merlin_magic.kleborate_version
    String? kleborate_docker = merlin_magic.kleborate_docker
    String? kleborate_key_resistance_genes = merlin_magic.kleborate_key_resistance_genes
    String? kleborate_genomic_resistance_mutations = merlin_magic.kleborate_genomic_resistance_mutations
    String? kleborate_mlst_sequence_type = merlin_magic.kleborate_mlst_sequence_type
    String? kleborate_klocus = merlin_magic.kleborate_klocus
    String? kleborate_ktype = merlin_magic.kleborate_ktype
    String? kleborate_olocus = merlin_magic.kleborate_olocus
    String? kleborate_otype = merlin_magic.kleborate_otype
    String? kleborate_klocus_confidence = merlin_magic.kleborate_klocus_confidence
    String? kleborate_olocus_confidence = merlin_magic.kleborate_olocus_confidence
    String? kleborate_virulence_score = merlin_magic.kleborate_virulence_score
    String? kleborate_resistance_score = merlin_magic.kleborate_resistance_score
    # Neisseria gonorrhoeae Typing
    File? ngmaster_tsv = merlin_magic.ngmaster_tsv
    String? ngmaster_version = merlin_magic.ngmaster_version
    String? ngmaster_ngmast_sequence_type = merlin_magic.ngmaster_ngmast_sequence_type
    String? ngmaster_ngmast_porB_allele = merlin_magic.ngmaster_ngmast_porB_allele
    String? ngmaster_ngmast_tbpB_allele = merlin_magic.ngmaster_ngmast_tbpB_allele
    String? ngmaster_ngstar_sequence_type = merlin_magic.ngmaster_ngstar_sequence_type
    String? ngmaster_ngstar_penA_allele = merlin_magic.ngmaster_ngstar_penA_allele
    String? ngmaster_ngstar_mtrR_allele = merlin_magic.ngmaster_ngstar_mtrR_allele
    String? ngmaster_ngstar_porB_allele = merlin_magic.ngmaster_ngstar_porB_allele
    String? ngmaster_ngstar_ponA_allele = merlin_magic.ngmaster_ngstar_ponA_allele
    String? ngmaster_ngstar_gyrA_allele = merlin_magic.ngmaster_ngstar_gyrA_allele
    String? ngmaster_ngstar_parC_allele = merlin_magic.ngmaster_ngstar_parC_allele
    String? ngmaster_ngstar_23S_allele = merlin_magic.ngmaster_ngstar_23S_allele
    # Neisseria meningitidis Typing
    File? meningotype_tsv = merlin_magic.meningotype_tsv
    String? meningotype_version = merlin_magic.meningotype_version
    String? meningotype_serogroup = merlin_magic.meningotype_serogroup
    String? meningotype_PorA = merlin_magic.meningotype_PorA
    String? meningotype_FetA = merlin_magic.meningotype_FetA
    String? meningotype_PorB = merlin_magic.meningotype_PorB
    String? meningotype_fHbp = merlin_magic.meningotype_fHbp
    String? meningotype_NHBA = merlin_magic.meningotype_NHBA
    String? meningotype_NadA = merlin_magic.meningotype_NadA
    String? meningotype_BAST = merlin_magic.meningotype_BAST
    # Acinetobacter Typing
    File? kaptive_output_file_k = merlin_magic.kaptive_output_file_k
    File? kaptive_output_file_oc = merlin_magic.kaptive_output_file_oc
    String? kaptive_version = merlin_magic.kaptive_version
    String? kaptive_k_locus = merlin_magic.kaptive_k_match
    String? kaptive_k_type = merlin_magic.kaptive_k_type
    String? kaptive_kl_confidence = merlin_magic.kaptive_k_confidence
    String? kaptive_oc_locus = merlin_magic.kaptive_oc_match
    String? kaptive_ocl_confidence = merlin_magic.kaptive_oc_confidence
    File? abricate_abaum_plasmid_tsv = merlin_magic.abricate_results
    String? abricate_abaum_plasmid_type_genes = merlin_magic.abricate_genes
    String? abricate_database = merlin_magic.abricate_database
    String? abricate_version = merlin_magic.abricate_version
    String? abricate_docker = merlin_magic.abricate_docker
    # Mycobacterium Typing
    File? tbprofiler_output_file = merlin_magic.tbprofiler_output_file
    File? tbprofiler_output_bam = merlin_magic.tbprofiler_output_bam
    File? tbprofiler_output_bai = merlin_magic.tbprofiler_output_bai
    File? tbprofiler_output_vcf = merlin_magic.tbprofiler_output_vcf
    String? tbprofiler_version = merlin_magic.tbprofiler_version
    String? tbprofiler_main_lineage = merlin_magic.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = merlin_magic.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = merlin_magic.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = merlin_magic.tbprofiler_resistance_genes
    File? tbp_parser_lims_report_csv = merlin_magic.tbp_parser_lims_report_csv
    File? tbp_parser_looker_report_csv = merlin_magic.tbp_parser_looker_report_csv
    File? tbp_parser_laboratorian_report_csv = merlin_magic.tbp_parser_laboratorian_report_csv
    File? tbp_parser_coverage_report = merlin_magic.tbp_parser_coverage_report
    Float? tbp_parser_genome_percent_coverage = merlin_magic.tbp_parser_genome_percent_coverage
    Float? tbp_parser_average_genome_depth = merlin_magic.tbp_parser_average_genome_depth
    File? clockwork_decontaminated_read1 = merlin_magic.clockwork_cleaned_read1
    File? clockwork_decontaminated_read2 = merlin_magic.clockwork_cleaned_read2
     # Legionella pneumophila typing
    File? legsta_results = merlin_magic.legsta_results
    String? legsta_predicted_sbt = merlin_magic.legsta_predicted_sbt
    String? legsta_version = merlin_magic.legsta_version
    # Staphylococcus aureus
    File? spatyper_tsv = merlin_magic.spatyper_tsv
    String? spatyper_docker = merlin_magic.spatyper_docker
    String? spatyper_repeats = merlin_magic.spatyper_repeats
    String? spatyper_type = merlin_magic.spatyper_type
    String? spatyper_version = merlin_magic.spatyper_version
    File? staphopiasccmec_results_tsv = merlin_magic.staphopiasccmec_results_tsv
    File? staphopiasccmec_hamming_distance_tsv = merlin_magic.staphopiasccmec_hamming_distance_tsv
    String? staphopiasccmec_types_and_mecA_presence = merlin_magic.staphopiasccmec_types_and_mecA_presence
    String? staphopiasccmec_version = merlin_magic.staphopiasccmec_version
    String? staphopiasccmec_docker = merlin_magic.staphopiasccmec_docker
    File? agrvate_summary = merlin_magic.agrvate_summary
    File? agrvate_results = merlin_magic.agrvate_results
    String? agrvate_agr_group = merlin_magic.agrvate_agr_group
    String? agrvate_agr_match_score = merlin_magic.agrvate_agr_match_score
    String? agrvate_agr_canonical = merlin_magic.agrvate_agr_canonical
    String? agrvate_agr_multiple = merlin_magic.agrvate_agr_multiple
    String? agrvate_agr_num_frameshifts = merlin_magic.agrvate_agr_num_frameshifts
    String? agrvate_version = merlin_magic.agrvate_version
    String? agrvate_docker = merlin_magic.agrvate_docker
    # Streptococcus pneumoniae Typing
    String? pbptyper_predicted_1A_2B_2X = merlin_magic.pbptyper_predicted_1A_2B_2X
    File? pbptyper_pbptype_predicted_tsv = merlin_magic.pbptyper_pbptype_predicted_tsv
    String? pbptyper_version = merlin_magic.pbptyper_version
    String? pbptyper_docker = merlin_magic.pbptyper_docker
    String? poppunk_gps_cluster = merlin_magic.poppunk_gps_cluster
    File? poppunk_gps_external_cluster_csv = merlin_magic.poppunk_gps_external_cluster_csv
    String? poppunk_GPS_db_version = merlin_magic.poppunk_GPS_db_version
    String? poppunk_version = merlin_magic.poppunk_version
    String? poppunk_docker = merlin_magic.poppunk_docker
    String? seroba_version = merlin_magic.seroba_version
    String? seroba_docker = merlin_magic.seroba_docker
    String? seroba_serotype = merlin_magic.seroba_serotype
    String? seroba_ariba_serotype = merlin_magic.seroba_ariba_serotype
    String? seroba_ariba_identity = merlin_magic.seroba_ariba_identity
    File? seroba_details = merlin_magic.seroba_details
    # Streptococcus pyogenes Typing
    String? emmtypingtool_emm_type = merlin_magic.emmtypingtool_emm_type
    File? emmtypingtool_results_xml = merlin_magic.emmtypingtool_results_xml
    String? emmtypingtool_version = merlin_magic.emmtypingtool_version
    String? emmtypingtool_docker = merlin_magic.emmtypingtool_docker
    # Haemophilus influenzae Typing
    String? hicap_serotype = merlin_magic.hicap_serotype
    String? hicap_genes = merlin_magic.hicap_genes
    File? hicap_results_tsv = merlin_magic.hicap_results_tsv
    String? hicap_version = merlin_magic.hicap_version
    String? hicap_docker = merlin_magic.hicap_docker
    # Vibrio Typing
    File? srst2_vibrio_detailed_tsv = merlin_magic.srst2_vibrio_detailed_tsv
    String? srst2_vibrio_version = merlin_magic.srst2_vibrio_version
    String? srst2_vibrio_ctxA = merlin_magic.srst2_vibrio_ctxA
    String? srst2_vibrio_ompW = merlin_magic.srst2_vibrio_ompW
    String? srst2_vibrio_toxR = merlin_magic.srst2_vibrio_toxR
    String? srst2_vibrio_biotype = merlin_magic.srst2_vibrio_biotype
    String? srst2_vibrio_serogroup = merlin_magic.srst2_vibrio_serogroup
  }
}