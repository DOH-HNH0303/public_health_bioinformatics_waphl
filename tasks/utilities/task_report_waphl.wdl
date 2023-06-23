version 1.0

task cdip_report {
  input {
    File assembly_tsv
    Array[File] mlst_tsvs
    Array[File]? clade_trees
    Array[File]? phylo_zip
    Array[File]? plot_roary
    File treefile
    String author
    String project_name = ""
    String cluster_name
    String docker = "hnh0303/seq_report_generator:1.0"
    Int threads = 6 
  }
  command <<<
    # date and version control
    date | tee DATE
    if [ -z ~{sep=' ' phylo_zip} ]; then
    for x in ~{sep=' ' phylo_zip}
    do  
        tar xzf "${x}" --one-top-level=$(basename "${x}" | cut -d. -f1)_phylo
    done;
    fi


    mkdir mlst_tsvs 
    
    for x in ~{sep=' ' mlst_tsvs}
    do
        echo "${x}"
        mv "${x}" mlst_tsvs
        echo "second ${x}"
    done;
    

    if [ -z ~{sep=' ' clade_trees} ]; then
    mkdir clade_trees
    for x in ~{sep=' ' clade_trees}
    do
        mv "${x}" clade_trees
    done;
    fi

    mv ~{assembly_tsv} assembly.tsv
    mv ~{treefile} file.tree

    echo "test2" 
    echo "~{sep=' ' plot_roary}"
    mkdir roary
    for x in ~{sep=' ' plot_roary}
    do
        echo  "${x}"
        mv "${x}" roary
        echo  "second ${x}"
    done;

    # cat /epi_reports/seq_report_generator.py
    python3<<CODE

    from fpdf import FPDF
    from PyPDF2 import PdfFileReader, PdfReader, PdfWriter
    import numpy as np
    import pandas as pd
    from Bio import Phylo
    import matplotlib.pyplot as plt
    from pathlib import Path
    import sys
    import os
    import re
    import datetime
    sys.path.append('/epi_reports')
    from seq_report_generator import add_dendrogram_as_pdf, add_page_header, add_section_header, add_paragraph, add_table, combine_similar_columns, create_dt_col, create_dummy_data, join_pdfs, remove_nan_from_list, unique, new_pdf, add_image


    for subdir, dirs, files in os.walk('.'):
      for file in files:
        
        if subdir == "./mlst_tsvs":
          if "mlst_df" in locals():

             with open(os.path.join(subdir, file)) as lines:
              count = 0
              #print("not first", lines)
              for line in lines:
                line.replace("\n", "")
                if count == 0:
                  cols = line.split("\t")
                  #print("not first, first line", line)
                else:  
                  line_count = 0
                  line = line.split()
                  row = line[:3]
                  allele_row = str(", ".join(line[3:]))
                  row.append(allele_row)
                  hold_df = pd.DataFrame([row], columns=cols)  
                  #print("hold_df", hold_df) 
                  #print("mlst_df before", mlst_df)
                  mlst_df = pd.concat([mlst_df, hold_df], axis=0)#.reset_index(drop=True, inplace=True)
                  #print("mlst_df after", mlst_df)
                count += 1
          else:
            with open(os.path.join(subdir, file)) as lines:
              count = 0
              #print("first", lines)
              for line in lines:
                line.replace("\n", "")
                if count == 0:
                  #print("first line", line)
                  cols = line.split("\t")
                else:  
                  line_count = 0
                  line = line.split()
                  row = line[:3]
                  allele_row = str(", ".join(line[3:]))
                  row.append(allele_row)
                  #print(row)
                  mlst_df = pd.DataFrame([row], columns=cols)    
                count += 1
    #print(mlst_df)
    mlst_df.to_csv('file1.tsv', sep="\t")
    df_mlst = mlst_df.copy()

    # Reformat mlst df to be a by-allele table
    id_list = []
    for i in df_mlst.columns[0]:
      id_list.append(i.split("_"))

    new_cols = ["Sequence ID", "Sequence Type"]
    mlst_count = 0
    alleles = df_mlst.columns[3]

    new_mlst = []
    #print(range(len(df_mlst[alleles].tolist())))
    for i in range(len(df_mlst[alleles].tolist())):
      if i == 0:
        test=re.sub("[\(\[].*?[\)\]]", "", df_mlst[alleles].tolist()[i])
        for allele in test.split(", "):
          new_cols.append(allele)
      hold = []
      for char in df_mlst[alleles].tolist()[i]:
          if char.isalpha():
            pass
          else:
            hold.append(str(char).split("(")[0].split(")")[0])

      hold = "".join(hold).split(", ")
      hold.insert(0,df_mlst[df_mlst.columns[2]].tolist()[i])
      hold.insert(0, df_mlst[df_mlst.columns[0]].tolist()[i].split("_")[0])
      #print(hold, type(hold))
      new_mlst.append(hold)
      
    df_mlst = pd.DataFrame(new_mlst, columns=new_cols) 


    df = pd.read_csv("assembly.tsv", sep="\t")
    df = df[df['assembly_fasta'].notna()]
    df.rename(columns={df.columns[0]: 'Seq ID', "ts_mlst_predicted_st": 'ST Type', "fastani_genus_species": "Species ID"},inplace=True)
    df = create_dt_col(df)

    amr_col = combine_similar_columns(df, ['abricate_amr_genes', 'amrfinderplus_amr_genes'])
    vir_col = combine_similar_columns(df,['abricate_virulence_genes', 'amrfinderplus_virulence_genes'] )

    df['AMR Genes'] = amr_col
    df['Virulence Genes'] = vir_col

    df = create_dummy_data(df, "Submitter", "Lorem ipsum")
    df = create_dummy_data(df, "Collection Date", "01.01.1000")
    df = create_dt_col(df)

    df_genes = df[['Seq ID', 'Collection Date', 'Submitter', 'Species ID', 'ST Type', 'AMR Genes', 'Virulence Genes', "DT"]]
    data = df_genes.copy()
    pdf_report = new_pdf()
    add_page_header(pdf_report, "Corynebacterium Sequencing Report")
    current_day = datetime.date.today()
    formatted_date = datetime.date.strftime(current_day, "%m-%d-%Y")
    print("before")
    header_table = pd.DataFrame([[formatted_date, "~{project_name}", "~{author}"]], columns=["Date", "Project Name", "Prepared By"]) 
    print("after")
    print(header_table)
    add_table(pdf_report, header_table, col_len_dict={})
    pdf_report.ln()
    #add_paragraph(pdf_report, text = "Are you going to want general outbreak info here? You have the option to pass me a text file to add here as a description")
    
    add_section_header(pdf_report, "General Sequencing Results")
    add_table(pdf_report, df_genes)
    
    pdf_report.ln()
    add_section_header(pdf_report, "MLST Results")
    add_table(pdf_report, df_mlst)  

    pdf_report.output("~{cluster_name}"+'_temp_report.pdf', 'F')
    add_dendrogram_as_pdf(pdf_report, tree_file="file.tree", output_filename="~{cluster_name}_tree.pdf")
    plots = new_pdf()
    plots.add_page()
    print("")
    for subdir, dirs, files in os.walk('.'):
      for file in files:

        #print("file", file, subdir, dir)
        if subdir == "./roary":
          print("in roary")
          add_image(plots, os.path.join(subdir, file))
    
    add_paragraph(plots, \
    "1.* Recombination events were not predicted in present gene\n\
    2.* Recombination is unique and is the only one on the gene in dataset\n\
    3.* Multiple recombination events on gene\n\
    4.* Recombination is NOT terminal and is the only recombination event on gene in dataset")

    plots.ln()
    plots.ln()
    plots.ln()

    add_section_header(plots, "Disclaimer")
    add_paragraph(plots,"The information included in this report may be used to support \
    infection prevention measures.  This report should not be used as a substitute for \
    diagnostic procedures, used to guide clinical decisions, for clinical use, nor \
    should it be included in patient records.  The performance characteristics of this \
    test were determined by the Washington State Public Health Laboratories.  \
    It has not been cleared or approved by the U.S. Food and Drug Administration.", italic=False) 

    plots.output("~{cluster_name}"+'_plots.pdf', 'F')
    join_pdfs(["~{cluster_name}_temp_report.pdf", "~{cluster_name}_tree.pdf", "~{cluster_name}_plots.pdf"], "~{cluster_name}_report.pdf")

    CODE
  >>>
  output {
    String date = read_string("DATE")
    File report = "~{cluster_name}_report.pdf"
    File compiled_mlst = "file1.tsv"  
    String split_clade_docker_image = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task plot_roary_waphl {
  input {
    String cluster_name
    File treefile
    File? recomb_gff
    File? pirate_aln_gff
    File? pirate_for_scoary_csv
    String docker = "hnh0303/plot_roary_waphl:1.0"
    Int threads = 6
    String? stripped = basename(treefile, ".tree") 
  }
  command <<<
    # date and version control
    # This task takes either a zipped input file or individual input files
    date | tee DATE
    python ../roary_plots_waphl.py \
    ~{'--recombinants ' + recomb_gff} \
    ~{treefile} \
    ~{pirate_for_scoary_csv} \
    ~{pirate_aln_gff}

    if  [[ "~{stripped}" != "None" ]]; then   
    mv pangenome_matrix.png ~{cluster_name}_matrix.png
    fi
    ls
  >>>
  output {
    String date = read_string("DATE")
    File? plot_roary_png = "~{cluster_name}_matrix.png"
    String plot_roary_docker_image = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
    continueOnReturnCode: "True"
  }
}

task save_output {
  input {
    String cluster_name
    File treefile
    File? recomb_gff
    File pirate_aln_gff
    File pirate_for_scoary_csv
    String cluster_name
    String docker = "hnh0303/plot_roary_waphl:1.0"
    Int threads = 6
    Int? snp_clade
  }
  command <<<
    # date and version control
    date | tee DATE

    python ../roary_plots_waphl.py \
    ~{'--recombinants' + recomb_gff} \
    ~{treefile} \
    ~{pirate_for_scoary_csv} \
    ~{pirate_aln_gff}

   
    mv pangenome_matrix.png ~{cluster_name}~{"_" + snp_clade}_matrix.png

  >>>
  output {
    String date = read_string("DATE")
    File plot_roary_png = select_first(["~{cluster_name}_matrix.png", "~{cluster_name}_~{snp_clade}_matrix.png"])
    String plot_roary_docker_image = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}