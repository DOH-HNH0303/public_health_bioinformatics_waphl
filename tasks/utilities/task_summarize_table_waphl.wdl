version 1.0

task summarize_string_data {
  input {
    String terra_project
    String terra_workspace
    String terra_table
    Array[String] samplenames
    Int disk_size = 100
  }
  command <<<   
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --entity_type ~{terra_table} --tsv_filename table-data.tsv 
    echo ~{sep=' ' samplenames}>list.txt
    cat list.txt

    python3<<CODE
    import pandas as pd
    import csv

    with open("list.txt") as file:
      tsv_list = list(csv.reader(file, delimiter=" "))[0]
      print(tsv_list)

    df = pd.read_csv("table-data.tsv", sep="\t")   
    df = df[df[df.columns[0]].isin(tsv_list)]

    df.to_csv("~{terra_table}_data.tsv", sep="\t", index=False)

    CODE    
  >>>
  output {
    File summarized_data = "~{terra_table}_data.tsv"
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}

task zip_files {
  input {
    Array[File]? clade_trees
    Array[File]? recomb_gff
    Array[File]? pirate_aln_gff
    Array[File]? pirate_gene_presence_absence
    String cluster_name
    String? cluster_tree
    Int disk_size = 100
  }
  command <<<   
    # when running on terra, comment out all input_table mentions
    ls
    mkdir ~{cluster_name}
    mv ~{sep=' ' clade_trees} ~{cluster_name}
    mv ~{sep=' ' recomb_gff} ~{cluster_name}
    mv ~{sep=' ' pirate_aln_gff} ~{cluster_name}
    mv ~{sep=' ' pirate_gene_presence_absence} ~{cluster_name}
    mv ~{cluster_tree} ~{cluster_name}

    cd ~{cluster_name}
    tar -cvzf ~{cluster_name}-archive.tar.gz *
    mv ~{cluster_name}-archive.tar.gz ../
    
  >>>
  output {
    File zipped_output = "~{cluster_name}-archive.tar.gz"
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}

task join_existing_tables {
  input {
    String terra_project
    String terra_workspace
    String? old_terra_table
    Array[String] filter = ['assembly_fasta']
    Array[String] samplenames
    Array[String] new_terra_tables
    Int? depth_cov = 50
    Float? busco = 95.0
    Int? avg_q = 30
    Int? contam = 5
    Int? perc_n = 5
    Int disk_size = 100
  }
  command <<<   
    # when running on terra, comment out all input_table mentions
    # Pull an existing table to add to *if* old_terra_table has been provided
    ~{'python3 /scripts/export_large_tsv/export_large_tsv.py --project '+terra_project+' --workspace '+terra_workspace+' --entity_type '+old_terra_table+' --tsv_filename data_old.tsv'} 
    echo ~{sep=' ' samplenames}>list.txt
    table_array=(~{sep=' ' new_terra_tables})

    for index in ${!table_array[@]}; do
      table=${table_array[$index]}
      python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --entity_type "${table}" --tsv_filename "${table}"-data.tsv 
    done

    cat list.txt
    
    python3<<CODE
    import pandas as pd
    import csv
    import os.path
    import datetime 

    
    currentDate = datetime.date.today()
    currentDate = currentDate.strftime("%d%m%Y")
    metrics = {"est_coverage":"~{depth_cov}", "busco_results":"~{busco}", "r1_mean_q":"~{avg_q}", "kraken2_clean_human":"~{contam}", "perc_n":"~{perc_n}"}

    def common_elements(list1, list2):
      return [element for element in list1 if element in list2]

    def pass_filter(data, metrics):\
      id_col = df.columns[0]
      data = data['sample_id'] = df[id_col].str[:9]
      data['perc_n'] = data['number_N']/data["number_Total"]
      data['busco'] = data['busco_results'].str[2:].split("%")[0]
      for filt in ~{filter}:
        data = data[data[filt].notna()]
      for metric in list(metrics.keys()):
                if metric in ['perc_n', 'contam']:
                  data = data[~(data[metric] <= metrics[metric])]
                else:
                  data = data[~(data[metric] >= metrics[metric])]
      return data
    

    #metrics = {"est_coverage":"~{depth_cov}", "busco_results":"~{busco}", "r1_mean_q":"~{avg_q}", "kraken2_clean_human":"~{contam}", "perc_n":"~{perc_n}"}

    check_file = os.path.isfile('./data_old.tsv')
    if check_file:
      df = pd.read_csv("data_old.tsv", sep="\t")
      df = pass_filter(data, metrics)

    for filename in os.scandir('./'):
      if filename.is_file():
        if "data_old.tsv" is in filename:
          pass
        else:
          if df not in locals():
            df = pd.read_csv(filename, sep="\t")
            df = pass_filter(df, metrics)
          else:
            temp_df = pd.read_csv(filename, sep="\t")
            temp_df = pass_filter(temp_df, metrics)
            common = common_elements(temp_df['sample_id'], df['sample_id'])
            for i in common:
              for metric in metrics:
                temp_val =temp_df.loc[temp_df['sample_id'] == i, metrics[metric]].iloc[0]
                df_val = df.loc[temp_df['sample_id'] == i, metrics[metric]].iloc[0]
                if metric in ['perc_n', 'contam']:
                  best = [temp_val, df_val].sort()[0]
                else:
                  best = [temp_val, df_val].sort()[-1]
                if temp_val == best:
                  df = df[df.metric != i]
                if df_val == best:
                  temp_df = temp_df[temp_df.metric != i]
            df = df.join(df_temp.set_index('sample_id'), on='sample_id')
            


    with open("list.txt") as file:
      tsv_list = list(csv.reader(file, delimiter=" "))[0]
      print(tsv_list)

    df = pd.read_csv("table-data.tsv", sep="\t")   
    df = df[df[df.columns[0]].isin(tsv_list)]

    df.to_csv("~{old_terra_table}_+"currentDate+"_dump.tsv", sep="\t", index=False)

    CODE    
  >>>
  output {
    File summarized_data = select_first(glob("*_dump.tsv"))
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}