## Functions for clustering and enrichment analyses Neat_RNA-Seq_1.0 script
## by Vered Chalifa-Caspi

##### Binary patterns #####

get_replicate_struct_from_coldata = function (col_data, GROUP) {
  #get a vector of the number of replicates in each treatment group, in the order that they appear in col_data
  #GROUP specifies the column name of the treatment group in col_data
  col_data[[GROUP]] %>% fct_inorder %>% table %>% as.numeric
}

zero_one_sequences = function(lenght){
  return(zero_one_sequences_helper(lenght, "", c()))
}

zero_one_sequences_helper = function(lenght, sequence, vec){
  if (nchar(sequence) == lenght){
    vec <- c(vec, sequence)
    return(vec)
  }
  
  vec <- zero_one_sequences_helper(lenght, paste(sequence, "0", sep = "", collapse = ""), vec)
  vec <- zero_one_sequences_helper(lenght, paste(sequence, "1", sep = "", collapse = ""), vec)
  return(vec)
}

get_all_possible_binary_patterns = function (pattern_length) {
  patterns <- zero_one_sequences(pattern_length)
  cat ("No. of binary patterns: ", length(patterns), " (", patterns[1:4], "... )") 
  return (patterns)
}

calc_pattern_1s_pass_min_counts = function (expr_data, patterns, replicates_structure, count_cutoff) {
  
  #create a matrix based on the norm counts.
  #every gene is marked "TRUE" for a pattern if all the samples that have to be 1 according to the pattern, have counts above a cutoff
  
  #count_cutoff = 500           #for testing - to run the function directly
  #expr_data = norm_counts_DE   #for testing - to run the function directly
  
  #create matrix with all results
  #in the matrix, an _ is appended to pattern name, otherwise it does not show well in Excel when exported
  
  pattern_length = nchar(patterns[1])
  
  #initialize results table by just taking row names
  passes = expr_data %>% as.data.frame %>% select(starts_with("999"))   
  
  calc_pass = function (x) {      #receives an expression pattern per gene
    bool_vec = x > count_cutoff   #vector showing for each element of the expression pattern, is it above count cutoff
    indexes = which (pat1 %in% 1) #indexes of samples which has to be 1 according to the pattern
    all (bool_vec[indexes])       #if all the samples which has to be 1 have true in the vector that checks expression above cutoff
  }
  
  for (pattern in patterns) {
    # neglect patterns  "0000" and "1111"
    if (pattern == paste(rep(0, pattern_length), collapse ="") || pattern == paste(rep(1, pattern_length), collapse ="")) next
    # convert pattern to vector, e.g. 1, 0, 1, 1
    pat = str_split(pattern, "") %>% unlist %>% strtoi  
    # match to replicate structure, e.g. 1 1 1 0 0 0 1 1 1 1 1 1
    pat1 = rep(pat, replicates_structure)                               
    # apply
    a = apply(expr_data, 1, calc_pass) %>% as.data.frame   #1 means: work on rows
    passes_col_name = paste0("_", pattern) 
    names(a) = passes_col_name
    passes = cbind(passes, a)
  }
  return (passes)  
}

calc_pattern_1s_and_0s_pass_counts_cutoffs = function(expr_data, patterns, replicates_structure, count_cutoff, count_cutoff_low) {
  
  # Create a matrix based on the norm counts.
  # Every gene is marked "TRUE" for a pattern if:
  # - All samples that have to be "1" according to the pattern have counts > count_cutoff
  # - All samples that have to be "0" according to the pattern have counts < count_cutoff_low
  
  pattern_length = nchar(patterns[1])
  
  # Initialize results table by just taking row names
  passes = expr_data %>% as.data.frame() %>% select(starts_with("999")) 
  
  calc_pass = function(x) {  #receives an expression pattern per gene      
    # Boolean vector for whether each element of the expression pattern meets the count_cutoff or count_cutoff_low condition
    bool_vec_high = x > count_cutoff
    bool_vec_low = x < count_cutoff_low
    
    # Indices for "1" and "0" positions in the pattern
    indexes_1 = which(pat1 == 1)  # indexes of samples which has to be 1 according to the pattern
    indexes_0 = which(pat1 == 0)  # indexes of samples which has to be 0 according to the pattern
    
    # Check conditions for both "1" and "0"
    #if all the samples which has to be 1 have true in the vector that checks expression above cutoff, and same for 0
    all(bool_vec_high[indexes_1]) && all(bool_vec_low[indexes_0])
  }
  
  for (pattern in patterns) {
    # Skip patterns "0000" and "1111"
    if (pattern == paste(rep(0, pattern_length), collapse ="") || pattern == paste(rep(1, pattern_length), collapse ="")) next
    
    # Convert pattern to vector, e.g., "101" -> c(1, 0, 1)
    pat = str_split(pattern, "") %>% unlist() %>% strtoi()
    
    # Match to replicate structure, e.g., "101" -> "111000111"
    pat1 = rep(pat, replicates_structure)
    
    # Apply function to check conditions for each gene
    a = apply(expr_data, 1, calc_pass) %>% as.data.frame()  # 1 means apply function on rows
    passes_col_name = paste0("_", pattern)
    names(a) = passes_col_name
    passes = cbind(passes, a)
  }
  
  return(passes)
}

calc_correlations_to_expression_patterns = function (expr_data, patterns, replicates_structure) {
  
  #create matrix with all correlation results
  #in the matrix, an _ is appended to pattern name, otherwise it does not show well in Excel when exported
  
  pattern_length = nchar(patterns[1])
  
  #initialize results table by just taking row names
  corrs = expr_data %>% as.data.frame %>% select(starts_with("999"))   
  
  for (pattern in patterns) {
    # neglect patterns  "0000" and "1111"
    if (pattern == paste(rep(0, pattern_length), collapse ="") || pattern == paste(rep(1, pattern_length), collapse ="")) next
    # convert pattern to vector, e.g. 1, 0, 1, 1
    pat = str_split(pattern, "") %>% unlist %>% strtoi  
    # match to replicate structure, e.g. 1 1 1 0 0 0 1 1 1 1 1 1
    pat1 = rep(pat, replicates_structure)                               
    # apply
    a = apply(expr_data, 1, cor, pat1) %>% as.data.frame
    corr_col_name = paste0("_", pattern) 
    names(a) = corr_col_name
    corrs = cbind(corrs, a)
  }
  
  return (corrs)
  
}

calc_best_pattern_per_gene = function (corrs, CORR_CUTOFF=0.8, result_file_name, stats_file_name, passes, passes1) {
  
  #create matrix indicating only correlations above cutoff, no. of patterns above cutoff, and the best pattern
  #best pattern is the pattern with the max. correlation, if the correlation is above the cutoff
  
  #best pattern pass counts cutoff shows the best pattern if counts in all samples considered 1 in the pattern have more than a certain no. of counts,
  #based on the "passes" object. The "passes" object is supposed to contain the same genes in the same order as in corrs
  
  get_pattern = function (a) {
    if(max(a)>CORR_CUTOFF){
      col_no = which.max(a)
      col_name = names(a)[col_no]
      col_name1 = str_remove(col_name, "Corr_with")
      return (col_name1)
    } else {
      return (NA)
    }
  }
  
  check_pass = function (a) {
    best_p = as.character (a["best_pattern"])
    if (!is.na (best_p) && isTRUE(as.logical(a[best_p]))) {
      return (best_p)
    } else {
      return (NA)
    }
  }
  
  corrs1 = as.matrix(corrs)
  corrs2 = as.data.frame (ifelse(corrs1 > CORR_CUTOFF, corrs1, NA))
  corrs2$nr_patterns  = apply(corrs, 1, function (x) sum(x>CORR_CUTOFF))
  corrs2$best_pattern = apply(corrs, 1, get_pattern)
  corrs2$best_pattern_pass_counts_cutoff = apply(cbind(passes,   best_pattern=corrs2$best_pattern), 1, check_pass)
  corrs2$best_pattern_pass_counts_cutoff1 = apply(cbind(passes1, best_pattern=corrs2$best_pattern), 1, check_pass)
  
  #print corrs2
  
  export_table(corrs2, result_file_name)
  
  #show and print stats
  
  corr2_stats = table(corrs2$best_pattern)
  cat ("\n")
  cat ("No. genes per pattern")
  print (corr2_stats)
  
  corr2_stats_with_counts_cutoff = table(corrs2$best_pattern_pass_counts_cutoff)
  cat ("No. genes passing counts cutoff")
  print (corr2_stats_with_counts_cutoff)

  corr2_stats_with_counts_cutoff1 = table(corrs2$best_pattern_pass_counts_cutoff1)
  cat ("No. genes passing counts cutoff1")
  print (corr2_stats_with_counts_cutoff1)
    
  #bad code since it uses cbind
  #corr2_stats_merged = as.data.frame(cbind(corr2_stats, corr2_stats_with_counts_cutoff, corr2_stats_with_counts_cutoff1))
  #corr2_stats_merged$Pattern = rownames(corr2_stats_merged)
  #corr2_stats_merged = corr2_stats_merged[,c(3,1,2)]
  #colnames(corr2_stats_merged) = c("Pattern", "No._genes", "No._genes_passing_counts_cutoff")
  
  #fixed code, using left_join (31.12.2024)

  # Convert named vectors to data frames
  df_corr2_stats                     <- data.frame(Pattern = names(corr2_stats),                     No._genes = as.numeric(corr2_stats), stringsAsFactors = FALSE)
  df_corr2_stats_with_counts_cutoff  <- data.frame(Pattern = names(corr2_stats_with_counts_cutoff),  No._genes_passing_counts_cutoff = as.numeric(corr2_stats_with_counts_cutoff), stringsAsFactors = FALSE)
  df_corr2_stats_with_counts_cutoff1 <- data.frame(Pattern = names(corr2_stats_with_counts_cutoff1), No._genes_passing_counts_cutoff1 = as.numeric(corr2_stats_with_counts_cutoff1), stringsAsFactors = FALSE)
  
  # Perform left joins
  corr2_stats_merged <- df_corr2_stats %>%
    left_join(df_corr2_stats_with_counts_cutoff, by = "Pattern") %>%
    left_join(df_corr2_stats_with_counts_cutoff1, by = "Pattern")
  
  # View the merged data frame
  print(corr2_stats_merged)
  
  write.table(x = corr2_stats_merged,
              file = stats_file_name,
              quote = F,
              sep = "\t",
              row.names = F,
              na = "0")
  
  return (list(corrs2=corrs2, binary_patterns_stats=corr2_stats_merged))
}

##### Partition clustering #####

perform_partition_clustering = function (mat2plot, col_data, clustering_dir) {
  
  #compute group means, then do z-scoring to scale expression data by row
  #assuming that col order in mat2plot is the same as row order in col_data

  mat2plot_means = compute_group_means (mat2plot, col_data)
  mat2plot_means_scaled = z_score(mat2plot_means)
  
  #clustering of all DE genes using group means
  #use group averages
  #use the gap algorithm to determine the number of clusters
  #this will run hierarchical clustering and then cut the tree to discrete clusters using the determined no. of clusters
  #it is possible to manually determine the no. of clusters using the k parameter
  #need to feed in scaled data by row, because the stand parameter in eclust scales the data by column which is bad
  
  if (!is.na(K_FIXED)) {
    #clustering with setting no. of clusters
    Clustering <- eclust(as.matrix(mat2plot_means_scaled),
                         stand = F,
                         FUNcluster= "hclust",
                         hc_metric = "pearson",
                         graph = F,  #should be FALSE
                         hc_method = "ward.D2",
                         #k.max =  k.max,
                         k = K_FIXED
    )
  } else {
    gap_maxSE_method = "firstSEmax"
    gap_maxSE_SE.factor = 1 #eclust parameter for determining the number of clusters [Not For Mclust]. higher number is less sensitive. Liron's default is 1
    
    Clustering <- eclust(as.matrix(mat2plot_means_scaled),
                         stand = F,
                         FUNcluster= "hclust",
                         hc_metric = "pearson",
                         graph = F,  #should be FALSE
                         hc_method = "ward.D2",
                         k.max =  K_MAX,
                         nboot = 10,          #no. bootstrap. consider increasing. Liron's default is 10
                         gap_maxSE = list(method= gap_maxSE_method, SE.factor = gap_maxSE_SE.factor)
    )
  }
  
  Clustering$nbclust
  Clustering$size
  
  head(Clustering$cluster)
  head(Clustering$data)
  
  #get "clusters": array where element names are genes and values are cluster no. not re-ordered by clustering.
  clusters=Clustering$cluster
  
  #hierarchical clustering at the group average level of each cluster separately
  #this is then concatenated to create clustering order in all the DE genes (by group averages)
  #the result is used by Liron to draw the heatmap of the averages, but I did not implement it.
  
  New_clusters=c()
  for (num in sort(unique(clusters))){
    res_red=unique(names(clusters))
    Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
    if (length(Genes)>1){
      
      heat_map=pheatmap(mat = as.matrix(mat2plot_means[Genes,]),   #mat2plot_means_Z_score
                        cluster_rows = T,
                        show_rownames=F,
                        clustering_distance_rows = "correlation" ,
                        clustering_method="ward.D2",
                        cluster_cols = F,
                        silent = F,   #to display the graph, comment this line and don't assign to Heat_map, just run pheatmap
                        scale = "row" #should be "row"
      )
      
      New_clusters=c(New_clusters,clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
    }else{
      New_clusters=c(New_clusters,clusters[Genes])
    }
  }
  
  return (New_clusters)
}

analyze_partition_clustering_stats = function  (partition_clusters, new_df, partition_clustering_dir) {
  
  # Summary stats of clusters from partition clustering (after z-scoring)
  # Summary statistics is computed per cluster. Note that FC is converted to linear
  # This will help determine cutoffs for manual clustering
  
  partition_clusters_df = make_data_frame_from_clusters(partition_clusters)
  new_df0 = merge (partition_clusters_df, new_df, by.x="gene", by.y="gene")
  row.names(new_df0) = new_df0$gene
  new_df0 = new_df0[partition_clusters_df$gene,]  #sort lines by order of genes in partition_clusters
  
  agg_mean =
    new_df0 %>%
    select (-gene) %>%
    aggregate.data.frame (by=list(Cluster=new_df0$cluster), FUN=mean) %>%
    mutate (B.vs.A = ifelse(B.vs.A>0,
                            yes = 2^B.vs.A,
                            no = -1/(2^B.vs.A))) %>%
    mutate (C.vs.B = ifelse(C.vs.B>0,
                            yes = 2^C.vs.B,
                            no = -1/(2^C.vs.B))) 
  
  agg_min =
    new_df0 %>%
    select (-gene) %>%
    aggregate.data.frame (by=list(Cluster=new_df0$cluster), FUN=min) %>%
    mutate (B.vs.A = ifelse(B.vs.A>0,
                            yes = 2^B.vs.A,
                            no = -1/(2^B.vs.A))) %>%
    mutate (C.vs.B = ifelse(C.vs.B>0,
                            yes = 2^C.vs.B,
                            no = -1/(2^C.vs.B))) 
  
  agg_max =
    new_df0 %>%
    select (-gene) %>%
    aggregate.data.frame (by=list(Cluster=new_df0$cluster), FUN=max) %>%
    mutate (B.vs.A = ifelse(B.vs.A>0,
                            yes = 2^B.vs.A,
                            no = -1/(2^B.vs.A))) %>%
    mutate (C.vs.B = ifelse(C.vs.B>0,
                            yes = 2^C.vs.B,
                            no = -1/(2^C.vs.B)))
  
  export_table(agg_mean, paste0 (partition_clustering_dir, "/Clusters_mean.txt"))
  export_table(agg_min,  paste0 (partition_clustering_dir, "/Clusters_min.txt"))
  export_table(agg_max,  paste0 (partition_clustering_dir, "/Clusters_max.txt"))
}

##### Manual clustering option 1 #####

perform_manual_clustering_option1 = function (new_df, cluster_def) {
  new_df1 =
    new_df %>%
    mutate (B.vs.A_direction = ifelse(B.vs.A > log_man_clust_FC_cutoff,
                                      yes="up",
                                      no=ifelse(B.vs.A < -log_man_clust_FC_cutoff,
                                                yes="down",
                                                no="nc"))) %>%
    mutate (C.vs.B_direction = ifelse(C.vs.B>log_man_clust_FC_cutoff,
                                      yes="up",
                                      no=ifelse(C.vs.B < -log_man_clust_FC_cutoff,
                                                yes="down",
                                                no="nc"))) %>%
    mutate (Pattern = paste0(B.vs.A_direction, "-", C.vs.B_direction)) %>%
    mutate (cluster = cluster_def[Pattern])
  
  #create a clustering object as a named vector
  man_clusters_opt1 <- setNames(as.numeric(new_df1$cluster), new_df1$gene) 
  
  return (man_clusters_opt1)
}

##### Ranked gene lists #####

create_ranked_genes_by_pval_wo_direction = function(stats_df, contrast) {
  
  #create a data frame: gene, -log10 of pval, sorted by pval (highest -log10(pval) on top)
  pval_col = paste0 ("pvalue.", contrast)
  contrast_df = stats_df %>%
    select (gene, {{pval_col}}) %>%
    setNames(c("gene", "pval")) %>%
    filter(!is.na(pval)) %>%
    mutate(pval = -log10(pval)) %>%
    mutate(pval = ifelse(is.na(pval), yes = 0, no = pval)) %>%
    arrange(desc(pval))
  
  #convert to a named vector and store
  ranked_genes = pull(contrast_df, pval)
  names(ranked_genes) = pull(contrast_df, gene)
  
  return(ranked_genes)
}

create_ranked_genes_by_pval_with_direction = function(stats_df, contrast) {
  
  #create a data frame: gene, log10 of pval with sign according to direction of change, sorted
  pval_col = paste0 ("pvalue.", contrast)
  fc_col   = paste0 ("linearFC.", contrast)
  contrast_df = stats_df %>%
    select (gene, {{fc_col}}, {{pval_col}}) %>%
    setNames(c("gene", "fc", "pval")) %>%
    filter (!is.na(fc)) %>%
    filter (!is.na(pval)) %>%
    mutate(pval = -log10(pval)) %>%
    mutate(pval = ifelse(is.na(pval), yes = 0, no = pval)) %>%
    mutate(pval = ifelse(is.na(fc),
                         yes = 0,
                         no = ifelse(fc>0,
                                     yes = pval,
                                     no = -pval
                         ))) %>%
    select (-fc) %>%
    arrange(desc(pval))
  
  #convert to a named vector
  ranked_genes = pull(contrast_df, pval)
  names(ranked_genes) = pull(contrast_df, gene)
  
  return(ranked_genes)
}

create_ranked_genes_by_fc = function(stats_df, contrast) {
  #I take the pval_col in order to enable filtering by p-val to get the same no. of genes as in the function create_ranked_genes_by_pval_with_direction
  pval_col = paste0 ("pvalue.", contrast)
  fc_col   = paste0 ("linearFC.", contrast)
  contrast_df = stats_df %>%
    select (gene, {{fc_col}}, {{pval_col}}) %>% 
    setNames(c("gene", "fc", "pval")) %>%
    filter (!is.na(fc)) %>%
    filter (!is.na(pval)) %>%
    mutate (fc = ifelse (fc > 0,
                         yes = fc,
                         no = -1/fc)) %>%
    mutate (fc = log2(fc)) %>%
    mutate (fc = signif(fc, digits = 4)) %>%
    arrange(desc(fc))
  
  #convert to a named vector
  ranked_genes = pull(contrast_df, fc)
  names(ranked_genes) = pull(contrast_df, gene)
  
  return(ranked_genes)
}

create_ranked_genes_by_fc_BKP = function(stats_df, contrast) {
  fc_col   = paste0 ("linearFC.", contrast)
  contrast_df = stats_df %>%
    select (gene, {{fc_col}}) %>%
    setNames(c("gene", "fc")) %>%
    filter (!is.na(fc)) %>%
    mutate (fc = ifelse (fc > 0,
                         yes = fc,
                         no = -1/fc)) %>%
    mutate (fc = log2(fc)) %>%
    mutate (fc = signif(fc, digits = 4)) %>%
    arrange(desc(fc))
  
  #convert to a named vector
  ranked_genes = pull(contrast_df, fc)
  names(ranked_genes) = pull(contrast_df, gene)
  
  return(ranked_genes)
}

get_ranked_genes_by_min_pval_any_contrast = function (stats_df) {
  
  #get top n DE genes based on min pval in any comparison, sorted by this pval
  
  pval_cols <- str_detect(string = names(stats_df),
                          pattern = "pvalue")
  
  if (sum(pval_cols, na.rm=TRUE) == 1) {  #if there is only one contrast
    min_pval_vec = stats_df[,pval_cols]
  } else {
    min_pval_vec = apply(X = stats_df[,pval_cols] ,   #else (two or more contrasts)
                         MARGIN=1,
                         FUN = function(x) ifelse(test = all(is.na(x)),
                                                  yes = NA,
                                                  no = min(x,na.rm = T)))
  }
  
  min_pval_df <-
    stats_df %>%
    # Add column with minimum pval
    cbind(min_pval=min_pval_vec) %>%
    select(gene, min_pval) %>%
    filter(!is.na(min_pval)) %>%
    mutate(min_pval = -log10(min_pval)) %>%
    mutate(min_pval = ifelse(is.na(min_pval), yes = 0, no = min_pval)) %>%
    arrange(desc(min_pval))
  
  #convert to a named vector
  ranked_genes = pull(min_pval_df, min_pval)
  names(ranked_genes) = pull(min_pval_df, gene)
  
  return(ranked_genes)
}

get_ranked_genes_by_min_pval_any_contrast_BAD = function (stats_df) {
  
  #this function sends an error when there is only one contrast in the experiment
  
  #get top n DE genes based on min pval in any comparison, sorted by this pval
  
  pval_cols <- str_detect(string = names(stats_df),
                          pattern = "pvalue")
  
  min_pval_df <-
    stats_df %>%
    # Add column with minimum pval
    cbind(min_pval=apply(X = stats_df[,pval_cols] ,
                         MARGIN=1,
                         FUN = function(x) ifelse(test = all(is.na(x)),
                                                  yes = NA,
                                                  no = min(x,na.rm = T)))) %>%
    select(gene, min_pval) %>%
    filter(!is.na(min_pval)) %>%
    mutate(min_pval = -log10(min_pval)) %>%
    mutate(min_pval = ifelse(is.na(min_pval), yes = 0, no = min_pval)) %>%
    arrange(desc(min_pval))
  
  #convert to a named vector
  ranked_genes = pull(min_pval_df, min_pval)
  names(ranked_genes) = pull(min_pval_df, gene)
  
  return(ranked_genes)
}

##### Enrichment analysis #####

Clusters_Enrichment_Test=function(Type, clusters,TERM2NAME,TERM2GENE, outDir, file_name, pvalueCutoff=0.05, pAdjustMethod='fdr', gene2ko=FALSE, maxCategory=1000){
  
  allRes0=list()
  
  #calculate enrichment for each cluster
  for (cluster_name in sort(unique(clusters))){
    
    genes_having_pathway=unique(TERM2GENE[,2])
    genes_in_cluster = get_genes_per_cluster(clusters, cluster_name)
    Genes = intersect(genes_in_cluster, genes_having_pathway)
    
    if (length(Genes)>0) {
      res=enricher(Genes,
                   TERM2GENE     = TERM2GENE, 
                   TERM2NAME     = TERM2NAME,
                   minGSSize     = 0,
                   #maxGSSize    = length(genes_having_pathway),
                   maxGSSize     = 10000,
                   pAdjustMethod = pAdjustMethod,
                   pvalueCutoff  = pvalueCutoff,
                   qvalueCutoff  = 1
      )
      
      #store results in allRes0 list
      #if (length(res)>0){  #6.1.2025 in case that there are enriched pathways but none had p.adj < 0.05, this condition will give TRUE while still the result table will have no results. This will cause error later in the code
      if (nrow(res@result) > 0 && nrow(res@result[res@result$p.adjust < pvalueCutoff,]) > 0)   { #changed on 6.1.2025. was told that simply writing nrow(res) is unreliable
        allRes0[[as.character(cluster_name)]]<-res
      }
    }
  } 
  
  #merge results for all clusters using clusterProfiler::merge_result
  
  if(length(allRes0) == 0) {
    return (list())
  } else {
    allRes=clusterProfiler::merge_result(enrichResultList = allRes0)
    if (Type=="KO"){
      allRes<-generat_urls(allRes, gene2ko)
    }
    
    #create and print enrichment results in excel
    enrichment_table = process_clusterprofiler_results_table (allRes@compareClusterResult)
    enrichment_table_file = file.path(outDir, paste0(file_name,'.csv'))
    write.csv(x         = enrichment_table,
              file      = enrichment_table_file,
              quote     = TRUE,
              row.names = TRUE)
    
    #add AllRes@fun (Vered: I don't know if/why this is necessary)
    if (Type=="GO"){
      allRes@fun<-"enrichGO"
    } else {
      allRes@fun<-"enrichKEGG"
    }
    
    #create and print GO simplify results in excel
    allRes_simplify = NULL
    enrichment_table_simplify = NULL
    if (Type=="GO"){
      allRes_simplify = clusterProfiler::simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min)
      enrichment_table_simplify = process_clusterprofiler_results_table (allRes_simplify@compareClusterResult)
      enrichment_table_simplify_file = file.path (outDir, paste0("Simplify_", file_name, ".csv"))
      write.csv(x         = enrichment_table_simplify,
                file      = enrichment_table_simplify_file,
                quote     = FALSE,
                row.names = TRUE)
    }
    
    #set font size for dot plot
    if (dim(allRes@compareClusterResult)[1]>50){
      font.size=4
    }else{
      font.size=9
    }
    
    #draw dot plot
    if (dim(allRes@compareClusterResult)[1]>0){
      dot_plot_file = file.path(outDir, paste0(file_name, ".pdf"))
      x = enrichplot::dotplot (allRes,
                               showCategory = maxCategory,
                               font.size=font.size)
      ggplot2::ggsave(filename = dot_plot_file,
                      dpi = 600,
                      device = "pdf",
                      width = 20,
                      height = 20)
    }
    
    #return: a list containing enrichment results for all clusters, simplify enrich. results
    return(list(allRes,
                allRes_simplify,
                enrichment_table,
                enrichment_table_simplify))			
  }
}

process_clusterprofiler_results_table = function(clusterprofiler_results_table) {
  
  MAX_NR_GENES_TO_SHOW = 1000
  text_to_show = paste0("Too many to show (>", MAX_NR_GENES_TO_SHOW, ")")
  
  enrichment_table = clusterprofiler_results_table %>%
    separate(col = GeneRatio, into = c("in_cluster_in_term", "in_cluster")) %>%
    separate(col = BgRatio,   into = c("in_term", "in_genome"))  %>%
    mutate_at (c('in_cluster_in_term', 'in_cluster', 'in_term', 'in_genome'), as.numeric) %>%
    mutate (Fold_enrichment = signif(((in_cluster_in_term/in_cluster)/(in_term/in_genome)),digits=2), .before = Count) %>%
    mutate (geneID = ifelse (test = Count <= MAX_NR_GENES_TO_SHOW,
                             yes  = geneID,
                             no   = text_to_show))
  return(enrichment_table)
}

plot_shared_genes<-function(allRes,outDir,file_name){
  
  heatmaps=list()
  if (length(allRes)>1){
    for (cluster in unique(allRes$Cluster)){
      TEMP = na.omit(allRes[allRes$Cluster==cluster,])   
      genes = unique(unlist(stringi::stri_split(paste(TEMP$geneID,sep = ,collapse = '/'),fixed = '/')))
      
      if ((length(TEMP$ID)>1) &&((length(genes)>1))){
        mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(genes))
        rownames(mat) = TEMP$ID
        colnames(mat) = genes
        
        
        for (x in TEMP$ID){
          for (y in genes){
            mat[x,y] = length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),y )))
            
          }
        }
        
        rownames(mat) = TEMP$Description
        
        if ((length(TEMP$ID)>20) || ((length(genes)>200))){
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              cluster_cols = T,
                              silent = T,
                              legend = F)
          
        }else{
          if (length(unique(as.vector(mat)))==2 ){
            color=colorRampPalette(c('white','pink'))(2)
            mybreaks=NA
          }else{
            color=colorRampPalette(c('pink'))(1)
            mybreaks=c(0,1)
          }
          
          heatmap_file = file.path(outDir, paste0('Cluster_',cluster,'_genes2term_',file_name,".pdf"))
          
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              width = 5,
                              height =5,
                              cellheight = 10,
                              fontsize_row =5,
                              fontsize_col = 1,
                              main = paste('Cluster ',cluster),
                              border_color = 'white',
                              cluster_cols = T,
                              breaks = mybreaks,
                              filename = heatmap_file,
                              color = color,
                              legend = F)
          
        }
        
        csv_file = file.path(outDir, paste0('Cluster_',cluster,'_genes2term_',file_name,".csv"))
        
        write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                  quote = T,
                  row.names = T,
                  file = csv_file)
        
        
        mat  = matrix(0, nrow = length(TEMP$ID), ncol = length(TEMP$ID))
        rownames(mat) = TEMP$ID
        colnames(mat) = TEMP$ID
        
        for (x in TEMP$ID){
          for (y in TEMP$ID){
            mat[x,y] =100*(2*length( unlist(intersect( unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/')),unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))))/(length(unlist(stringi::stri_split(TEMP[TEMP$ID==x,'geneID'],fixed = '/'))) +  length(unlist(stringi::stri_split(TEMP[TEMP$ID==y,'geneID'],fixed = '/'))))
            
          }
        }
        
        rownames(mat) = TEMP$Description   
        colnames(mat) = TEMP$Description   
        if (length(unique(as.vector(mat)))>1 ){
          color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu") ))( 100)
          mybreaks=seq(0,1,0.01)
        }else{
          color=colorRampPalette(c('red'))(1)
          mybreaks=c(0,1)
        }
        
        heatmap_file1 = file.path(outDir, paste0('Cluster_', cluster, '_term2term_', file_name,".pdf"))
        
        if (length(TEMP$ID)>20){
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              main = paste('Cluster ',cluster),
                              treeheight_col = 0,
                              treeheight_row = 0,
                              cluster_cols = T,
                              silent = T,
                              legend = T)
          
        }else{
          heat_map = pheatmap(mat = mat,cluster_rows = T,
                              width = 5,
                              height =5,
                              treeheight_col = 0,
                              treeheight_row = 0,
                              breaks = mybreaks,
                              fontsize_row =5,
                              fontsize_col = 5,
                              main = paste('Cluster ',cluster),
                              border_color = 'white',
                              cluster_cols = T,
                              color = color,
                              filename = heatmap_file1,
                              legend = T)
        }
        
        csv_file1 = file.path(outDir,paste0('Cluster_',cluster,'_term2term_',file_name,".csv"))
        
        write.csv(x = mat[heat_map$tree_row$order,heat_map$tree_col$order],
                  quote = T,
                  row.names = T,
                  file = csv_file1)
        
        heatmap_title=htmltools::h4(paste(unlist(stringi::stri_split(str = paste('Cluster',cluster,file_name,' Terms Genes Overlap',collapse = ""),fixed='_')),collapse = " "))
        heatmaps = c(heatmaps ,list(heatmap_title))
        mat = mat[heat_map$tree_row$order,heat_map$tree_col$order]
        heatmaps = c(heatmaps ,list(plot_ly(z = mat,
                                            x=colnames(mat),
                                            y=rownames(mat),
                                            colors = color,
                                            type = "heatmap")
        )
        )
      }
    }
    
  }
  return(heatmaps)
}

generat_urls<-function(allRes,gene2ko){
  
  #not sure this function is necessary (but it is used in other functions above)
  
  temp_table=allRes@compareClusterResult
  temp_table$URL=apply(X = temp_table,MARGIN = 1,FUN = function(x) paste(c("http://www.kegg.jp/kegg-bin/show_pathway?",stringi::stri_replace_all(str = x["ID"],replacement = "",regex = "path:",collapse = ""),"/", paste(sapply( unlist(stringi::stri_split(str = x["geneID"],regex = "/")),FUN = function(x) gene2ko[gene2ko$V1==x,"V2"]),collapse = "+") ),collapse = ""))
  allRes@compareClusterResult=temp_table
  return(allRes)
}

ridgeplot_edited <- function(x, showCategory=30, fill="p.adjust",
                             core_enrichment = TRUE, label_format = 30,
                             orderBy = "NES", decreasing = FALSE, x_axis_title, file,
                             values_for_x_axis=c(),
                             xlimits=c()) {
  #function modified from:
  #https://rdrr.io/github/GuangchuangYu/enrichplot/src/R/ridgeplot.R
  
  #the "fill" parameter can be one of "pvalue", "p.adjust", "qvalue"
  
  #the "values_for_x_axis" parameter is used for cases where the GSEA was done using p-values for ranking,
  #but the user wants the ridgeplot to show FC values on the X-axis.
  #values_for_x_axis has the form of a named vector, with gene IDs as names and log2FC values as elements
  #(the same format as the input to GSEA)
  
  #preparations
  
  if (!is(x, "gseaResult"))
    stop("currently only support gseaResult")
  
  ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
  if (fill == "qvalue") {
    fill <- "qvalues"
  }
  if (!fill %in% colnames(x@result)) {
    stop("'fill' variable not available ...")
  }
  
  ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')
  if (orderBy !=  'NES' && !orderBy %in% colnames(x@result)) {
    message('wrong orderBy parameter; set to default `orderBy = "NES"`')
    orderBy <- "NES"
  }
  
  #create a list of top n pathways, where each pathway points to all/core genes in that pathway.
  
  n <- showCategory
  if (core_enrichment) {
    gs2id <- geneInCategory(x)[seq_len(n)]
  } else {
    gs2id <- x@geneSets[x$ID[seq_len(n)]]
    ridgeplot_data_file = str_replace (file, pattern = "\\.png$", replacement = ".csv")
    file = str_replace (file, pattern = "\\.png$", replacement = "_allG.png")  #allG is allGenes
  }
  
  #create a list of top n pathways, where each pathway points to a named vector with genes as names
  #and ranking values (FC, p-value or signed p-value) as elements.
  #the ranking values are taken either from the GSEA result (stored in x@geneList, which contains the input to GSEA)
  #or a named vector that was passed as argument to the ridgeplot_edited function.
  
  #original code:
  # gs2val <- lapply(gs2id, function(id) {
  #   res <- x@geneList[id]
  #   res <- res[!is.na(res)]
  # })
  #Vered 27.12.2023, to enable plotting FC on the X-axis:
  
  if (length(values_for_x_axis) > 0) {
    genes2values = values_for_x_axis
  } else {
    genes2values = x@geneList
  }
  
  gs2val <- lapply(gs2id, function(id) {
    res <- genes2values[id]
    res <- res[!is.na(res)]
  })
  
  #create vector of pathway titles
  
  #x$ID is the list of pathways in x, which may be larger or smaller than nn.
  #nn length is the same as showCategory
  nn <- names(gs2val)
  i <- match(nn, x$ID)   #similar to nn %>% x$ID 
  nn <- x$Description[i] #vector of pathway titles with length of showCategory
  
  #order pathways (actually, pathway indexes) by orderBy (default: by NES)
  # j <- order(x$NES[i], decreasing=FALSE)
  j <- order(x@result[[orderBy]][i], decreasing = decreasing)
  #create a named vector where each pathway points to the no. of (core) genes in that pathway
  len <- sapply(gs2val, length)
  
  #create a table (data frame) with columns: category (pathway), color, value, where row names are pathway.gene
  gs2val.df <- data.frame(category = rep(nn, times=len),         #here the pathways are still the first 30 pathways according to pvalue       
                          color = rep(x[i, fill], times=len),
                          value = unlist(gs2val))                #Vered, need to check here
  
  #convert p.adjust to -log10(p.adjust)
  gs2val.df = mutate(gs2val.df, color = -log10(color))
  
  #this changes the column name "color" to the argument "fill" (which is either "pvalue", "p.adjust" or "qvalue")
  colnames(gs2val.df)[2] <- fill
  
  #??
  gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])  #this probably orders the levels of the factor
  
  #label_format
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  #print data behind the ridgeplot
  
  gs2val_for_print.df = gs2val.df %>% rownames_to_column(var="name") %>% separate(col=name, into=c("pathway", "gene"), sep="\\.") 
  ridgeplot_data_file = str_replace (file, pattern = "\\.png$", replacement = ".csv")
  print (ridgeplot_data_file)
  write.csv(gs2val_for_print.df, file=ridgeplot_data_file, row.names = F)  #here the pathways are still the first 30 pathways according to pvalue, ordered by pvalue  
  
  #print ridgeplot
  
  fill_title = paste0("-log10(", fill, ")")
  
  print (file)
  png(filename = file, width = 700, height = 700, units = "px")
  
  plot = ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
    ggridges::geom_density_ridges() +
    ## scale_x_reverse() +
    scale_fill_continuous(low="blue", high="red", name = fill_title,
                          guide=guide_colorbar(reverse=F)) +
    scale_y_discrete(labels = label_func) +
    ## scale_fill_gradientn(name = fill, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ## geom_vline(xintercept=0, color='firebrick', linetype='dashed') +
    xlab(x_axis_title) + ylab(NULL) +  theme_dose()
  
  if (length(xlimits)>0) {
    plot = plot + xlim(xlimits[1], xlimits[2])
  }
  
  print(plot)
  dev.off()
  
}

ridgeplot_edited1 <- function(x, showCategory=30, fill="p.adjust",
                              core_enrichment = TRUE, label_format = 30,
                              orderBy = "NES", decreasing = FALSE, x_axis_title, file,
                              values_for_x_axis=c(),
                              xlimits=c()) {
  
  #this function should replace the function ridgeplot_edited. Just, before replacing test it on GSEA results
  
  #function modified from:
  #https://rdrr.io/github/GuangchuangYu/enrichplot/src/R/ridgeplot.R
  
  #the "fill" parameter can be one of "pvalue", "p.adjust", "qvalue"
  
  #the "values_for_x_axis" parameter is used for cases where the GSEA was done using p-values for ranking,
  #but the user wants the ridgeplot to show FC values on the X-axis.
  #values_for_x_axis has the form of a named vector, with gene IDs as names and log2FC values as elements
  #(the same format as the input to GSEA)
  
  #preparations
  
  if (!fill %in% colnames(x@result)) {
    stop("'fill' variable not available ...")
  }
  
  if (!orderBy %in% colnames(x@result)) {
    stop ("'orderBy' variable not available ...")
  }
  
  #create a list of top n pathways, where each pathway points to all/core genes in that pathway.
  
  n <- showCategory
  if (core_enrichment) {
    gs2id <- geneInCategory(x)[seq_len(n)]
    if(length(x$ID) < n) {   #Vered 4.2.2024 necessary for when res is a result of enrichResult (simple enrichment)
      gs2id = gs2id[x$ID]
    }
  } else {
    gs2id <- x@geneSets[x$ID[seq_len(n)]]
    ridgeplot_data_file = str_replace (file, pattern = "\\.png$", replacement = ".csv")
    file = str_replace (file, pattern = "\\.png$", replacement = "_allGenes.png")
  }
  
  #create a list of top n pathways, where each pathway points to a named vector with genes as names
  #and ranking values (FC, p-value or signed p-value) as elements.
  #the ranking values are taken either from the GSEA result (stored in x@geneList, which contains the input to GSEA)
  #or a named vector that was passed as argument to the ridgeplot_edited function.
  
  #original code:
  # gs2val <- lapply(gs2id, function(id) {
  #   res <- x@geneList[id]
  #   res <- res[!is.na(res)]
  # })
  #Vered 27.12.2023, to enable plotting FC on the X-axis:
  
  if (length(values_for_x_axis) > 0) {
    genes2values = values_for_x_axis
  } else {
    genes2values = x@geneList
  }
  
  gs2val <- lapply(gs2id, function(id) {
    res <- genes2values[id]
    res <- res[!is.na(res)]
  })
  
  #create vector of pathway titles
  
  #x$ID is the list of pathways in x, which may be larger or smaller than nn.
  #nn length is the same as showCategory
  nn <- names(gs2val)
  i <- match(nn, x$ID)   #similar to nn %>% x$ID 
  nn <- x$Description[i] #vector of pathway titles with length of showCategory
  
  #order pathways (actually, pathway indexes) by orderBy (default: by NES)
  # j <- order(x$NES[i], decreasing=FALSE)
  j <- order(x@result[[orderBy]][i], decreasing = decreasing)
  #create a named vector where each pathway points to the no. of (core) genes in that pathway
  len <- sapply(gs2val, length)
  
  #create a table (data frame) with columns: category (pathway), color, value, where row names are pathway.gene
  gs2val.df <- data.frame(category = rep(nn, times=len),         #here the pathways are still the first 30 pathways according to pvalue       
                          color = rep(x[i, fill], times=len),
                          value = unlist(gs2val))                #Vered, need to check here
  
  #convert p.adjust to -log10(p.adjust)
  gs2val.df = mutate(gs2val.df, color = -log10(color))
  
  #this changes the column name "color" to the argument "fill" (which is either "pvalue", "p.adjust" or "qvalue")
  colnames(gs2val.df)[2] <- fill
  
  #??
  gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])  #this probably orders the levels of the factor
  
  #label_format
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  #print data behind the ridgeplot
  
  gs2val_for_print.df = gs2val.df %>% rownames_to_column(var="name") %>% separate(col=name, into=c("pathway", "gene"), sep="\\.") 
  ridgeplot_data_file = str_replace (file, pattern = "\\.png$", replacement = ".csv")
  print (ridgeplot_data_file)
  write.csv(gs2val_for_print.df, file=ridgeplot_data_file, row.names = F)  #here the pathways are still the first 30 pathways according to pvalue, ordered by pvalue  
  
  #print ridgeplot
  
  fill_title = paste0("-log10(", fill, ")")
  
  print (file)
  png(filename = file, width = 700, height = 700, units = "px")
  
  plot = ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
    ggridges::geom_density_ridges() +
    ## scale_x_reverse() +
    scale_fill_continuous(low="blue", high="red", name = fill_title,
                          guide=guide_colorbar(reverse=F)) +
    scale_y_discrete(labels = label_func) +
    ## scale_fill_gradientn(name = fill, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    ## geom_vline(xintercept=0, color='firebrick', linetype='dashed') +
    xlab(x_axis_title) + ylab(NULL) +  theme_dose()
  
  if (length(xlimits)>0) {
    plot = plot + xlim(xlimits[1], xlimits[2])
  }
  
  print(plot)
  dev.off()
  
}

draw_ridgeplot_per_cluster = function (clust_method, clust_round, cluster_name, values_for_x_axis, ridgeplot_file, TERM2NAME,TERM2GENE, pvalueCutoff=0.05, pAdjustMethod='fdr') {
  clusters = gene_lists[[clust_method]][[clust_round]]
  genes_having_pathway=unique(TERM2GENE[,2])
  genes_in_cluster = get_genes_per_cluster(clusters, cluster_name)
  Genes = intersect(genes_in_cluster, genes_having_pathway)
  
  res=enricher(Genes,
               TERM2GENE     = TERM2GENE, 
               TERM2NAME     = TERM2NAME,
               minGSSize     = 0,
               #maxGSSize    = length(genes_having_pathway),
               maxGSSize     = 10000,
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff  = pvalueCutoff,
               qvalueCutoff  = 1
  )
  
  if(length(res$ID) > 0) {
    #show DE genes only
    ridgeplot_edited1 (res, x_axis_title = x_axis_title, file = ridgeplot_file, values_for_x_axis = values_for_x_axis, orderBy = "pvalue", decreasing = T)  # showCategory = 20 (default is 30)
    #show all genes
    ridgeplot_edited1 (res, x_axis_title = x_axis_title, file = ridgeplot_file, values_for_x_axis = values_for_x_axis, orderBy = "pvalue", decreasing = T, core_enrichment = F)  
  } else {
    print ("No enrichment")
  }
}

ep_str_wrap <- function(string, width) {
  #utility function from
  #https://rdrr.io/bioc/enrichplot/src/R/utilities.R
  
  #' ep_str_wrap internal string wrapping function
  #' @param string the string to be wrapped
  #' @param width the maximum number of characters before wrapping to a new line
  #' @noRd
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}


default_labeller <- function(n) {
  #utility function from
  #https://rdrr.io/bioc/enrichplot/src/R/utilities.R
  #' default_labeller
  #'
  #' default labeling function that uses the
  #' internal string wrapping function `ep_str_wrap`
  #' @noRd
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}


##### Graphics #####

plot_expression_heatmap_for_precomputed_clusters = function (mat2plot, clusters, col_data, heatmap_file_name,
                                                             row_distance_measure="correlation", plot_title = NA,
                                                             plot_width=15, plot_height=15, row_annot=T, col_annot=T, show_sample_names=T, show_dend=T,
                                                             scale=1.5,
                                                             print_file_per_cluster=T) {
  
  #code adjusted from Liron's DeSeq2 module
  
  # show or hide row dendogram
  if (show_dend) {
    dendogram_height = 50
  } else {
    dendogram_height = 0
  }
  
  New_clusters=c()
  for (num in sort(unique(clusters))){  #for each cluster
    #get a list of the genes in this cluster
    res_red=unique(names(clusters))
    Genes=res_red[res_red  %in% names(clusters[clusters %in% num])]
    
    if (length(Genes)>1){
      
      heatmap_title = paste0("Cluster ", num, " (", length(Genes), " genes)")
      
      #create a heatmap of the genes in the cluster
      heat_map=pheatmap(mat = mat2plot[Genes,],
                        cluster_rows = T,
                        show_rownames = F,
                        clustering_distance_rows = "correlation" ,
                        clustering_method = "ward.D2",
                        main = heatmap_title,
                        cluster_cols = F,
                        silent = F,   
                        scale = "row",
                        show_colnames = show_sample_names,
                        treeheight_row = dendogram_height)
      
      
      if (print_file_per_cluster) {
        file_per_cluster = str_replace(heatmap_file_name, ".png", paste0("_cluster", num, ".png"))
        ggsave(file_per_cluster,heat_map$gtable,dpi = 600, width = plot_width, height = plot_height, units = "cm",scale = 1.5)
      }
      
      
      #add row order as computed by the hierarchical clustering
      #an array where element names are the genes and values are cluster no.
      New_clusters=c(New_clusters,clusters[heat_map$tree_row$labels[heat_map$tree_row$order]])
    }else{
      New_clusters=c(New_clusters,clusters[Genes])
    }
  }
  
  #prepare for final clustering  
  
  clusters=New_clusters
  
  #order the expression matrix by the order from the clustering (each cluster after the other,
  #and within each cluster, the order from the sample-wise hierarchical clustering)
  
  mat2plot_ordered = mat2plot[names(clusters),]
  
  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot_ordered <- mat2plot_ordered[complete.cases(mat2plot_ordered), ]
  
  #create column annotation (sample names) for the heatmap
  #creates a data frame where row names are sample names and there is one column called X_AXIS with the biol groups
  
  #choose one of the options. if GROUP is not equal to GROUP1 you may need to change the order
  col_annotation = col_data[GROUP] #to show only the main effect (e.g. Diet)
  col_annotation = col_data %>% select(all_of(EFFECTS))  #to show all effects
  
  #create row annotation (cluster no.) for the heatmap
  #create a data frame where row names are genes and there is one column called Clusters with the cluster no. clusters are factors
  if (row_annot) {
    annotation_row = data.frame("Clusters" = factor(clusters))
  } else {
    annotation_row = NA
  }
  
  #define row gaps
  row_gaps = cumsum(aggregate.data.frame(x = as.data.frame(clusters),by =list(clusters),FUN = length )[,"clusters"] )
  
  #define col gaps
  if (GROUP == GROUP1) {
    col_gaps = cumsum(aggregate.data.frame(x = col_annotation[GROUP],by =col_annotation[GROUP],FUN = length)[,-1])
  } else {
    #check this! may need to switch between GROUP and GROUP1. Note, this required that the samples are originally ordered as expected
    #(otherwise, see code from Liron how to reorder them according the needs of this hiererchical clustering)
    col_gaps = c(cumsum(rev(aggregate.data.frame(x = col_annotation[GROUP1],by =col_annotation[c(GROUP1,GROUP)],FUN = length))[,1]),
                 rep(cumsum(rev(aggregate.data.frame(x = col_annotation[GROUP],by =col_annotation[c(GROUP)],FUN = length))[,1]),2))    
  }
  
  #show or hide column annotation
  
  if (col_annot) {
    annotation_col = col_annotation
  } else {
    annotation_col = NA
  }
  
  
  
  #do clustering
  
  #create heatmap
  heat_map=pheatmap(mat = mat2plot_ordered ,
                    cluster_rows= F,   #don't cluster rows. they are already ordered.
                    show_rownames=F,
                    cluster_cols=F,
                    scale = "row",
                    silent = F,
                    annotation_row=annotation_row,
                    #annotation_colors = annotation_row_COLORS,
                    #annotation_names_row=T,
                    annotation_col=annotation_col,
                    show_colnames = show_sample_names,
                    gaps_row = row_gaps,
                    gaps_col = col_gaps,
                    treeheight_row = dendogram_height,
                    border_color=NA
  )
  
  ggsave(heatmap_file_name,heat_map$gtable,dpi = 600, width = plot_width, height = plot_height, units = "cm",scale = scale)
  
  
  #for saving the graph in a long format:
  #ggsave(heatmap_file_name,heat_map$gtable,dpi = 600, width = 23, height = 40, units = "cm")
  
  return(New_clusters)
}


plot_cluster_profiles<-function(clusters,mat2plot,col_data, color_group=c('Type','Time'),titles=c("Type","Time","Normalized counts"),split_by=c(),smooth=T,X_AXIS_ORDER=NA, out_file = ""){
  
  #code adjusted from Liron's DeSeq2 module
  
  #arguments hard coded (for testing)
  
  # color_group = c(GROUP1,GROUP)
  # titles = c(GROUP1,GROUP,"Normalized counts")
  # split_by=c()
  # smooth=T
  # X_AXIS_ORDER=unique(col_data[GROUP])
  # out_file = cluster_profiles_file
  
  vs_ano = cbind(col_data[c(GROUP1, GROUP)], t(mat2plot))  #should be col_data[color_group]
  
  plot_list  = list()
  title_list = list()
  count=1
  for (i in sort(unique(clusters))){
    Df=data.frame()
    genes=names(clusters[clusters==i])
    for (j in genes){
      temp=subset.data.frame(x =vs_ano,,c(color_group,j))
      colnames(temp)=c("Type","Time","EXP")
      Df=rbind(Df,temp)
    }
    if (color_group[1]==color_group[2]){
      Df$Type='Trend'
    }
    if (length(X_AXIS_ORDER)>1){
      Df$Time <-factor(Df$Time,levels=intersect(X_AXIS_ORDER,unique(Df$Time)))
    }
    if (length(split_by)>0){
      title_list[count]=list(c(i,length(genes)) ) 
      plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                              facet_wrap(~split_by, ncol=1,scales = "free_y",strip.position = c("right"))+
                              theme(strip.text.y = element_text(size = 7, colour = "black"))+#,face="bold" ,angle = 90))+
                              #theme(strip.background = element_blank())+
                              ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                              #theme(plot.title = element_text(colour = "black", size = 12)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width =unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }else{
      title_list[count]=list(c(i,length(genes)) ) 
      plot_list[count]=list(ggplot(data=Df, aes(x=Time, y=EXP, group=Type)) + 
                              ggtitle( paste(paste("Cluster ", i  ,sep="") ,length(genes) ,sep="\n")) +
                              #theme(plot.title = element_text(colour = "black", size = 12)) + 
                              theme(legend.position  ="none")+
                              xlab(titles[2]) +
                              ylab(titles[3]) +
                              theme(axis.title   = element_text(colour = "black", size = 10) )+
                              theme(legend.title = element_text(colour = "black", size = 10) )+
                              theme(legend.text  = element_text(colour = "black", size = 8) )+
                              theme(axis.text.y  = element_text(colour = "black", size = 6))+
                              theme(axis.text.x  = element_text(colour = "black", size = 6))+
                              theme(plot.margin  = unit(c(0.5,0.2,0.5,0.2),"cm")  )+
                              #scale_color_discrete(name=titles[1])+
                              #scale_linetype_manual(name=titles[1],values=c(1,5))+
                              scale_x_discrete(limits=unique(Df$Time))+
                              theme(legend.key.width = unit(3,"line")) 
                            #scale_y_continuous(breaks=c(seq(0,10,by=2)) )
      )
    }
    
    
    if (smooth){
      if ( (dim(Df)[1]<30000) && (length(unique(Df$Time))>3)){
        old_plot = plot_list[[count]] 
        res <- try( list(plot_list[[count]] + 
                           stat_smooth(method="loess",
                                       fullrange=F,
                                       size=0.5,
                                       span=1,
                                       aes(linetype=Type,color=Type)) ),
                    silent = T)
        if (inherits(res,"try-error")){
          plot_list[count]=list(old_plot +
                                  stat_smooth(method='glm',
                                              formula = y ~ poly(x, length(unique(Df$Time))-1),
                                              fullrange=F,
                                              size=0.5,
                                              aes(linetype=Type,color=Type)) 
          )
          # plot_list[count]=list(plot_list[[count]] + 
          # stat_smooth(method="loess", span = 0.1 ,fullrange=F, size=0.5 ,aes(linetype=Type,color=Type)) 
          # )
        }else{
          plot_list[count]=res
          
        }
      }else{
        plot_list[count]=list(plot_list[[count]]+
                                stat_smooth(method='glm',
                                            formula = y ~ poly(x, length(unique(Df$Time))-1),
                                            fullrange=F,
                                            size=0.5,
                                            aes(linetype=Type,color=Type)) 
        )
      }
      #below is Liron's code which gave bad results for me
      # gp2=plot_list[count]
      # res<-try(get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1]))) ,silent = T)
      # if (inherits(res,"try-error")){
      #   plot_list[count]=list(plot_list[[count]]+
      #                           stat_summary(fun.y=mean, geom="point", size = 3,aes(color=Type))+
      #                           stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
      #   )
      # }
    }else{
      plot_list[count]=list(plot_list[[count]]+
                              stat_summary(fun.y=mean, geom="line", size = 1.3,aes(linetype=Type,color=Type))+
                              stat_summary(fun.data = "mean_se", geom = "errorbar", width = .3,size = 1,aes(color=Type),show.legend = F)
      )
    }
    
    out_file1 = str_replace(out_file, ".pdf", paste0("_cluster", i, ".png"))
    ggsave(out_file1,plot_list[[count]],dpi = 600, width = 3, height =3, units = "in")
    
    count=count+1  
  }
  
  gp2=plot_list[1]
  num_of_plots=length(plot_list)
  if (num_of_plots <= 2) {
    grob_rows=1
    grob_cols=num_of_plots
  } else if (num_of_plots <= 4) {
    grob_rows = grob_cols = 2
  } else {
    grob_rows=2
    grob_cols=3
  }
  if (color_group[1]==color_group[2]){
    Clustering_Plot<-marrangeGrob(plot_list, nrow=grob_rows, ncol=grob_cols, top=NULL)
  } else {
    legend <- get_legend(gp2[[1]] + theme(legend.position  ="right")+guides(color=guide_legend(title=titles[1]),linetype=guide_legend(title=titles[1])))
    Clustering_Plot<-marrangeGrob(plot_list, nrow=grob_rows, ncol=grob_cols,right=legend, top=NULL)
  }
  
  
  #Clustering_Plot = clusters_plots[[1]]
  
  ggsave(out_file,Clustering_Plot,dpi = 600, width = 6.99, height =3.99, units = "in")
  
  #export data behind plot
  
  out_file2 = str_replace(out_file, ".pdf", "_data.txt")
  
  for (i in sort(unique(clusters))){
    out_file3 = str_replace(out_file, ".pdf", paste0("_cluster", i, "_data.txt"))
    Df1 = plot_list[[i]]$data %>% select(Group = Time, Exp = EXP)
    export_table(Df1, out_file3)  #optional
    se_by_group    = Df1 %>% group_by(Group) %>% summarise(Mean = mean(Exp), SE = sd(Exp)/sqrt(length(Exp)), Mean_SE = mean_se(Exp))%>% as.data.frame
    se_by_group = cbind (Cluster = i, se_by_group)
    if (i==1) {
      summary_df = se_by_group
    } else {
      summary_df = rbind(summary_df, se_by_group)
    }
  }
  
  write.table(x = summary_df,
              file = out_file2,
              quote = F,
              sep = "\t",
              na = "",
              row.names = F,
              col.names = T)
  
  #return(list(Clustering_Plot))
  return(Clustering_Plot)
}

##### Utilities #####

compute_group_means = function (mat2plot, col_data) {
  
  #compute group means
  #assuming that col order in mat2plot is the same as row order in col_data
  
  mat2plot_means = mat2plot %>%
    t %>%
    aggregate.data.frame(by = col_data[GROUP],FUN = mean) %>%
    remove_rownames %>%
    column_to_rownames(var=GROUP) %>%
    as.data.frame() %>%
    t
  return (mat2plot_means)
}

prepare_generic_df_for_manual_clustering = function (mat2plot, col_data) {
  
  mat2plot_means = compute_group_means (mat2plot, col_data)
  mat2plot_means_scaled <- z_score(mat2plot_means)
  
  new_df = mat2plot_means_scaled[,COLUMNS_FOR_MAN_CLUSTERING] #take first 3 columns of mat2plot_means (you may change it if needed)
  colnames(new_df) = c("A", "B", "C")
  
  new_df = 
    new_df %>%
    as.data.frame %>%
    rownames_to_column(var="gene") %>%
    mutate ("B.vs.A" = B - A) %>%
    mutate ("C.vs.B" = C - B) 
  
  return (new_df)
}

get_average_expression_per_cluster = function (clusters, mat2plot, manual_clusters_dir) {
  
  #The function is not yet ready.currently it just prints data to file for further manipulation in excel
  
  for (num in sort(unique(clusters))){  #for each cluster
    #get a list of the genes in this cluster
    Genes = get_genes_per_cluster(clusters, num)
    
    if (length(Genes)>1){
      expression_matrix = mat2plot[Genes,]
      out_file = paste0 (manual_clusters_dir, "/", "cluster_", num, "_expression.csv")
      write.csv(expression_matrix, file=out_file)
    }
  }
}

create_clusters_object = function (genes, clusters) {
  #receives 2 matching vectors 
  #you may wish to do as.numeric(clusters) before submitting to the function
  clusters <- setNames(clusters, genes)
  return (clusters)
}

make_data_frame_from_clusters = function (clusters) {
  if (length(clusters)> 0) {
    clusters_df = data.frame(gene=names(clusters), cluster=clusters)
  } else {
    clusters_df = data.frame(matrix(nrow=0, ncol=2)) %>% setNames(c("gene", "cluster")) #empty df
  }
  return (clusters_df)
}

make_clusters_from_data_frame = function (df) {
  #clusters <- setNames(as.numeric(df$cluster), df$gene)  #not sure we need as.numeric
  clusters <- setNames(df$cluster, df$gene)
  return(clusters)
}

get_genes_per_cluster = function (clusters, cluster_name) {
  return (names(clusters[clusters == cluster_name]))
}

direction = function (FC, cutoff) {
  
  #I haven't seen any place in the code where this function is used
  
  if (FC>cutoff) {
    return ("up")
  } else if (FC<-cutoff) {
    return ("down")
  } else {
    return ("nc")
  }
}