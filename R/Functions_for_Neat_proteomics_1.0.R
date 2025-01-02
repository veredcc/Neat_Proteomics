## Functions for Analysis of ProteinGroups results file from MaxQuant
## author: "Vered Chalifa-Caspi"
## Modified from: "Tony Lin, UCSF, https://github.com/YHTLin/DataScience-tutorial/tree/master/proteomics-series"

# Data Imputation Explanation
# Since missing values are associated with proteins with low levels of expression, we can substitute the missing values
# with numbers that are considered ?small? in each sample. We can define this statistically by drawing from a normal
# distribution with a mean that is down-shifted from the sample mean and a standard deviation that is a fraction of the
# standard deviation of the sample distribution.
# read more: Lazar et al J. Proteome Res 2016 https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00981
# 
# Vered: the defaults he puts in the function arguments, width = 0.3, downshift = 1.8, are the same as those in Perseus
# Vered: he uses set.seed(1) I checked it, and this causes the imputation to produce the same imputed data every run
# I commented it out in order to be able to run the imputation several times and each time get different results


##### Preparations #####

import_data = function(data_file) {
  
  #read file into data frame
  data = read.delim(data_file, stringsAsFactors = FALSE, colClasses = "character")
  
  return(data)
}

import_samples = function(samples_file) {
  #read file into data frame
  samples = read.delim(samples_file, stringsAsFactors = FALSE, colClasses = "character")
  #make sample names valid names (otherwise they may not fit the sample names that were read from the proteinGroups file)
  samples = data.frame(lapply(X=samples, FUN=make.names))
  
  return(samples)
}

import_col_data = function(experiment_design_file) {
  
  #read file into data frame
  design = read.delim(experiment_design_file, stringsAsFactors = TRUE, colClasses = "factor")
  
  #re-factoring all expr_design factors
  design = data.frame(lapply(X = design, FUN = factor))
  
  #set row names to SampleID               #new Vered 11.2.2024
  rownames(design) <- design$SampleID
  
  return (design)
}

##Gil

# import_go_pathways = function(GO_pathways_file) {
#   Proteins2GOs = read.delim(file = GO_pathways_file,sep = '\t',col.names = c('Gene','GO'))
#   
#   return (Proteins2GOs)
# }

clean_raw_data = function(raw_data) {
  
  # Filter false hits
  # Potential.contaminant: Yishai does not filter them out but if I retain them
  #there is a problem with parsing protein names because there is no _HUMAN suffix to some of them
  # Only.identified.by.site: Keren does not remove these proteins, but Gideon said we should remove
  # them. See also in the presentation of Gideon project 06. 
  df = raw_data %>% 
    filter_at(vars(Potential.contaminant, Reverse, Only.identified.by.site), all_vars(. != "+"))
  
  # Summary of Q values #Q values indicate the confidence in protein identification.
  # It represents the probability that the protein is a false hit
  # A typical cutoff is set at 0.01. Fortunately, MaxQuant takes care of this operation
  # and ensures that all Q values are below the threshold.
  
  dim(df)  
  summary(as.numeric(df$Q.value))  #according to Tony Lin, a typical cutoff for Q is set at 0.01.
  
  return (df)
  
}

#extract_protein_name = function(df, reg_expr = "(?<=\\s).*(?=.OS)") {
#  
#  # default reg_expr starts from the first space ('\\s') until 'OS'
#  
#  # we use regular expression to extract the protein names
#  regex = regexpr(reg_expr, df$Fasta.headers, perl = TRUE)
#  
#  # return the matching header for each protein
#  return (regmatches(df$Fasta.headers, regex))
#}

extract_protein_name <- function(df, REF_DB) {  #new, 26.12.2024. need to be tested on uniprot and genbank
  
  # Define regex patterns for each type
  regex_patterns <- list(
    uniprot = "(?<=\\s).*(?=.OS)",   # Try also (?<=\\s).*(?=\\sOS=), perhaps it is better
    genbank = "(?<=\\|)[^|]+$",      # NOT YET TESTED Adjust if needed for GenBank format specifics
    refseq = "(?<=\\s).*?(?=\\s\\[)" # Extract text after first space until before the species info in brackets
  )
  
  # Select the appropriate regex based on the fasta_title_type
  if (!REF_DB %in% names(regex_patterns)) {
    stop("Invalid REF_DB. Choose from: 'uniprot', 'genbank', 'refseq'.")
  }
  
  reg_expr <- regex_patterns[[REF_DB]]
  
  # Apply the regex to extract the protein name
  regex <- regexpr(reg_expr, df$Fasta.headers, perl = TRUE)
  matches <- regmatches(df$Fasta.headers, regex)
  
  # Return the matches
  return(matches)
}




#get measurments (intensity, LFQ, nr. peptides) per sample

get_measurements_per_sample = function (df, measurement_prefix) {
  
  #measurement_prefix = "LFQ.intensity."
  
  #store Interanl.ID column
  #ids = pull(df, Internal.ID)
  #
  
  #create a data frame of the required measurements
  df_m = dplyr::select (df, starts_with(measurement_prefix))
  df_m = setNames(df_m, str_replace (names(df_m), pattern = measurement_prefix, replacement = ""))
  
  #log2 transformation
  if (grepl("intensity", measurement_prefix, ignore.case=T)) {
    df_m = as.data.frame(log2(sapply(df_m, as.numeric)), row.names = row.names(df_m))
  }
  
  #set new sample names (according to samples data frame)
  
  #remove leading X from sample orig names (if exists). These Xs are added during import of the Samples.txt file if sample names start with digits.
  sample_orig_names = samples$Sample_orig_name
  if (all(grepl("^X", sample_orig_names))) {
    # Remove the leading "X" from each value
    sample_orig_names <- sub("^X", "", sample_orig_names)
  }
  
  #set new sample names
  colnames(df_m) = samples$SampleID[match(colnames(df_m), sample_orig_names)]
  
  #reorder samples according to experiment design (according to $sampleID column in experiment design data frame)
  
  #add leading Xs to SampleID from col_data if df_m col names start with Xs
  
  col_data_sampleID = make.names(col_data$SampleID)
  
  #now reorder samples
  df_m1 = df_m[,as.character(col_data_sampleID)]
  
  return (df_m1)
}

remove_samples_from_col_data = function (col_data, outlier_samples) {
  col_data[col_data$SampleID %in% setdiff(col_data$SampleID, outlier_samples),]
}

remove_samples_from_dataset = function (df, outlier_samples) {
  select(df, -contains(outlier_samples))
}

match_df_to_col_data = function (df, col_data) { ###Vered 8.2.2022 The code is from the DESeq2 script. Need to adjust it to the proteomics!!!
  for(tbl in c("abundance", "counts", "length")) {
    df[[tbl]] <- df[[tbl]][,col_data$SampleID]
  }
  
  # Check col_data and df agree in column headers:
  for(tbl in c("abundance", "counts", "length")) {
    print(tbl)
    identical(rownames(col_data),   # SampleID order in col_data  ###need to fix. exr_design row names are simple numbers
              colnames(df[[tbl]])) %>% print
  }
  
  return (df)
}

##### Statistical analysis #####

# Generate design model for limma
generate_design_model_Bad = function(col_data, APPLY_BATCH_CORR) {  ###Vered: the manipulations may cause a wrong design matrix, prob. because of wrong automatic ordering of the cond levels
  
  # create the condition factors vector for the model
  cond_factors = factor(as.character(col_data[,"Condition"]))
  
  # create the batch factors vector for the model
  if (isTRUE(APPLY_BATCH_CORR))
    batch_factors = factor(as.character(col_data[,BATCH_FACTOR]))
  
  # create the matrix model
  if (APPLY_BATCH_CORR == FALSE) {
    design <- model.matrix(~ 0 + cond_factors)
    colnames(design) <- unique(col_data[,"Condition"])
  } else {
    design <- model.matrix(~ 0 + cond_factors + batch_factors)  #Vered 3.2.2022
    cols = colnames(design)                                     #need to find
    cols = sub("batch_factors","", cols)                        #a more elegant way
    cols = sub("cond_factors", "", cols)                        #to do the col name modification
    colnames(design) <- cols
  }
  
  return (design)
}

# Generate design model for limma  ###vered: here I go back to my original code. should be corrected, to use the experiment design file
generate_design_model_BAD1 = function(samples_to_groups, conditions, APPLY_BATCH_CORR) {
  if (APPLY_BATCH_CORR == FALSE) {
    design <- model.matrix(~ 0+factor(c(samples_to_groups)))  
    colnames(design) <- conditions
  } else {
    # batch_factors = as.character(col_data[,BATCH_FACTOR])
    # design <- model.matrix(~ 0+factor(c(samples_to_groups))+factor(batch_factors))  
    # colnames(design) <- conditions  #Vered: gives error!
    
    batch = factor(as.character(col_data[,BATCH_FACTOR]))
    cond  = factor(as.character(col_data[,"Condition"]), levels=conditions)
    cond  = factor((c(samples_to_groups)))
    design <- model.matrix(~ 0+cond+batch)  
    colnames(design) <- conditions  #Vered: gives error!
  }
  
  return (design)
}

# Generate design model for limma 
generate_design_model = function(col_data, DESIGN_TYPE = 'simple') {
  
  #conditions are defined in the main program. the levels parameters makes sure that the levels are in the order they
  #appear in the experimental design, which is supposed to be the same as in the proteinGroups file.
  #Otherwise, the design matrix may become wrong.
  
  cond  = factor(as.character(col_data[,"Condition"]), levels=conditions) 
 
  if (DESIGN_TYPE %in% c('simple', 'batch-multi_level')) {
    
    design <- model.matrix(~ 0+cond)  
    colnames(design) <- conditions
        
  } else if (DESIGN_TYPE == 'batch') {
    
    batch_values = col_data[,BATCH_FACTOR] %>%
      as.character %>%
      make.names
    batch = factor(batch_values, levels=batches)  #batches are defined in the main program
    
    design <- model.matrix(~ 0+cond+batch)  
    
    cols = colnames(design)
    cols = gsub('cond|batch', '', cols)
    colnames(design) <- cols

  }
   
  return (design)
}

## Filter data for "valid" (expressed) proteins only
filter_valids = function(df, conditions, min_count, Log2.names) {
  
  # df         = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping, e.g. c("C", "S", "SH")
  # min_count  = a numeric vector of the same length as "conditions" indicating the minimum 
  #              number of valid values for each condition for retention, e.g. c(3, 3, 3)
  
  cond.filter = sapply(1:length(Log2.names), function(i) { #for each condition (e.g. C, S, SH)
    df2 = df[Log2.names[[i]]]      # Extract columns of interest
    df2 = as.matrix(df2)           # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })

  df$KEEP = apply(cond.filter, 1, any)
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}

# Generate the contrasts matrix
generate_contrasts = function(contrasts_data, design) {
  
  #create contrasts vector
  contrasts <- vector()
  contrast.names <- contrasts_data[,'Contrast_name']
  
  #loop all contrasts from the file
  for (i in 1:nrow(contrasts_data))
    #contrast.names <- c(contrast.names, contrasts_data[i,'Contrast_name'])
    contrasts <- c(contrasts, sprintf("%s - %s",contrasts_data[i,'Numerator'],contrasts_data[i,'Denominator']))
  
  #create the contrasts matrix
  contrast.matrix <- makeContrasts(contrasts=contrasts, levels=design)
  colnames(contrast.matrix) <- contrast.names
  
  return (contrast.matrix)
}

# present-absent analysis
count_nonzero_values = function(df, conditions, Log2.names) {
  
  present.absent = dplyr::select(df)  #take just row names
  
  #count number of non-zero values per condition
  for (condition in conditions) {
    cols.per.cond = Log2.names[[match(condition,conditions)]]
    
    col_name_sum = paste ("No.nonzero.values.", condition, sep="")
    #col_name_pa  = paste ("Present.Absent.",    condition, sep="")
    
    pa_df = df %>%
      as_tibble %>%
      dplyr::select(all_of(cols.per.cond)) %>%
      replace_with_na_all (condition = ~.x == -Inf) %>%  #replace -Inf with NA.
      mutate(!!col_name_sum := rowSums(!is.na(.))) %>%
      as.data.frame
    present.absent = cbind(present.absent, pa_df)
  }
  
  #leave only count columns
  present.absent = present.absent %>%
    as_tibble %>%
    dplyr::select(starts_with("No.nonzero.values")) %>%
    as.data.frame
  
  return (present.absent)
}

# up-down analysis
analyze_updown = function(present.absent, contrasts) {
  up.down = dplyr::select(present.absent)  #take just row names
  
  for (contrast in contrasts) {
    contrast_print = gsub("\ ", "", contrast) #for column headers remove spaces, e.g. instead of SH - C, write SH-C
    conds = unlist(strsplit(contrast_print, "-", fixed=T))  #split to conditions
    conds = unlist(strsplit(contrast_print, ".vs.", fixed=T))  #split to conditions
    cond1 = conds[2]
    cond2 = conds[1]
    cond1_col = paste ("No.nonzero.values.", cond1, sep="")
    cond2_col = paste ("No.nonzero.values.", cond2, sep="")
    cond1_col <- rlang::sym(cond1_col)
    cond2_col <- rlang::sym(cond2_col)
    
    
    col_name_pa_updown = paste ("pa_upDown.", contrast_print, sep="")
    col_name_pa_updown <- rlang::sym(col_name_pa_updown)
    col_name_pa_pass = paste ("pa_pass.", contrast_print, sep="")
    
    pa_df = present.absent %>%
      as_tibble %>%
      mutate(!!col_name_pa_updown := ifelse (test=  !!cond1_col <= MAX_COUNT_FOR_ABSENT & !!cond2_col >=MIN_COUNT_FOR_PRESENT, yes = "pa_up", no = 
                                               ifelse(test = !!cond2_col <= MAX_COUNT_FOR_ABSENT & !!cond1_col >=MIN_COUNT_FOR_PRESENT, yes = "pa_down", no = NA))) %>%
      mutate(!!col_name_pa_pass := ifelse (!is.na(!!col_name_pa_updown), yes = 1, no = NA)) %>%
      dplyr::select(starts_with("pa_")) %>%
      as.data.frame
    up.down = cbind(up.down, pa_df)
  }
  
  return (up.down)
}

# Perform Limma analysis
run_Limma = function(df.FI, design, contrast.matrix, DESIGN_TYPE='simple') {
  
  #creating a data matrix of the relevant expression data for Limma
  Log2.data = df.FI

  if (DESIGN_TYPE == 'batch-multi_level') {
    
    batch = factor(as.character(col_data[,BATCH_FACTOR]), levels=batches)
    
    corfit <- duplicateCorrelation(Log2.data, design, block=batch)
    corfit$consensus
    fit <- lmFit(Log2.data, design, block=batch, correlation=corfit$consensus)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    #topTable(fit2, coef="DiseasedvsNormalForTissueA")
    
  } else {
    
    #fit a linear model
    fit <- lmFit(Log2.data, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
  }
  
  return (fit2)
}

compute_DE = function(fit2) {
  
  # perform DE testing for all contrasts
  if (FDR_ADJ) { 
    results <- decideTests(fit2, adjust.method="BH", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  } else {
    results <-decideTests(fit2, adjust.method="none", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF)) 
  }
  
  #print summary of DE genes to screen
  print(summary(results))
  
  # generate final results table for all contrasts
  stats_df0 = output_contrasts_results(df.FI, results, fit2, FDR_ADJ)
  return (stats_df0)
} 

# Perform Limma analysis
run_Limma_wo_contrasts = function(df, design, FDR_ADJ) {  #vered 3.2.2022 STILL NEED to check and generalize this

  #for example, for this model, which does not have 0 in the beginning, and everything is compared to the "intercept":
  #design <- model.matrix(~ batch_factors + cond_factors)
    
  #creating a data matrix of the relevant expression data for Limma
  Log2.data = dplyr::select(df, unlist(Log2.names))
  
  #fit a linear model
  fit <- lmFit(Log2.data, design)
  fit2 <- eBayes(fit)
  
  # perform DE testing for all contrasts - Vered 3.2.2022 - NOT GOOD for wo contrasts
  # if (FDR_ADJ) { 
  #   results <- decideTests(fit2, adjust.method="BH", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  # } else {
  #   results <-decideTests(fit2, adjust.method="none", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  # }
  
  #Instead, in the specific example of Dima's data:  ###Vered 3.2.2022 Need to make this general!!!
  if (FDR_ADJ) {
    results = topTable(fit2, coef="NAM7mM", number=Inf, adjust.method="BH", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  } else {
    results = topTable(fit2, coef="NAM7mM", number=Inf, adjust.method="none", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  }
  
  #print summary of DE genes to screen
  print(summary(results))
  
  #Print results to a text file  ###Vered 6.2.2022 Modify to be more general!!!
  results_file = paste (results_dir, "/datasets/topTable_Results.txt", sep="")
  export_table(results, "Results_NAM7mM_1/Datasets/topTable_Results.txt", "")
  
  # generate final results table for all contrasts  ###Vered 3.2.2022 This still does not work for wo contrasts
  stats_df0 = output_contrasts_results(df, results, fit2, FDR_ADJ)
  return (stats_df0)
} 

# Extract the contrasts results
output_contrasts_results = function(df.FI, results, fit2, FDR_ADJ) {
  
  stats_df0 = df.FI
  
  #Loop over topTable for each contrast, and add to results matrix
  for (contrast in contrasts) {
    contrast_print = gsub("\ ", "", contrast) #for column headers remove spaces, e.g. instead of SH - C, write SH-C
    stat_res = topTable(fit2, coef=contrast, adjust="BH", sort.by="none", number=Inf)
    stat_res = stat_res %>%
      as_tibble %>%
      mutate(linearFC = ifelse(is.na(logFC),
                               yes = NA,
                               no = ifelse(logFC>0,
                                           yes = 2^logFC,
                                           no = -1/(2^logFC))%>% 
                                 signif(digits = 3))) %>%
      mutate(linearRatio = ifelse(is.na(logFC),    #necessary for averaging over imputations
                                  yes = NA,
                                  no = 2^logFC%>% 
                                    signif(digits = 5))) %>%
      mutate(pvalue = ifelse(test = is.na(P.Value),
                             yes = NA,
                             no = signif(P.Value,digits = 3))) %>%
      mutate(padj   = ifelse(test = is.na(adj.P.Val),
                             yes = NA,
                             no = signif(adj.P.Val,digits = 3)))
    if (FDR_ADJ) {
      stat_res = stat_res %>% 
        mutate(pass = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                               padj <= P_CUTOFF &
                               !is.na(padj),
                             yes = ifelse(test = as.numeric(linearFC)>0,
                                          yes="up",
                                          no="down"),
                             no = "")) %>%
        mutate(pass1 = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                                padj <= P_CUTOFF &
                                !is.na(padj),
                              yes = 1,
                              no = NA))
    } else {
      stat_res = stat_res %>% 
        mutate(pass = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                               pvalue <= P_CUTOFF &
                               !is.na(pvalue),
                             yes = ifelse(test = as.numeric(linearFC)>0,
                                          yes="up",
                                          no="down"),
                             no = "")) %>%
        mutate(pass1 = ifelse(test = abs(as.numeric(linearFC)) >= LINEAR_FC_CUTOFF & 
                                pvalue <= P_CUTOFF &
                                !is.na(pvalue),
                              yes = 1,
                              no = NA))
    }
    stat_res = stat_res %>%
      mutate(manual_cutoffs = manual_cutoff_formula)
    #select(linearFC,pvalue,padj,pass,manual_cutoffs) %>%
    stats_df0 = stat_res %>%
      dplyr::select(linearFC,linearRatio, pvalue,padj,pass, pass1, manual_cutoffs) %>% 
      dplyr::rename(!!paste0("linearFC.",    contrast_print) := linearFC) %>%
      dplyr::rename(!!paste0("linearRatio.", contrast_print) := linearRatio) %>%
      dplyr::rename(!!paste0("pvalue.",      contrast_print) := pvalue) %>%
      dplyr::rename(!!paste0("padj.",        contrast_print) := padj) %>%
      dplyr::rename(!!paste0("pass.",        contrast_print) := pass) %>%
      dplyr::rename(!!paste0("pass1.",       contrast_print) := pass1) %>%
      dplyr::rename(!!paste0("manual_cutoffs.",contrast_print) := manual_cutoffs) %>%
      as.data.frame %>% 
      cbind(stats_df0,.) 
  }
  
  return(stats_df0)
}

summarize_limma_pass1_results = function(results_mult_imput) {
  
  #first, initialize a data frame to which we will add the summary of imputation results
  stat_pass1_summary_imputs = dplyr::select(results_mult_imput)  #just leave row names
  
  for (contrast in contrasts) {
    contrast_print = gsub("\ ", "", contrast) #for column headers remove spaces, e.g. instead of SH - C, write SH-C
    #pass1_col_name = paste ("pass1", contrast_print, sep=".")
    pass1_col_name = sprintf("pass1.%s.",contrast_print)
    col_name_sum = paste ("sum.pass.", contrast_print, sep="")
    col_name_pass = paste ("pass.imputs.", contrast_print, sep="")
    
    res_df = results_mult_imput %>%
      as_tibble %>%
      dplyr::select(contains(pass1_col_name)) %>%
      mutate(!!col_name_sum  := rowSums(., na.rm=TRUE)) #%>%
    res_df = res_df %>%
      mutate(!!col_name_pass := ifelse(get(col_name_sum,res_df) >= MIN_NO_PASSED, yes = 1, no = NA)) %>%
      as.data.frame
    stat_pass1_summary_imputs = cbind(stat_pass1_summary_imputs, res_df)
  }
  
  return (stat_pass1_summary_imputs)
}

summarize_limma_avgFC_imputs = function(results_mult_imput) {
  
  #first, initialize a data frame to which we will add the summary of the results
  stat_avgFC_over_imputs =  dplyr::select(results_mult_imput)  #take just row names
  
  for (contrast in contrasts) {
    contrast_print = gsub("\ ", "", contrast) #for column headers remove spaces, e.g. instead of SH - C, write SH-C
    
    Ratio_col_name = paste ("linearRatio", contrast_print, sep=".")  #col name in results_mult_imput
    pval_col_name  = paste ("pvalue",      contrast_print, sep=".")  #col name in results_mult_imput
    padj_col_name  = paste ("padj",        contrast_print, sep=".")  #col name in results_mult_imput
    col_name_Ratio = paste ("linearRatio.imputs.", contrast_print, sep="")  #col name in summary_FC_p_padj_over_imputs
    col_name_pval  = paste ("pvalue.imputs.",      contrast_print, sep="")  #col name in summary_FC_p_padj_over_imputs
    col_name_padj  = paste ("padj.imputs.",        contrast_print, sep="")  #col name in summary_FC_p_padj_over_imputs
    col_name_FC    = paste ("linearFC.imputs.", contrast_print, sep="")
    col_name_Ratio <- rlang::sym(col_name_Ratio)
    col_name_pval  <- rlang::sym(col_name_pval)
    col_name_padj  <- rlang::sym(col_name_padj)
    
    fc_df = results_mult_imput %>%
      as_tibble %>%
      select(contains(Ratio_col_name)) %>%
      mutate(!!col_name_Ratio := rowSums(., na.rm=TRUE)/NO_REPETITIONS) %>%
      mutate(!!col_name_FC := ifelse((!!col_name_Ratio)>=1,
                                     yes = (!!col_name_Ratio),
                                     no = -1/(!!col_name_Ratio))) %>%
      select(starts_with("linearFC.imputs")) %>%
      signif(digits = 4) %>%
      as.data.frame
    stat_avgFC_over_imputs = cbind(stat_avgFC_over_imputs, fc_df)
    
    pval_df = results_mult_imput %>%
      as_tibble %>%
      select(contains(pval_col_name)) %>%
      mutate(!!col_name_pval := apply (., MARGIN=1, FUN = function(x) 
        quantile(x,MIN_NO_PASSED/NO_REPETITIONS))) %>%
      #mutate(!!col_name_pval := apply (., MARGIN=1, FUN = median)) %>%
      select(starts_with("pvalue.imputs.")) %>%
      as.data.frame
    stat_avgFC_over_imputs = cbind(stat_avgFC_over_imputs, pval_df)
    
    padj_df = results_mult_imput %>%
      as_tibble %>%
      select(contains(padj_col_name)) %>%
      mutate(!!col_name_padj := apply (., MARGIN=1, FUN = function(x) 
        quantile(x,MIN_NO_PASSED/NO_REPETITIONS))) %>%
      #mutate(!!col_name_padj := apply (., MARGIN=1, FUN = median)) %>%
      select(starts_with("padj.imputs.")) %>%
      as.data.frame
    stat_avgFC_over_imputs = cbind(stat_avgFC_over_imputs, padj_df)
  }
  
  return (stat_avgFC_over_imputs)
}

summarize_upDown_limma_results = function(df) {
  stat_upDown_summary_imputs = dplyr::select(df) # take just row names
  
  for (contrast in contrasts) {
    contrast_print = gsub("\ ", "", contrast)
    FC_col_name = paste ("linearFC.imputs", contrast_print, sep=".")
    FC_col_name <- rlang::sym(FC_col_name)
    pass_col_name = paste ("pass.imputs", contrast_print, sep=".")
    pass_col_name <- rlang::sym(pass_col_name)
    updown_col_name = paste ("upDown.imputs", contrast_print, sep=".")
    
    upDown_df = df %>%
      as_tibble %>%
      mutate(!!updown_col_name := ifelse(!is.na(!!pass_col_name),
                                         yes = ifelse(test = as.numeric(!!FC_col_name)>=0,
                                                      yes="up",
                                                      no="down"),
                                         no = "")) %>%
      dplyr::select(starts_with("upDown.imputs")) %>%
      as.data.frame #%>%
    stat_upDown_summary_imputs = cbind (stat_upDown_summary_imputs, upDown_df)
  }
  
  return (stat_upDown_summary_imputs)
}

summarize_DE_proteins = function(results_mult_imput, stat_upDown_summary_imputs, stat_passany_summary_imputs) {
  DE_df = stat_upDown_summary_imputs %>% 
    mutate(dplyr::select(stat_passany_summary_imputs, pass_any_contrast))
  colnames(DE_df) = gsub("upDown.imputs", "pass", colnames(DE_df))
  return (DE_df)
}

de_summary_stats = function (DE_df, DE_proteins_stats_file) {
  c ("No. of DE proteins in any comparison: ", sum(DE_df$pass_any_contrast==1))
  DE_proteins_stats =
    apply(X = DE_df[,str_detect(string = names(DE_df),
                                pattern = "pass")],
          MARGIN = 2,
          FUN = function(x) c(up=sum(x=="up"),down=sum(x=="down"),any=sum(x!=""))) %>%
    data.frame %>%
    t
  file <- sprintf("%s_%s.txt",
                  sub(".txt","",DE_proteins_stats_file),
                  ifelse(FDR_ADJ==TRUE,
                         sprintf("PADJ_%s",P_CUTOFF),
                         sprintf("P_%s",P_CUTOFF))
                  )
  write.table (DE_proteins_stats, file=file, quote = F, sep = "\t",
               eol = "\n", na = "", dec = ".", row.names = T,
               col.names = TRUE)
  return (DE_proteins_stats)
}

create_stats_df = function (contrasts, stat_avgFC_over_imputs, stat_upDown_summary_imputs, stat_passany_summary_imputs) {
  
  #this function prepares stats_df object to be the same as in the Neat_RNA-Seq script and functions
  
  #create a dataframe with one column which is row names of stat_avgFC_over_imputs
  stats_df    = data.frame(gene=rownames(stat_avgFC_over_imputs))  
  
  #add per-contrast results
  for (contrast in contrasts) {                  
    contrast_print = gsub("\ ", "", contrast)
    stats_df = cbind(stats_df,
                     dplyr::select(stat_avgFC_over_imputs,
                                   sprintf("linearFC.imputs.%s", contrast_print),
                                   sprintf("pvalue.imputs.%s",contrast_print),
                                   sprintf("padj.imputs.%s",contrast_print)),
                     dplyr::select(stat_upDown_summary_imputs, sprintf("upDown.imputs.%s",
                                                                       contrast_print)))
  }
  
  #add pass_any_contrast per gene
  stats_df = cbind(stats_df, dplyr::select(stat_passany_summary_imputs, pass_any_contrast))
  
  #modify column names
  new_col_names = names(stats_df)                                         %>%
    str_replace (pattern = "imputs.",           replacement = "")         %>%
    str_replace (pattern = "upDown",            replacement = "pass")     %>%
    str_replace (pattern = "pass_any_contrast", replacement = "pass_any")
  stats_df = setNames(stats_df, new_col_names)
  
  return(stats_df)
}

##### Data Imputation #####

impute_data = function(df.F, width = 0.3, downshift = 1.8) {
  # df.F = data frame containing filtered LFQ values
  # Assumes missing data (in df.F) follows a narrowed and downshifted normal distribution
  
  col.names = colnames(df.F)
  impute.names = sub("^", "impute.", col.names)  #create column names "impute.C1", "impute.C2"... (columns data will be TRUE/FALSE flags)
  
  #generate an imputed dataset with imputed data and with logical values denoting whether it was imputed or not
  
  df.FI0 = df.F
  
  # Create new columns indicating whether the values are imputed
  df.FI0[impute.names] = lapply(col.names, function(x) !is.finite(df.FI0[, x]))
  
  # Imputation
  #set.seed(1)       #remove comment if you want the imputed results to be the same every run
  df.FI0[col.names] = lapply(col.names,
                           function(x) {
                             temp = df.FI0[[x]]               #Vered: for every LOG2.sample columns (I think)
                             temp[!is.finite(temp)] = NA    #Vered: if missing value, replace with NA (I think)
                             
                             temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width. Vered: na.rm tells the sd function if missing values should be removed
                             temp.mean = mean(temp, na.rm = TRUE) - 
                               downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
                             
                             n.missing = sum(is.na(temp))
                             temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                             return(temp)
                           })
  return(df.FI0)
}

##### Gene Ontology #####

unlist_pathways = function(df, index = "Gene", subject = "GO", sep = "; ") {
  l1=apply(X = df,MARGIN = 1, function(x) {
    m=as.data.frame( x = stringi::stri_split(str = x[subject],regex = sep),col.names = c("GO"))
    m["Gene"]<-x[index]
    return(m[c("Gene","GO")])
  })
  return(do.call(what = "rbind",args = l1))
}

run_enrichGO = function(DE_GO_pathways, OrgDb, GO_CUTOFF, ont = "BP") {
  go_enrich <- enrichGO(gene = DE_GO_pathways, 
                        OrgDb=OrgDb, 
                        keyType = "GO", 
                        readable = TRUE,
                        ont = ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff=GO_CUTOFF, 
                        qvalueCutoff=0.10)
  return (go_enrich)
}

createGO2Object = function(object, ont) {
  GOobject = GO2Name[GO2Name$Ontology == ont, c('GO',ifelse(object == "Gene",'Gene','Term'))]
  GOobject = GOobject[!duplicated.data.frame(GOobject),]
  
  if (!dir.exists("GO_map_files"))
    dir.create(GO_map_files)
  write.table(GOobject,
              file = file.path(getwd(),sprintf("GO_map_files/GO2%s_%s.tab",object,ont)),
              quote = F,
              row.names = F,
              sep = "\t")
  
  return (GOobject)
}

loadGO2Object = function(object, ont) {
  GOobject = read.delim(file = file.path(getwd(),sprintf("GO_map_files/GO2%s_%s.tab",object,ont)),
                        header = TRUE,
                        sep = "\t")
  return (GOobject)
}

clusters_enrichment_test = function(clusters, enricher, OrgDb, ont, Genes2GOs, 
                                    TERM2GENE, TERM2NAME, GO_CUTOFF, MIN_GENES) {
  num_of_clusters = length(sort(unique(clusters)))
  if (num_of_clusters > 1 | enricher == "enricher") {
    allRes = list()
    cluster_names = list()
    count = 1
  }
  file_name <- sprintf("%s/GO Enrichments/%s/GO_%s_%s_%s_Enrichment.csv", results_dir,
                       ifelse(enricher == "enrichGO",
                              deparse(substitute(OrgDb)),
                              "GOMap"),
                       enricher, 
                       ifelse(enricher == "enrichGO", 
                              sub("_OrgDb","",deparse(substitute(OrgDb))), 
                              "GOMap"),
                       switch(ont, 
                              "BP" = "Biological_Process", 
                              "MF" = "Molecular_Function", 
                              "CC" = "Cellular_Component"))
  for (i in sort(unique(clusters))){
    if (enricher == "enrichGO")
      enrichResult = enrichGO_enrichment_test(clusters, c(i), OrgDb, Genes2GOs, ont, MIN_GENES)
    if (enricher == "enricher")
      enrichResult = universal_enrichment_test(clusters, c(i), TERM2NAME, TERM2GENE, GO_CUTOFF, MIN_GENES)
    if (length(enrichResult) > 0 & (num_of_clusters > 1 | enricher == "enricher")){
      allRes[count] <- enrichResult
      cluster_names[count] = i
      count=count + 1
    }
  }
  if (num_of_clusters > 1 | enricher == "enricher") {
    names(allRes) = cluster_names
    allRes = clusterProfiler::merge_result(enrichResultList = allRes)
    allRes@compareClusterResult$Description <- str_replace_all(
      allRes@compareClusterResult$Description, ',', ';')
    write.csv(x = allRes@compareClusterResult,
              file = file_name,
              quote = TRUE,
              row.names = TRUE)
    allRes@fun = "enrichGO"
    return (allRes)
  } else {
    #enrichResult@result@Description <- str_replace_all(enrichResult@result@Description, ',', ';')
    write.csv(x = enrichResult@result,
              file = file_name,
              quote = TRUE,
              row.names = TRUE)
    return (enrichResult)
  }
}

enrichGO_enrichment_test <- function(clusters, num, OrgDb, Genes2GOs, ont, MIN_GENES) {
  res_red=unique(Genes2GOs[,1])
  Genes=res_red[res_red %in% names(clusters[clusters %in% num])]
  GOs = Genes2GOs %>% dplyr::filter(Gene %in% Genes) %>% dplyr::select(GO) %>% unlist
  
  res = enrichGO(gene = GOs, 
                 OrgDb = OrgDb,
                 keyType = "GO", 
                 ont = ont,
                 minGSSize     = MIN_GENES,
                 maxGSSize     = length(res_red),
                 pAdjustMethod = "BH",
                 pvalueCutoff  = GO_CUTOFF) 
  
  return(res)
}

universal_enrichment_test <- function(clusters, num, TERM2NAME, TERM2GENE, 
                                      GO_CUTOFF, MIN_GENES, pAdjustMethod='BH') {
  res_red=unique(TERM2GENE[,2])
  Genes=res_red[res_red %in% names(clusters[clusters %in% num])]
  
  res=enricher(Genes, TERM2GENE = TERM2GENE, 
               TERM2NAME = TERM2NAME,
               minGSSize     = MIN_GENES,
               maxGSSize     = length(res_red),
               pAdjustMethod = pAdjustMethod,
               pvalueCutoff  = GO_CUTOFF) 
  
  return(res)
}

simplify_enrichment_test <- function(clusters, enrichGO_object, enricher, OrgDb, ont) {
  num_of_clusters = length(sort(unique(clusters)))
  
  file_name <- sprintf("%s/GO Enrichments/%s/Simplify/GO_%s_%s_%s_Enrichment_Simplify.csv", results_dir,
                       ifelse(enricher == "enrichGO",OrgDb,"GOMap"),
                       enricher, 
                       ifelse(enricher == "enrichGO", sub("_OrgDb","",OrgDb), "GOMap"),
                       switch(ont, 
                              "BP" = "Biological_Process", 
                              "MF" = "Molecular_Function", 
                              "CC" = "Cellular_Component"))
  
  enrichResult_simp <- try(clusterProfiler::simplify(enrichGO_object, 
                                                     cutoff = 0.7, 
                                                     by = "p.adjust", 
                                                     select_fun = min,
                                                     measure = "Wang",
                                                     semData = NULL),
                           silent = T)
  if (!inherits(enrichResult_simp,"try-error")){
    if (num_of_clusters > 1 | enricher == "enricher") {
      enrichResult_simp@compareClusterResult$Description <- str_replace_all(
        enrichResult_simp@compareClusterResult$Description, ',', ';')
      write.csv(x = enrichResult_simp@compareClusterResult,
                file = file_name,
                quote = FALSE,
                row.names = TRUE)
      return (enrichResult_simp)
    } else {
      write.csv(x = enrichResult_simp@result,
                file = file_name,
                quote = TRUE,
                row.names = TRUE)
      return (enrichResult_simp)
    }
  }
  return (NULL)
}

##### Filtering #####

get_proteins_from_file = function (proteins_list_file, has_header=FALSE) {
  
  # get protein list from file, assuming proteins are in the first column
  data = read.delim(proteins_list_file, as.is=T, header=has_header)
  protein_list_from_file = data[[1]]
  
  return(protein_list_from_file)
}

filter_expression_matrix_by_protein_list = function(LFQ_counts, DE_proteins, proteins_list) {
  
  #this function produces an object with three gene expression matrices
  
  LFQ_counts = replace_NAs_with_min_counts(LFQ_counts)
  
  proteins <- rownames(LFQ_counts) %in% DE_proteins
  mats2plot <- list(LFQ_counts,LFQ_counts[proteins,])
  if (isTRUE(proteins_list)) {   #Vered 21.11.2023 this was not checked or tested
    proteins_file <- rownames(LFQ_counts) %in% protein_list_from_file
    mats2plot[[3]] <- LFQ_counts[proteins_file,]
  }
  
  return (mats2plot)
}

filter_expression_matrix_by_gene_list = function (expr_data, gene_list) {
  
  #this function is from the Neat_RNA-Seq script, and it produces one gene expression matrix
  
  #assuming row names are gene names
  genes <- rownames(expr_data) %in% gene_list
  return (expr_data[genes,])
}

replace_NAs_with_min_counts = function(LFQ_counts) {
  for (col in colnames(LFQ_counts))
    LFQ_counts[,col] <- replace(LFQ_counts[,col],
                                is.na(LFQ_counts[,col]),
                                min(LFQ_counts[!is.na(LFQ_counts[,col]),col]))
  return (LFQ_counts)
}

get_DE_genes_list_per_contrast = function (contrast, stats_df) {
  
  #create a named vector where element names are genes and values are clusters ("up", "down")
  
  clusters = c() #start with an empty vector.
  pass_col = paste0("pass.", contrast)
  up   = stats_df %>% filter (get({{pass_col}}) == "up") %>% row.names
  down = stats_df %>% filter (get({{pass_col}}) == "down") %>% row.names
  up_cluster   = setNames(rep("up",   length(up)),   up)
  down_cluster = setNames(rep("down", length(down)), down)
  clusters = c(clusters, up_cluster, down_cluster)
  
  return (clusters)
}

get_genes_with_corr_to_best_pattern = function (corrs2, pattern, sort_by_corr=T) {
  pattern_name = paste0("_", pattern)
  corrs2_f = corrs2 %>% rownames_to_column("gene") %>% filter (best_pattern == pattern_name)
  if (sort_by_corr == TRUE) {
    pattern_name <- rlang::sym(pattern_name)
    corrs2_f = arrange(corrs2_f, desc(!!pattern_name))
  }
  return(corrs2_f$gene)
}

##### Graphics #####

# Draw LFQ histograms
draw_LFQ_histograms = function(df_m, plots_dir) {
  
  result_file = file.path(plots_dir, "LFQ_histograms_summary.png")
  
  Log2.data = df_m %>% gather("Condition", "Intensity")
  Log2.data = Log2.data[is.finite(Log2.data$Intensity),]
  
  #remove leading X from Condition (if exists).
  if (all(grepl("^X", names(df_m)))) {
    # Remove the leading "X" from each value oc Condition
    Log2.data$Condition <- sub("^X", "", Log2.data$Condition)
  }
    
  # Create labels
  Log2.data = Log2.data %>%
    mutate(replicate = str_extract(Condition, ".$")) %>%
    mutate(Condition = get_conditions_from_samples(Condition))
  
  plot <- ggplot(Log2.data, aes(x = Intensity, fill=Condition)) +
    geom_histogram(alpha = 0.3, binwidth = 0.4, position = "identity") +
    labs(x = expression("log"[2]*"-transformed LFQ Intensity"), y = "Frequency") +
    ggtitle("LFQ Histograms summary") +
    facet_grid(replicate ~ Condition) + 
    theme_stata()
  ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 6) 
  
  for (cond in names(df_m)) {
    result_file = file.path(plots_dir, "LFQ_Histograms", paste0(cond, ".png"))
    cond_data = dplyr::select(df_m, cond) %>% gather("Condition", "Intensity")
    cond_data = cond_data[is.finite(cond_data$Intensity),]
    plot <- ggplot(cond_data, aes(x = Intensity)) +
      geom_histogram(alpha = 0.3, binwidth = 0.4, position = "identity", fill="green", color="black") +
      labs(x = expression("log"[2]*"-transformed LFQ Intensity"), y = "Frequency") +
      ggtitle(sprintf("%s Histogram", cond)) +
      theme_stata()
    ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 6)
  }
}

# Draw LFQ boxplots
draw_LFQ_boxplots = function(df_m, plots_dir) {
  
  result_file = file.path(plots_dir, "LFQ_boxplots.png")
  
  Log2.data = df_m %>% gather("Condition", "Intensity")
  Log2.data = Log2.data[is.finite(Log2.data$Intensity),]
  
  #remove leading X from Condition (if exists).
  if (all(grepl("^X", names(df_m)))) {
    # Remove the leading "X" from each value oc Condition
    Log2.data$Condition <- sub("^X", "", Log2.data$Condition)
  }
  
  Log2.data = mutate(Log2.data, Condition_name = Condition) %>%
    mutate(Condition = get_conditions_from_samples(Condition_name))
  
  plot <- ggplot(Log2.data, aes(Condition_name, Intensity)) + 
    geom_boxplot(aes(color=Condition), coef=5) +
    ggtitle("LFQ Boxplots") +
    xlab("Condition") +
    ylab(expression("log"[2]*"(Intensity)")) +
    theme_stata()
  
  ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 6)
}

# Draw distribution of imputed data
draw_imputed_histograms = function(df.FI0, plots_dir) {
  
  result_file = file.path(plots_dir, "imputed_histograms_summary.png")
  
  LOG2.df   = select (df.FI0, -starts_with("impute"))
  impute.df = select (df.FI0, starts_with("impute"))
  
  # Reshape data into key-value pairs
  LOG2.df = gather(LOG2.df, "sample", "intensity")
  impute.df = gather(impute.df, "sample", "impute")
  
  # Combine data
  combine.df = bind_cols(LOG2.df, impute.df["impute"])
  
  # Create labels
  combine.df = combine.df %>%
    mutate(replicate = str_extract(sample, ".$")) %>%
    mutate(sample = sub(".$", "", sample))
  
  plot = ggplot(combine.df, aes(x = intensity, fill = impute)) +
    geom_histogram(alpha = 0.3, binwidth = 0.4, position = "identity") +
    labs(x = expression("log"[2]*"-transformed LFQ Intensity"), y = "Frequency") +
    facet_grid(replicate ~ sample) +
    ggtitle(sprintf("Imputed data histograms summary (WIDTH=%s, DOWNSHIFT=%s)", WIDTH, DOWNSHIFT)) +
    theme_stata() +
    scale_fill_discrete(name = "Imputed",
                        breaks = c("FALSE", "TRUE"),
                        labels = c("-", "+"))
  
  ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 7)
}

draw_sample_correlation_matrix = function (drawing_df, result_file) {
  
  sampleDists <- dist(t(drawing_df))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  png(filename = result_file)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
}

perform_pca = function (drawing_df, plots_dir, file_extension="") {
  
  #compute PCA
  vsd_pca <- prcomp(t(drawing_df))                      #compute PCA
  pca_percentVar <- vsd_pca$sdev^2/sum(vsd_pca$sdev^2)  #compute percent variation in each PC
  names(pca_percentVar) <- colnames(vsd_pca$x)
  
  #draw_2D_plots
  result_file = file.path(plots_dir, paste0("PCA", file_extension))
  draw_pca (vsd_pca, pca_percentVar, result_file, 1, 2) #draw PCA 1, 2
  draw_pca (vsd_pca, pca_percentVar, result_file, 1, 3) #draw PCA 1, 3
  
  #draw_3D_plot
  result_file = paste0("PCA_3D_", analysis_name, file_extension, ".html")
  draw_pca_3d (vsd_pca, pca_percentVar, result_file)
}

draw_pca = function (vsd_pca, pca_percentVar, result_file, pc_x=1, pc_y=2) {
  
  axis1 = paste0("PC", pc_x)
  axis2 = paste0("PC", pc_y)
    
  color_factor = EFFECTS[1]            
  shape_factor = EFFECTS[2]
  
  result_file1 = paste0(result_file, "_", axis1, ".vs.",axis2, ".png")
  
  pca2plot <- data.frame(vsd_pca$x[,c(axis1,axis2)],
                         Color=col_data[,color_factor],
                         Sample=col_data$SampleID)
  if(!is.na(shape_factor))
    pca2plot <- pca2plot %>% mutate(Shape=col_data[,shape_factor])
  
  png(filename = result_file1)
  
  if(!is.na(shape_factor)) {
    plot = ggplot(pca2plot,aes(x=get(axis1),
                               y=get(axis2),
                               col=Color,
                               shape=Shape))
  } else {
    plot = ggplot(pca2plot,aes(x=get(axis1),
                               y=get(axis2),
                               col=Color))
  }
  plot = plot + geom_point(size=5) +
    xlab(sprintf("%s: %1.2f%% of variance", axis1, pca_percentVar[axis1]*100)) + 
    ylab(sprintf("%s: %1.2f%% of variance", axis2, pca_percentVar[axis2]*100)) + 
    ggtitle(sprintf("PCA: %s.vs.%s",axis1,axis2)) +
    theme_economist() +
    coord_fixed()
  print(plot)
  dev.off()
}

draw_pca_3d = function (vsd_pca, pca_percentVar, result_file) {

  color_factor = EFFECTS[1]
  shape_factor = EFFECTS[2]
  
  pca2plot3d <- data.frame(vsd_pca$x[,c("PC1","PC2","PC3")],
                           Color=col_data[,color_factor],
                           Sample=col_data$SampleID)
  if(!is.na(shape_factor))
    pca2plot3d <- pca2plot3d %>% mutate(Shape=col_data[,shape_factor])
  
  symbols2use <- c('diamond','triangle-down','square','circle','x','o')
  
  if(!is.na(shape_factor)) {
    p<-plot_ly(pca2plot3d) %>% 
      add_markers(x = ~PC1, y = ~PC2, z = ~PC3, 
                  colors = rainbow(4),     #c("red","blue"),
                  text=~Sample, 
                  color = ~Color,
                  symbol = ~Shape,
                  symbols = symbols2use)
  } else {
    p<-plot_ly(pca2plot3d) %>% 
      add_markers(x = ~PC1, y = ~PC2, z = ~PC3, 
                  colors = rainbow(4),     #c("red","blue"),
                  text=~Sample, 
                  color = ~Color)
  }
  p <- p %>% layout(scene = list(xaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC1",
                                                              pca_percentVar["PC1"]*100)),
                                 yaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC2",
                                                              pca_percentVar["PC2"]*100)),
                                 zaxis = list(title = sprintf("%s: %1.2f%% of variance", "PC3",
                                                              pca_percentVar["PC3"]*100))))
  
  #saveWidget(p, file = file.path(getwd(),result_file), selfcontained = T, title=sprintf("Principal Component Analysis"))  #Vered 6.2.2022 This does not work
  saveWidget(p, file = result_file, selfcontained = T, title=sprintf("Principal Component Analysis"))                    #Vered 6.2.2022 Need to place the file in an internal dir !!!
}

draw_ma_plots = function (fit2, plots_dir) {
  
  # perform DE testing for all contrasts
  if (FDR_ADJ) { 
    results <- decideTests(fit2, adjust.method="BH", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  } else {
    results <-decideTests(fit2, adjust.method="none", p.value = P_CUTOFF, lfc = log2(LINEAR_FC_CUTOFF))
  }
  
  for (i in 1:length(contrasts)) {
    contrast = contrasts[i]
    MA_Plot_file = file.path(plots_dir, paste0("MA_Plot_", contrast, ".png"))
    png(filename = MA_Plot_file)
    plotMD(fit2,coef=i,
           status=results[,i],values=c(1,-1),
           hl.col="blue", hl.cex=0.5,
           bg.col="darkgray", bg.cex=0.5, 
           main=sprintf("MA-Plot: %s (Fold Change cutoff: %s)",contrast, LINEAR_FC_CUTOFF),
           xlab="Mean of normalized counts",
           ylab="log fold change",
           legend = FALSE)
    #ylim=range(topTable(fit2, coef=contrast, adjust="BH", sort.by="none", number=Inf)$logFC, na.rm=TRUE))
    abline(h = 0, col = "dimgrey", lty=1, lwd=5)
    abline(h = c(LINEAR_FC_CUTOFF,-LINEAR_FC_CUTOFF),col="green")
    
    dev.off()
  }
}

draw_volcano_plots = function (fit2, plots_dir) {
  
  for (contrast in contrasts) {
    
    data2plot = topTable(fit2, coef=contrast, adjust="BH", sort.by="none", number=Inf) %>%
      rownames_to_column(var="protein")
    if (FDR_ADJ) {
      data2plot = data2plot %>%
        mutate(highFC=(abs(logFC)>log2(LINEAR_FC_CUTOFF)),
               lowPADJ=(adj.P.Val<=P_CUTOFF),
               signif=((abs(logFC)>log2(LINEAR_FC_CUTOFF)) & (adj.P.Val<=P_CUTOFF)),
               neglog10=-log10(P.Value),
               ratio=2^logFC,
               FC=ifelse(logFC>0,
                         2^logFC,
                         -1/2^logFC)) 
    } else {
      data2plot = data2plot %>%
        mutate(highFC=(abs(logFC)>log2(LINEAR_FC_CUTOFF)),
               lowP=(P.Value<=P_CUTOFF),
               signif=((abs(logFC)>log2(LINEAR_FC_CUTOFF)) & (P.Value<=P_CUTOFF)),
               neglog10=-log10(P.Value),
               ratio=2^logFC,
               FC=ifelse(logFC>0,
                         2^logFC,
                         -1/2^logFC)) 
    }
       
    volcano_Plot_file = file.path(plots_dir, paste0("Volcano_Plot_", contrast, ".png"))
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=logFC, y=neglog10)) +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(log[2]~fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-log2(LINEAR_FC_CUTOFF), log2(LINEAR_FC_CUTOFF)),color = "blue", lty = 2) +
      geom_abline(slope = 0, intercept = ifelse(FDR_ADJ,-log10(P_CUTOFF),-log10(P_CUTOFF)),color = "blue", lty = 2) +
      ggtitle(sprintf("Volcano-Plot: %s (%s = %s)",contrast, 
                      ifelse(FDR_ADJ,"PADJ cutoff","p-value cutoff"),
                      ifelse(FDR_ADJ,P_CUTOFF,P_CUTOFF)))
    #volcano_plot + geom_label(label = "Look at this", x=0, y=6)
    
    ggsave(plot = volcano_plot, device="png",
           filename = volcano_Plot_file, width = 7, height = 7)
 
    #print volcano plot with point labels on DE genes
    
    volcano_Plot_file1 = file.path(plots_dir, paste0("Volcano_Plot_", contrast, "_w_labels.png"))
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=logFC, y=neglog10, label=protein)) +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      geom_text(data = subset(data2plot, signif == T), hjust = 0, nudge_x = 0.15, size = 2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(log[2]~fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-log2(LINEAR_FC_CUTOFF), log2(LINEAR_FC_CUTOFF)),color = "blue", lty = 2) +
      geom_abline(slope = 0, intercept = ifelse(FDR_ADJ,-log10(P_CUTOFF),-log10(P_CUTOFF)),color = "blue", lty = 2) +
      ggtitle(sprintf("Volcano-Plot: %s (%s = %s)",contrast, 
                      ifelse(FDR_ADJ,"PADJ cutoff","p-value cutoff"),
                      ifelse(FDR_ADJ,P_CUTOFF,P_CUTOFF)))
    
    ggsave(plot = volcano_plot, device="png",
           filename = volcano_Plot_file1, width = 7, height = 7)
       
    #print volcano plot with x axis in linear scale
    
    volcano_Plot_file2 = file.path(plots_dir, paste0("Volcano_Plot_", contrast, "_FC_linScale.png"))
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=FC, y=neglog10)) +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-LINEAR_FC_CUTOFF, LINEAR_FC_CUTOFF),color = "blue", lty = 2) +
      geom_abline(slope = 0, intercept = ifelse(FDR_ADJ,-log10(P_CUTOFF),-log10(P_CUTOFF)),color = "blue", lty = 2) +
      ggtitle(sprintf("Volcano-Plot: %s (%s = %s)",contrast, 
                      ifelse(FDR_ADJ,"PADJ cutoff","p-value cutoff"),
                      ifelse(FDR_ADJ,P_CUTOFF,P_CUTOFF)))
    
    ggsave(plot = volcano_plot, device="png",
           filename = volcano_Plot_file2, width = 7, height = 7)
    
    #print volcano plot with X axis in linear scale with point labels on DE proteins
    
    volcano_Plot_file3 = file.path(plots_dir, paste0("Volcano_Plot_", contrast, "_FC_linScale_w_labels.png"))
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=FC, y=neglog10, label=protein)) +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      geom_text(data = subset(data2plot, signif == T), hjust = 0, nudge_x = 1, size = 2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(-LINEAR_FC_CUTOFF, LINEAR_FC_CUTOFF),color = "blue", lty = 2) +
      geom_abline(slope = 0, intercept = ifelse(FDR_ADJ,-log10(P_CUTOFF),-log10(P_CUTOFF)),color = "blue", lty = 2) +
      ggtitle(sprintf("Volcano-Plot: %s (%s = %s)",contrast, 
                      ifelse(FDR_ADJ,"PADJ cutoff","p-value cutoff"),
                      ifelse(FDR_ADJ,P_CUTOFF,P_CUTOFF)))
    
    ggsave(plot = volcano_plot, device="png",
           filename = volcano_Plot_file3, width = 7, height = 7)
    
    #print volcano plot with X axis in log scale but axis lables in linear scale, with point labels on DE proteins
    
    volcano_Plot_file4 = file.path(plots_dir, paste0("Volcano_Plot_", contrast, "_FC_linScale1_w_labels.png"))
    
    volcano_plot <- 
      ggplot(data2plot, aes(x=ratio, y=neglog10, label=protein)) +
      scale_x_continuous(trans='log2') +
      geom_point(aes(colour = signif), 
                 size=2.5) +
      geom_text(data = subset(data2plot, signif == T), hjust = 0, nudge_x = 0.15, size = 2.5) +
      scale_colour_manual(values = c("TRUE"= "red", "FALSE"= "black")) + 
      xlab(expression(fold~change)) +
      ylab(expression(-log[10]~pvalue)) +
      geom_vline(xintercept = c(1/LINEAR_FC_CUTOFF, LINEAR_FC_CUTOFF),color = "blue", lty = 2) +
      geom_abline(slope = 0, intercept = ifelse(FDR_ADJ,-log10(P_CUTOFF),-log10(P_CUTOFF)),color = "blue", lty = 2) +
      ggtitle(sprintf("Volcano-Plot: %s (%s = %s)",contrast, 
                      ifelse(FDR_ADJ,"PADJ cutoff","p-value cutoff"),
                      ifelse(FDR_ADJ,P_CUTOFF,P_CUTOFF)))
    
    ggsave(plot = volcano_plot, device="png",
           filename = volcano_Plot_file4, width = 7, height = 7)
    
    
  }
}

create_heatmaps = function (mats2plot, DE_df, col_data, graphics_dir, file_extension="") {
  
  file_path_start = file.path(graphics_dir, paste0("Heatmap", file_extension, "_"))
  
  ### All proteins heatmap ###
  title = sprintf("Hierarchical clustering of all proteins (n=%s)", nrow(mats2plot[[1]]))
  file_name = paste0(file_path_start, "all_proteins.png")
  pheatmap_data_All_proteins = plot_expression_heatmap(mats2plot[[1]], DE_df, col_data, file_name, plot_title=title) # without columns clustering
  
  file_name = paste0(file_path_start, "all_proteins_columns.png")
  pheatmap_data_All_proteins = plot_expression_heatmap(mats2plot[[1]], DE_df, col_data, file_name, plot_title=title, cluster_cols = TRUE)  # with columns clustering
  
  ### DE proteins heatmap ###
  title = sprintf("Hierarchical clustering of DE proteins (n=%s)", nrow(mats2plot[[2]]))
  file_name = paste0(file_path_start, "DE_proteins.png")
  pheatmap_data_DE_proteins = plot_expression_heatmap(mats2plot[[2]], DE_df, col_data, file_name, plot_title=title)
  
  ### Proteins from example list heatmap ###
  if (exists("protein_list_from_file")) {
    title = sprintf("Hierarchical clustering of desired proteins (n=%s)", nrow(mats2plot[[3]]))
    file_name = paste0(file_path_start, "desired_proteins.png")
    pheatmap_data_desired_DE_proteins = plot_expression_heatmap(mats2plot[[3]], DE_df,col_data, file_name, plot_title=title)
  }
  
  return (pheatmap_data_DE_proteins)
}

plot_expression_heatmap = function (mat2plot, DE_df, col_data, file_name, color_range="red2blue", row_distance_measure="correlation", plot_title = NA, cluster_cols = FALSE) {
  
  # heatmap row annotation
  row_annotation <- DE_df[, str_detect(names(DE_df),"pass") & !names(DE_df)=="pass_any",drop=F]
  colnames(row_annotation) <- colnames(row_annotation) %>% str_replace("pass.","")
  #rownames(row_annotation) <- row.names(DE_df)  #Vered: commented this out on 22.11.2023. it seems unnecessary
  row_annotation <- row_annotation[row.names(mat2plot),,drop=F]  #this filters the rows to only contain the genes in mat2plot (which may contain either all proteins or only DE proteins)
  row_annotation[row_annotation==""] <- NA
  #remove row_annotation columns (a.k.a. comparisons) with all NA (0 DE genes), otherwise pheatmap throws error (Vered)
  row_annotation = row_annotation %>% select_if (function(x) any(!is.na(x)))
  
  # heatmap column annotation
  col_annotation = col_data %>% dplyr::select(all_of(EFFECTS))
  
  # as_tibble(rownames = "sample") %>%
  # arrange(Tissue,Type) %>%                  # Enter columns to sort by here
  # as.data.frame()
  # rownames(col_annotation) <- col_annotation$sample
  col_annotation$sample <- NULL
  rownames(col_annotation) = colnames(mat2plot)
  
  
  # Z-scoring: scale expression data by row
  mat2plot <- mat2plot %>% t %>% scale %>% t
  
  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot <- mat2plot[complete.cases(mat2plot), ]
  
  #make color scale
  hmcol = switch (color_range,
                  "red2green"  =colorRampPalette(c("red", "black", "green"))(100),
                  "yellow2blue"=colorRampPalette(viridis(n=3) %>% rev)(255),
                  "red2blue"   =colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
                  "red2blue_W" =colorRampPalette(c("red", "white", "blue"))(100)
  )
  
  png(filename = file_name,
      width = 28,
      height = 20,
      units = "cm",
      res=300)
  pheatmap_data <-
    pheatmap(mat2plot,
             clustering_method = "complete",
             clustering_distance_rows = row_distance_measure,
             clustering_distance_cols = row_distance_measure,
             col = rev(hmcol),
             cluster_cols = cluster_cols,
             annotation_row = row_annotation,
             annotation_col = col_annotation,
             main = plot_title,
             show_rownames = F,
             width = 15, height = 10)                       # Comment these two lines for stdout
  dev.off()
  return (pheatmap_data)
}

plot_expression_heatmap1 = function (mat2plot, EFFECTS, stats_df, col_data, file_name,
                                    color_range="red2green", row_distance_measure="correlation", plot_title = NA,
                                    plot_width=28, plot_height=20, row_annot=T, col_annot=T, show_sample_names=T, show_dend=T) {
  
  #this function is from the Neat_RNA-Seq script
  
  # heatmap row annotation
  if (row_annot) {
    row_annotation <- stats_df[, str_detect(names(stats_df),"pass") & !names(stats_df)=="pass_any" & !names(stats_df)=="pass_combined",drop=F]
    colnames(row_annotation) <- colnames(row_annotation) %>% str_replace("pass.","")
    rownames(row_annotation) <- stats_df$gene
    row_annotation <- row_annotation[row.names(mat2plot),,drop=F]
    row_annotation[row_annotation==""] <- NA
    #remove row_annotation columns (a.k.a. comparisons) with all NA (0 DE genes), otherwise pheatmap throws error (Vered)
    row_annotation = row_annotation %>% select_if (function(x) any(!is.na(x)))
    #if nothing was left, don't include row annotation in the plot
    if(ncol(row_annotation) == 0) {
      row_annotation = NA
    }
  } else {
    row_annotation = NA
  }
  
  
  # heatmap column annotation
  if (col_annot) {
    col_annotation = col_data %>% select(all_of(EFFECTS))
    
    # as_tibble(rownames = "sample") %>%
    # arrange(Tissue,Type) %>%                  # Enter columns to sort by here
    # as.data.frame()
    # rownames(col_annotation) <- col_annotation$sample
    col_annotation$sample <- NULL
  } else {
    col_annotation = NA
  }
  
  # show or hide row dendogram
  if (show_dend) {
    dendogram_height = 50
  } else {
    dendogram_height = 0
  }
  
  # Z-scoring: scale expression data by row
  mat2plot <- mat2plot %>% t %>% scale %>% t
  
  #Vered 24.9.2020: remove rows which contain NaN values (otherwise pheatmap throws an error)
  mat2plot <- mat2plot[complete.cases(mat2plot), ]
  
  #make color scale
  hmcol = switch (color_range,
                  "red2green"  =colorRampPalette(c("red", "black", "green"))(100),
                  "red2green.via.gray"  =colorRampPalette(c("brown2", "grey", "seagreen4"))(100),
                  "yellow2blue"=colorRampPalette(viridis(n=3) %>% rev)(255),
                  "red2blue"   =colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(100),
                  "red2blue_W" =colorRampPalette(c("red", "white", "blue"))(100)
  )
  
  png(filename = file_name,
      width  = plot_width,
      height = plot_height,
      units  = "cm",
      res=300)
  pheatmap_data <-
    pheatmap(mat2plot,
             clustering_method = "complete",
             clustering_distance_rows = row_distance_measure,
             col = rev(hmcol),
             cluster_cols = F,
             annotation_row = row_annotation,
             annotation_col = col_annotation,
             main = plot_title,
             show_rownames = F,
             show_colnames = show_sample_names,
             treeheight_row = dendogram_height,
             width = 15, height = 10                       # Comment these two lines for stdout
    )
  dev.off()
  
  return (pheatmap_data)
  
}

GO_barplot = function(enrichGO_object, enricher, OrgDb, ont, Simplify = FALSE, showCategory = 10) {
  
  result_file <- sprintf("%s/Graphics/GO Plots/%s/GO_%s_%s_%s_Enrichment%s_barplot.png",
                         results_dir,
                         ifelse(enricher == "enrichGO",
                                ifelse(Simplify == TRUE, sprintf("%s/Simplify",OrgDb), OrgDb),
                                ifelse(Simplify == TRUE, sprintf("%s/Simplify","GOMap"), "GOMap")),
                         enricher, 
                         ifelse(enricher == "enrichGO", sub("_OrgDb","",OrgDb), "GOMap"),
                         switch(ont, 
                                "BP" = "Biological_Process", 
                                "MF" = "Molecular_Function", 
                                "CC" = "Cellular_Component"),
                         ifelse(Simplify == TRUE,"_Simplify",""))
  
  title <- sprintf("GO %s Pathways",
                   switch(ont, 
                          "BP" = "Biological", 
                          "MF" = "Molecular", 
                          "CC" = "Cellular"))
  
  plot <- barplot(enrichGO_object, 
                  x = 'GeneRatio',
                  showCategory = showCategory, 
                  color = "p.adjust",
                  title = title,
                  font.size = 10.5)
  plot$labels$y = "GeneRatio"
  ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 6)
}

GO_dotplot = function(enrichGO_object, enricher, OrgDb, ont, Simplify = FALSE, showCategory = 10) {
  
  result_file <- sprintf("%s/Graphics/GO Plots/%s/GO_%s_%s_%s_Enrichment%s_dotplot.png",
                         results_dir,
                         ifelse(enricher == "enrichGO",
                                ifelse(Simplify == TRUE, sprintf("%s/Simplify",OrgDb), OrgDb),
                                ifelse(Simplify == TRUE, sprintf("%s/Simplify","GOMap"), "GOMap")),
                         enricher, 
                         ifelse(enricher == "enrichGO", sub("_OrgDb","",OrgDb), "GOMap"),
                         switch(ont, 
                                "BP" = "Biological_Process", 
                                "MF" = "Molecular_Function", 
                                "CC" = "Cellular_Component"),
                         ifelse(Simplify == TRUE,"_Simplify",""))
  
  title <- sprintf("GO %s Pathways",
                   switch(ont, 
                          "BP" = "Biological", 
                          "MF" = "Molecular", 
                          "CC" = "Cellular"))
  
  plot <- dotplot(enrichGO_object, 
                  #x = 'GeneRatio',
                  showCategory = showCategory, 
                  color = "p.adjust",
                  title = title,
                  font.size = 10.5)
  #plot$labels$y = "GeneRatio"
  ggsave(result_file, plot, device="png", scale = 1, width = 9, height = 6)
  
}

##### Utilities #####

create_dir = function (directory) {
  #create a new directory if it does not yet exist (path is relative to the current script)
  if(!dir.exists(directory)) dir.create(directory)
  return(directory)
}

remove_LFQ_from_all_colnames_if_exist = function (x) {
  #if all col names have leading LFQ.LOG2, remove it
  if (all(str_detect(colnames(x), "^LFQ.LOG2"))) {
    colnames(x) = gsub("LFQ.LOG2.", "", colnames(x))
  }
  return (x)
}

remove_intensity_from_all_colnames_if_exist = function (x) {
  #if all col names have leading LFQ.intensity, remove it
  if (all(str_detect(colnames(x), "^LFQ.intensity"))) {
    colnames(x) = gsub("LFQ.intensity.", "", colnames(x))
  }
  return (x)
}

get_conditions_from_samples = function(samples) {
  conditions <- vector()
  for (sample in unique(samples)) {
    condition <- col_data %>%
      filter(SampleID == sample) %>%
      select(Condition) %>%
      unlist %>%
      as.character
    conditions <- c(conditions,rep(condition,sum(str_count(samples,sample))))
  }
  return (conditions)
}

export_imputed_data = function (df.FI, file_name) {
  df.output = select(df.FI, -starts_with("impute"))
  write.table(df.output, file =file_name, quote = F, sep = "\t",
              eol = "\n", na = "", dec = ".", row.names = T,
              col.names = TRUE)
}
import_imputed_data = function(df, Log2.names, import_directory) {
  mult_df.FI = list()
  
  col.names = unlist(Log2.names)
  impute.names = sub("^LFQ.LOG2", "impute", col.names)  #create column names "impute.C1", "impute.C2"... (columns data will be TRUE/FALSE flags)
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] = lapply(col.names, function(x) !is.finite(df[, x]))
  
  for (n in 1:length(dir(path = file.path(import_directory,"Datasets","Imputed repetitions"), pattern = ".txt"))) {
    LFQ.imputs <- read.delim(file = file.path(import_directory,
                                              "Datasets",
                                              "Imputed repetitions",
                                              sprintf("Imputed.%s.txt",n)),
                             header = TRUE, stringsAsFactors = FALSE)
    df <- dplyr::mutate(df, dplyr::select(LFQ.imputs,all_of(col.names)))
    mult_df.FI[[n]] = df
  }
  return (mult_df.FI)
}

prepare_data2write = function(final_results) {
  results_to_excel = final_results
  
  # convert manual_cutoffs columns to 'formula' class
  formula_columns <- results_to_excel %>% names %>% str_detect(pattern = "manual_cutoffs") %>% which
  for (col in formula_columns)
    class(results_to_excel[,col]) <- c(class(results_to_excel[,col]), "formula")
  
  return (results_to_excel)
}

prepare_data2write_DE = function(results_to_excel, pheatmap_data_DE_proteins) {
  results_to_excel_DE <- results_to_excel[(!is.na(results_to_excel$pass_any_contrast) & 
                                             results_to_excel$pass_any_contrast == "1"),]
  diff_proteins_order <- cbind(Protein=DE_proteins[pheatmap_data_DE_proteins$tree_row$order],
                               order = seq(from=1,by=1,along.with = DE_proteins)) %>% 
    as_tibble()
  
  results_to_excel_DE <- results_to_excel_DE %>% 
    as_tibble %>% 
    left_join(y = diff_proteins_order,by = c("Protein"))
  
  # Adding Z-score data:
  zscore <- mats2plot[[2]] %>% t %>% scale %>% t
  colnames(zscore) <- zscore %>% colnames %>% paste0(.,".zscore")
  zscore <- zscore %>% as_tibble(rownames = "Protein") 
  
  results_to_excel_DE <- results_to_excel_DE %>% 
    left_join(zscore,by="Protein")
  
  dim(results_to_excel_DE)
  results_to_excel_DE <- results_to_excel_DE %>% unique()
  dim(results_to_excel_DE)
  
  return (results_to_excel_DE)
}

output_main_excel_table = function(results_to_excel, grouping_headers, 
                                   contrast_headers, FDR_ADJ, all_results_with_DE) {
  # Create the excel workbook
  excel_wb = createWorkbook()
  
  metadata_cols = EFFECTS
  metadata_rows <- length(metadata_cols)+1
  
  # index in grouping_header for which to produce column info above the table (2=counts data)
  grouping_header_meta <- c(3)
  
  # Add sheet: Cutoffs
  addWorksheet(excel_wb, sheetName = "Cutoffs", gridLines = TRUE)
  
  # Write data: Cutoffs
  writeData(excel_wb, "Cutoffs", 
            #x = c(ifelse(FDR_ADJ,"Adjusted pvalue (FDR)","p-value"),"linear Fold Change (linearFC)"),
            x = c("p-value","Adjusted pvalue (FDR)","linear Fold Change (linearFC)"), 
            startCol = 2, 
            startRow = 4)
  writeData(excel_wb, "Cutoffs", 
            #x = c(ifelse(FDR_ADJ,PADJ_CUTOFF,P_CUTOFF),LINEAR_FC_CUTOFF),
            x = c(ifelse(FDR_ADJ,"",P_CUTOFF),
                  ifelse(FDR_ADJ,P_CUTOFF,""),
                  LINEAR_FC_CUTOFF), 
            startCol = 3, 
            startRow = 4)
  
  createNamedRegion(excel_wb, "Cutoffs", cols=3, rows=4, "PVAL_CO")
  createNamedRegion(excel_wb, "Cutoffs", cols=3, rows=5, "FDR_CO")
  createNamedRegion(excel_wb, "Cutoffs", cols=3, rows=6, "LFC_CO")
  
  style_COs_Used <- 
    createStyle(border = "TopBottomLeftRight",
                borderStyle = "thick", 
                fgFill = "green",
                halign = "center")
  style_COs_NotUsed <- 
    createStyle(border = "TopBottomLeftRight",
                borderStyle = "thick", 
                fgFill = "red",
                halign = "center")
  if (FDR_ADJ) {
    addStyle(excel_wb, "Cutoffs", style_COs_NotUsed, rows=4, cols=3, gridExpand = FALSE, stack = FALSE)
    addStyle(excel_wb, "Cutoffs", style_COs_Used, rows=5:6, cols=3, gridExpand = FALSE, stack = FALSE)
  } else {
    addStyle(excel_wb, "Cutoffs", style_COs_NotUsed, rows=5, cols=3, gridExpand = FALSE, stack = FALSE)
    addStyle(excel_wb, "Cutoffs", style_COs_Used, rows=c(4,6), cols=3, gridExpand = FALSE, stack = FALSE)
  }
  setColWidths(excel_wb, "Cutoffs", cols = 2, widths = "auto")
  
  
  # Add sheet: All proteins
  sheetName = "All proteins"
  addWorksheet(excel_wb, sheetName = sheetName, gridLines = TRUE)
  
  # Write data: All proteins
  writeDataTable(excel_wb, 
                 sheet = sheetName, 
                 x = results_to_excel, # add this if for testing: %>% head(n=2150),
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  addStyle(excel_wb, 
           sheet = sheetName, 
           style =  createStyle(numFmt = "0.000"), 
           rows = (metadata_rows+2):(metadata_rows+1+dim(results_to_excel)[1]), 
           cols = (grouping_headers[1]+1):sum(grouping_headers[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  # Adding additional column information above column headers
  writeData(excel_wb, 
            sheet = sheetName, 
            x = contrast_headers %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  # Adding col_data where required:
  for (i in cumsum(grouping_headers)[grouping_header_meta-1]){
    writeData(excel_wb,
              sheet = sheetName,
              x = col_data[,metadata_cols] %>% as.matrix %>% t,
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
  }
  
  mycolors = rep(brewer.pal(12, name = "Paired"),3)
  colorscheme <- data.frame(colnum = grouping_headers,
                            fgcol= mycolors[1:length(grouping_headers)],
                            stringsAsFactors = F)
  
  colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                yes = "dimgray",
                                no = "white")
  
  for (i in 1:length(grouping_headers)) {
    cat(sprintf("%s of %s\n",i,length(grouping_headers)))
    # Add title color:
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(excel_wb, 
             sheet = sheetName, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1,cumsum(grouping_headers)[i-1]+1):cumsum(grouping_headers)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator:
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(excel_wb, 
             sheet = sheetName, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(results_to_excel)[1]), 
             cols = cumsum(grouping_headers)[i], 
             gridExpand = FALSE, stack = TRUE)
    
  }
  if (isFALSE(all_results_with_DE)) {
    # openXL(wb)
    file <- sprintf("%s_%s.xlsx",sub(".xlsx","",Final_excel_file),
                    ifelse(FDR_ADJ==TRUE,sprintf("PADJ_%s",P_CUTOFF),
                           sprintf("P_%s",P_CUTOFF)))
    # Writing workbook:
    print(sprintf("Writing file: %s", file))
    saveWorkbook(excel_wb, 
                 file = file, 
                 overwrite = TRUE) ## save to working directory
  }
  return (excel_wb)
}

output_DE_excel_table = function(results_to_excel_DE, grouping_DE_headers, 
                                 contrast_DE_headers, FDR_ADJ, all_results_with_DE) {
  sheetName = "DE Proteins"
  metadata_cols = EFFECTS
  metadata_rows <- length(metadata_cols)+1
  grouping_header_meta_DE <- c(2,length(grouping_DE_headers))
  
  # Create the excel workbook
  excel_wb_DE = createWorkbook()
  
  # Add sheet: DE Proteins
  addWorksheet(excel_wb_DE, sheetName = sheetName, gridLines = TRUE)
  
  # Add data: DE Proteins
  writeDataTable(excel_wb_DE, 
                 sheet = sheetName, 
                 x = results_to_excel_DE, # add this if for testing: %>% head(n=2150),
                 startRow = metadata_rows+1,
                 colNames = TRUE,
                 rowNames = FALSE,
                 tableStyle = "TableStyleLight9")
  
  addStyle(excel_wb_DE, 
           sheet = sheetName, 
           style =  createStyle(numFmt = "0.000"), 
           rows = (metadata_rows+2):(metadata_rows+1+dim(results_to_excel_DE)[1]), 
           cols = (grouping_DE_headers[1]+1):sum(grouping_DE_headers[1:2]), 
           gridExpand = TRUE, stack = FALSE)
  
  writeData(excel_wb_DE, 
            sheet = sheetName, 
            x = contrast_DE_headers %>% as.matrix %>% t,
            startRow = metadata_rows,
            colNames = F,
            rowNames = F)
  
  # Adding col_data where required:
  for (i in cumsum(grouping_DE_headers)[grouping_header_meta_DE-1]){
    writeData(excel_wb_DE, 
              sheet = sheetName, 
              x = col_data[,metadata_cols] %>% as.matrix %>% t , 
              startRow = 1,
              startCol = i,
              colNames = FALSE,
              rowNames = TRUE)
    
    mycolors = rep(brewer.pal(12, name = "Paired"),3)
    colorscheme <- data.frame(colnum = grouping_DE_headers,
                              fgcol= mycolors[1:length(grouping_DE_headers)],
                              stringsAsFactors = F)
    
    colorscheme$textcol <- ifelse((as(hex2RGB(colorscheme$fgcol),"polarLUV"))@coords[,1] > 65, 
                                  yes = "dimgray",
                                  no = "white")
  }
  
  for (i in 1:length(grouping_DE_headers)) {
    cat(sprintf("%s of %s\n",i,length(grouping_DE_headers)))
    # Add title color:
    hs1 <- createStyle(fgFill = colorscheme$fgcol[i], 
                       fontColour = colorscheme$textcol[i])
    addStyle(excel_wb_DE, 
             sheet = sheetName, 
             style =  hs1, 
             rows = metadata_rows+1, 
             cols = ifelse(i==1,1:grouping_DE_headers[1],
                           cumsum(grouping_DE_headers)[i-1]+1):cumsum(grouping_DE_headers)[i], 
             gridExpand = FALSE, stack = TRUE)
    # Add group separator:
    hs2 <- createStyle(border = "Right",
                       borderColour = "black",
                       borderStyle = "thick")
    addStyle(excel_wb_DE, 
             sheet = sheetName, 
             style = hs2, 
             rows = (metadata_rows+2):(metadata_rows+1+dim(results_to_excel_DE)[1]), 
             cols = cumsum(grouping_DE_headers)[i], 
             gridExpand = FALSE, stack = TRUE)
  }
  if (isFALSE(all_results_with_DE)) {
    # openXL(wb)
    file <- sprintf("%s_%s.xlsx",sub(".xlsx","",Final_excel_DE_file),
                    ifelse(FDR_ADJ==TRUE,sprintf("PADJ_%s",P_CUTOFF),
                           sprintf("P_%s",P_CUTOFF)))
    # Writing workbook:
    print(sprintf("Writing file: %s", file))
    saveWorkbook(excel_wb_DE, 
                 file = file, 
                 overwrite = TRUE) ## save to working directory
  }
  return (excel_wb_DE)
}

output_AH_query = function(query, AH_specie) {
  if (length(query) == 0) {
    print(sprintf("query for '%s' found 0 records", AH_specie))
    return (NULL)
  }
  if (length(query) > 1) {
    print(sprintf("query for '%s' found %s records", AH_specie, length(query)))
    for(i in 1:length(query))
      print(sprintf("Record %s: %s (specie: %s)", i, query$description[i], query$species[i]))
    print("Please update 'AH_Specie' parameter to be more specific")
    return (NULL)
  }
  if (length(query) == 1) {
    print(sprintf("query for '%s' found 1 record", AH_specie))
    print(sprintf("Downloading record '%s' for specie: %s",query$title, query$species))
    print(sprintf("Record description: %s", query$description))
    AH_OrgDb <- query[[query$ah_id]]
    return (AH_OrgDb)
  }
}

z_score = function (expr_matrix) {
  expr_matrix %>% t %>% scale %>% t
}

export_table = function (data, file_name, na_symbol = "NaN") {
  data <- cbind(Name = rownames(data), data)
  write.table(x = data,
              file = file_name,
              quote = F,
              sep = "\t",
              na = na_symbol,
              row.names = F,
              col.names = T)
}