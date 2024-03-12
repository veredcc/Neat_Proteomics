# Global settings for Neat_proteomics script

# Update this file with the parameters and input files of your current project, before executing "Neat_proteomics_1.0.Rmd".

##### Global parameters #####

##### General ######

project_name = "Proj_A"

analysis_round = "Analysis_01"

EXPORT_INTERMEDIATE_FILES = TRUE #(T or F)

##### Filtering #####

#define data filtering conditions
#min_count is the required minimum number of samples having non-zero intensities per condition
#it needs to be met in at least one of the conditions

min_count     = c(2,2) 

##### Experiment effects #####

#define experiment effects

EFFECTS = c("Batch", "Condition") ###vered: need to write a solution for cases that there is only one effect of interest (instead of 2)

##### Imputation #####

#define cutoffs for imputation

WIDTH     = 0.2 #Perseus and our function default is 0.3
DOWNSHIFT = 1.6 #Perseus and our function default is 1.8

#define imputation parameters

NO_REPETITIONS = 10 #number of imputed datasets to create
MIN_NO_PASSED  = 8  #minimum no. of imputed datasets which have to pass the DE cutoff (per contrast)

##### DE analysis #####

#define cutoffs for DE analysis

FDR_ADJ = FALSE #TRUE means p-value cutoff is after FDR. FALSE means use unadjusted p-value
P_CUTOFF = 0.05
LINEAR_FC_CUTOFF = 1.5

DESIGN_TYPE = 'batch' #options: 'simple', 'batch', 'batch-multi_level'
BATCH_FACTOR = "Batch"

#note: to include interaction, add corresponding contrasts to the contrasts file,
#as explained in Limma User Guide, section "Interaction Models: 2 X 2 Factorial Designs", subsection "Analysing as for a Single Factor"
#for example, see the last contrast in the contrasts.txt file here:
#Contrast_name	Factor	Numerator	Denominator
#WT:LPS-Cont	Condition	WT_LPS	WT_Cont
#KO:LPS-Cont	Condition	KO_LPS	KO_Cont
#Cont:KO-WT	Condition	KO_Cont	WT_Cont
#LPS:KO-WT	Condition	KO_LPS	WT_LPS
#Interaction	Condition	(KO_LPS-KO_Cont)	(WT_LPS-WT_Cont)




#define cutoffs for present-absent analysis
#"count" refers to number of non-missing values per condition (per protein)

MIN_COUNT_FOR_PRESENT = 9
MAX_COUNT_FOR_ABSENT = 2

##### Visualization #####

# APPLY_BATCH_CORR = T  #modify expression values for display purposes, e.g. in hierarchical clustering

#batch correction  ##new code 20.12.2023

#note, in the DE analysis section we included a variable called BATCH_FACTOR. That parameter is used for the experiment design in the Limma DE analysis.
#here we define BATCH_EFFECTS - the batch effect(s) for batch correction. We do not use the same variable as above because in limma it is possible to correct for two batches, and we want to allow that.

REMOVE_BATCH_EFFECT = T     #TRUE or FALSE
BATCH_EFFECTS = c('Batch')  #indicate colunm name(s) from the Experiment_design file. sva can use only one. limma can use up to two.
BATCH_CORR_METHOD = 'sva'   #sva or limma

##### Clustering #####

##new code 20.12.2023

#set group factor for the binary pattern calculations and for partition clustering
#(this is usually the factor that was used for the contrasts)

GROUP = "Condition"

#set correlation cutoff for binary pattern calculation

CORR_CUTOFF = 0.8

#set min counts cutoff for 1s in the binary pattern calculation

COUNTS_CUTOFF = 0 #in proteomics use 0 (requiring that all samples which are 1 in the pattern will be non-zero, in RNA-Seq use no. of min counts in unnorm data, e.g. 500)

#show only top n DE genes in hierarchical clustering?

SHOW_TOP_GENES_IN_HEATMAP = FALSE   #TRUE: show only top DE genes in heatmap. FALSE: show all DE genes in heatmap
NR_TOP_GENES = 5000  #if SHOW_TOP_GENES_IN_HEATMAP is set to FALSE, this parameter will be neglected

#partition clustering

K_FIXED = NA #no. of requested clusters. if NA, K_MAX will be used
K_MAX = 20

GROUP1 = "Condition"   #usually GROUP will be the main treatment (e.g. diet), and GROUP1 will be another biological factor, e.g. age

#manual clustering
PERFORM_MANUAL_CLUSTERING = T  #TRUE or FALSE
COLUMNS_FOR_MAN_CLUSTERING = c(1:3)  #currently manual clustering can only work on three biological groups. specify column indexes on mat2plot
#for manual clustering option 1
MAN_CLUST_FC_CUTOFF = 1.2 #orig was 1.5
#for manual clustering option 2
MAN_CLUST_DESEQ_FC_CUTOFF = 1.3
MAN_CLUST_DESEQ_P_CUTOFF = 0.05   #wo FDR

##### Enrichment analysis #####

##old code from Gil

# RUN_GO_ENRICHMENT = T

# define GO related parameters 
# GO_CUTOFF = 0.05
# GO_MIN_GENES = 4
# Simplify_CUTOFF = 0.7 # The recommended cutoff for simplify analysis is 0.7
# USE_BioC_OrgDb = FALSE # TRUE means the organism is listed in "http://bioconductor.org/packages/release/BiocViews.html#___OrgDb"
# Organism = "org.Hs.eg.db" # Organism will be ignored if USE_BioC_OrgDb is FALSE
# USE_AH_OrgDb = FALSE # TRUE means to search for OrgDb record in AnnotationHub
# AH_Specie = "SPECIE" # Specie search term in AnnotationHub

##new code 20.12.2023

FILTER_KEGG_PATHWAYS_BY_TAXON = NA  #provide KEGG taxon (name or ID?) to filter the pathways. for no filtering use NA

#passed to clusterProfiler::enrichr
ENRICHMENT_PVAL_CUTOFF = 0.05
ENRICHMENT_PADJ_METHOD = 'fdr' #options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#passed to enrichplot::dotplot as argument showCategory
MAX_TERMS_IN_DOTPLOT = 30  #number of enriched terms to display in dotplot

##### Output to Excel#####

#define excel related parameters

all_results_with_DE = FALSE # TRUE means to create 1 workbook for all data (All proteins & DE proteins)
INCL_PA_IN_EXCEL    = FALSE # TRUE means to include present-absent analysis in the final results file


##### Input files and directories #####

proteinGroups_file        = "proteinGroups.txt"

samples_file              = "Samples_wo_CTL.1.txt"

experiment_design_file    = "Experiment_design_wo_CTL.1.txt"  # Assumes first column is 'SampleID'. The experimental groups should be in column "Condition"
                                                              # Currently the script assumes only one treatment ("Condition"), and optional batch(es)

contrasts_file            = "Contrasts.txt"

proteins_list_file        = "protein_list_example.txt"  #external list of proteins of interest

#directory with data required for enrichment analyses
#kegg data are taken from Vered Perl script.
#GO data and the Annotation.tab file are from the R script.
#The perl script is located on Vered’s PC, at
#Programming/Perl_workspace/Assaf_Rudich_03_HFD_KEGG_from_Vered_script/KEGG_2_ensembl_conversion8.pl

functional_annot_dir     = "Func_annot_data"  