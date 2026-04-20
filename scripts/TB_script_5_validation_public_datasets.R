# This script contains validation of the 4-gene signature in public datasets
# Validating 4-gene signature in publicly available datasets
# We validate in TB datasets and also datasets including multiple lung diseases
# We calculate mean of z-scored expression of the gene signature (plot = boxplots of signature scores) and also use signature scores to create glm model where score predicts disease, to calculate roc curves (plot = roc)
# As of 09/03, no longer using GSVA method, using mean of z-scored expression instead

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

# Mac
my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("/Volumes/One Touch/RBMB/rlibrary")


library("readxl")
library("GEOquery")
library("edgeR")
library("GSVA")
library("pROC")
library("R.utils")
library("ggplot2")
library("rstatix")
library("ggpubr")
library("DESeq2")
library("limma")
library("stringr")
library("tidyverse")

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")
if(!exists(output.dir)) dir.create(output.dir)

gc()


# ================================================================================== #
# 1. TB DATASETS ===================================================================
# ================================================================================== #

# GSE89403  -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE89403"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))


this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

#Ran in HPC
# # load counts table from GEO
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE89403", "file=GSE89403_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
# 
# # load gene annotations
# apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
# annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
# rownames(annot) <- annot$GeneID

# all(row.names(tbl) == annot$GeneID) #39,376 genes and 914 columns


## 1) Load public data  ----------------------------------------------------
# Ran in R
raw_counts <- read.delim("raw_counts.txt", sep = " ")
gene_annot <- read.table("gene_annot.txt")
raw_metadata <- read.table("raw_metadata.txt")


raw_clinical <- raw_metadata[,c("sample_code.ch1",
                                "disease.state.ch1",
                                "characteristics_ch1.2",
                                "characteristics_ch1.4",
                                "characteristics_ch1.5",
                                "characteristics_ch1.6",
                                "treatmentresult.ch1")]

colnames(raw_clinical) <- c("sample_id",
                            "disease",
                            "subject",
                            "treatmentresult",
                            "time",
                            "timetonegativity",
                            "treatmentresult")

raw_clinical[3:6] <- lapply(raw_clinical[,3:6], function(x) gsub(".*: ", "", x))


all(row.names(raw_clinical) == colnames(raw_counts))

clinical <- raw_clinical


## 2) Get groups to be compared --------------------------------------------

clinical[which(clinical$disease == "Healthy Controls"),"disease"] <- "Healthy" #community control?
clinical[which(clinical$disease == "Lung Dx Controls"),"disease"] <- "Lungdx_ctrl" #lung disease
clinical[which(clinical$disease == "MTP Controls"),"disease"] <- "MTP_ctrl" #MTP
clinical[which(clinical$disease == "TB Subjects"),"disease"] <- "TB"
if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}

clinical$group <- paste0(clinical$disease,"_",clinical$time)
table(clinical$group)


#1) For controls, disregard week ??? but keep the controls seperate
clinical[which(clinical$disease == "Healthy"), "group"] <- "Healthy"
clinical[which(clinical$disease == "Lungdx_ctrl"), "group"] <- "Lungdx_ctrl"
clinical[which(clinical$disease == "MTP_ctrl"), "group"] <- "MTP_ctrl"

counts_vst <- vst(as.matrix(raw_counts))
counts_vst <- counts_vst[,row.names(clinical)]

write.table(counts_vst, file.path("counts_vst.txt"))
write.table(clinical, file.path("clinical.txt"))

counts_norm <- counts_vst

gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }


## 3) Mean of z-scored expression ---------------------------
mean_zscore_func <- function(){

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- counts_norm[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

if(all(row.names(mean_sig_zscore) != row.names(clinical))){ stop("Row names do not match", this.accession.no)}
return(mean_sig_zscore)


} #close mean_sig_zscore function

#Run the function and save the score
mean_sig_zscore <- mean_zscore_func()


## 3.1) Boxplot ---------------------------
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                        disease = clinical$disease,
                                    group = clinical$group,
                                    treatment_outcome = clinical$treatmentresult))



boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 




# Specify order that we want the x-axis variables in for the plot
x_order <- c("Healthy", "Lungdx_ctrl", "MTP_ctrl", "TB_DX", "TB_day_7", "TB_week_4", "TB_week_24")

# Make boxplot function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot_func <- function(outcome){

#get all comparisons
my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)

boxplot$group <- factor(boxplot$group, levels = x_order)
boxplot$score <- as.numeric(boxplot$score)

# calculate p values
stat.table<- boxplot  %>%
  wilcox_test(score ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]
lowest_bracket <- max(boxplot$score) + 0.05*(max(boxplot$score))
stat.table$y.position <- seq(lowest_bracket, by= 0.2, length.out = nrow(stat.table))


ggplot(boxplot, aes(
  x = factor(group, level = x_order),
  y = score,
  group = group)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  
  stat_pvalue_manual(stat.table,
                     label = "p",
                     tip.length = 0.01,
                     vjust = 0.3,
                     size = 4)+
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  theme(axis.text.x = element_text(size = 15))+
  ylab (label = "Signature Score") +
  xlab (label = "Disease")

} #end boxplot function



#Run boxplot function for all, cured and non-cured
for (outcome in c("all_outcomes", "not_cured", "cured")){
  
    if(outcome == "all_outcomes"){
  boxplot <- boxplot_all}
  
  if(outcome == "not_cured"){
  boxplot <- boxplot_all[which(boxplot_all$treatment_outcome == "Not Cured" | boxplot_all$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    boxplot <- boxplot_all[which(boxplot_all$treatment_outcome == "Definite Cure" | boxplot_all$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
    boxplotfig <- boxplot_func(outcome = outcome)
    boxplotfig <- boxplotfig +   
      labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Raw counts were VST normalised (DESeq2 1.42.1)\n",
       "P values from Mann-Whitney U test shown")) 
    
    ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, "_", outcome, ".png")),
       width = 3500,
       height = 3600,
       units = "px" )

}


## 4) Validation  ------------------------------------------------------
# Options for ROC
## a) multiclass.roc(): This is a one-vs-rest approach (macro-average). Where one group i selected as the refernece/positive class and the rest are grouped into a negative class. It takes turns with each group being the reference class, giving us AUCs for each reference vs rest model. Then finally it calculates the average of the AUCs, a single value to classify overall performance
# probbly not great because if our one vs rest model is healthy vs rest, that means rest includes TB day 7, week 4 and week 24. TB week 24 would theoretically have expression profiles that are closer to healthy because patients have been treated. so the model would be confused
## b) Pairwise ROC: run one vs one (binary) ROC fr each pair of timepoints - seperate binary logistic regression models for each comparison

#Decision : Use pairwise ROC! multiclass won't work well in this scenario (best if the variables aren't longitudinal / are unrelated)


# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Healthy", "TB_DX"),
  c("MTP_ctrl", "TB_DX"),
  c("TB_DX", "TB_day_7"),
  c("TB_DX", "TB_week_4"),
  c("TB_DX", "TB_week_24"),
  c("Healthy", "TB_week_24"),
  c("MTP_ctrl", "TB_week_24")
)

# Define for plotting
comparison_levels <- c(
  "HealthyvsTB_DX",
  "MTP_ctrlvsTB_DX",
  "TB_DXvsTB_day_7",
  "TB_DXvsTB_week_4",
  "TB_DXvsTB_week_24",
  "HealthyvsTB_week_24",
  "MTP_ctrlvsTB_week_24"
)


comparison_plotlabel_levels <- c(
  "Healthy vs TB_T0",
  "MTP Controls vs TB_T0",
  "TB_T0 vs TB_Day7",
  "TB_T0 vs TB_Wk4",
  "TB_T0 vs TB_Wk24",
  "Healthy vs TB_Wk24",
  "MTP Controls vs TB_Wk24"
)


#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

roc_func <- function(outcome){

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical_treat[clinical_treat$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # Mean-zscore

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)  #default 95% CI is computed with 2000 stratified bootstrap replicates.
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  
  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1,"vs",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
} # close pair
  


# 5) ROC Curves & Forest plot-----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison, levels = comparison_levels)

roc_data$Comparison_plotlabel <- roc_data$Comparison
levels(roc_data$Comparison_plotlabel) <-  comparison_plotlabel_levels

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")

return(list(
  res_table = res_table,
  forestplot_res_table = forestplot_res_table,
  roc_data = roc_data
))

} #close roc func


#Run function
for (outcome in c("all_outcomes", "not_cured","cured")){
    if(outcome == "all_outcomes"){
  clinical_treat <- clinical
  }
  
  if(outcome == "not_cured"){
  clinical_treat <- clinical[which(clinical$treatmentresult == "Not Cured" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    clinical_treat <- clinical[which(clinical$treatmentresult == "Definite Cure" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
    }
    
      roc_func_res <- roc_func(outcome = outcome)
res_table <- roc_func_res$res_table
forestplot_res_table <- roc_func_res$forestplot_res_table

write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))

}



disease_roc_subset <- c(
"Healthy vs TB_T0",
"MTP Controls vs TB_T0",
"Healthy vs TB_Wk24",
"MTP Controls vs TB_Wk24"
)


timepoint_roc_subset <- c(
"TB_T0 vs TB_Day7",
"TB_T0 vs TB_Wk4",
"TB_T0 vs TB_Wk24"
  )

roc_plot_width = 3000
roc_plot_height = 3200

disease_roc_plot_func <- function(legend_nrow = 2){
  
# Disease plot
roc_data <- roc_func_res$roc_data
disease_roc_data <- roc_data[which(roc_data$Comparison_plotlabel %in% disease_roc_subset),]

ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = paste0("ROC - Control vs TB", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature:  TAP1, GBP5, GBP2, FCGR1CP") 
}

# disease_roc <- disease_roc_plot_func()
# ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
#        width = roc_plot_width, 
#        height = roc_plot_height, 
#        units = "px")

#Timepoint plot

timepoint_roc_plot_func <- function(legend_nrow = 2){

roc_data <- roc_func_res$roc_data
timepoint_roc_data <- roc_data[which(roc_data$Comparison_plotlabel %in% timepoint_roc_subset),]

ggplot(timepoint_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = paste0("ROC - TB treatment timepoints", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR(Sensitivity)",
    color = "Comparison",
    caption = "Signature:TAP1, GBP5, GBP2, FCGR1CP")
}


# timepoint_roc <- timepoint_roc_plot_func()
# ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome, ".png")), 
#        width = roc_plot_width, 
#        height = roc_plot_height, 
#        units = "px")


## Forest plot -------
forestplot_func <- function(){
forestplot_res_table <- roc_func_res$forestplot_res_table
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = comparison_levels)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 
auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  

  
  return(panel_forest)
} #close function
# 
# panel_forest <- forestplot_func()
#   
# ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
#          width = 15, height = 20, units = "cm",   bg = "white"  )
# 


for (outcome in c("all_outcomes", "not_cured","cured")){
    if(outcome == "all_outcomes"){
  clinical_treat <- clinical
  }
  
  if(outcome == "not_cured"){
  clinical_treat <- clinical[which(clinical$treatmentresult == "Not Cured" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    clinical_treat <- clinical[which(clinical$treatmentresult == "Definite Cure" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
    }
    
  #run roc_func again to get the results for this outcome type
roc_func_res <- roc_func(outcome = outcome)

# disease plot
disease_roc <- disease_roc_plot_func()
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

#timepoint plot
timepoint_roc <- timepoint_roc_plot_func()
ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome, ".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

#forestplot
panel_forest <- forestplot_func()
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1CP"), 
                       size = 12, hjust = 0, x = 0)
  )
ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

}



#Run function
for (outcome in c("all_outcomes", "not_cured","cured")){
    if(outcome == "all_outcomes"){
  clinical_treat <- clinical
  }
  
  if(outcome == "not_cured"){
  clinical_treat <- clinical[which(clinical$treatmentresult == "Not Cured" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    clinical_treat <- clinical[which(clinical$treatmentresult == "Definite Cure" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
    }
    
      roc_func_res <- roc_func(outcome = outcome)
res_table <- roc_func_res$res_table
forestplot_res_table <- roc_func_res$forestplot_res_table

write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))

}


# GSE193777 -----------------------------------------------------------------------------------------------------------------------------------
#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c("gene_annot",lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE193777"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)


## 1) Load public data -------------------------------------------------------------
#Ran in HPC
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE193777", "file=GSE193777_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

gse=getGEO(filename="GSE193777_series_matrix.txt.gz", getGPL = FALSE)

raw_counts <- read.delim("raw_counts.txt", sep = " ")
raw_metadata <- gse@phenoData@data

colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "distinct.stages.of.tb.ch1",
                                "gender.ch1",
                                "mtb.culture.ch1", #collecting blood and culturing the Mycobacterium tuberculosis (MTB) bacteria
                                "qft.result.ch1")] #blood test used to diagnose tuberculosis (TB) by measuring the amount of interferon-gamma (IFN-γ) released by a person's immune cells when exposed to specific TB antigens

colnames(raw_clinical) <- c("sample_id",
                            "disease",
                            "sex",
                            "disease_mtb",
                            "disease_qft")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))

## 2) Get groups to be compared --------------------------------------------
table(clinical$disease)
#NOTE: all patients are household contacts, none are community controls
clinical[which(clinical$disease == "HH controls_baseline"),"group"] <- "Healthy" 
clinical[which(clinical$disease == "Active TB_baseline"),"group"] <- "Active_TB" 
clinical[which(clinical$disease == "Incipient TB/Non-Progressor_baseline"),"group"] <- "Icp_TB_bl" 
clinical[which(clinical$disease == "Incipient TB/Non-Progressor_followup"),"group"] <- "Icp_TB_fu"
clinical[which(clinical$disease == "Sub-clinical TB/Progressor_baseline"),"group"] <- "Sub_TB_bl"
clinical[which(clinical$disease == "Sub-clinical TB/Progressor_followup"),"group"] <- "Sub_TB_fu"

if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}

counts_vst <- vst(as.matrix(raw_counts))
counts_vst <- counts_vst[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

write.table(counts_vst, file.path("counts_vst.txt"))
write.table(clinical, file.path("clinical.txt"))

counts_norm <- counts_vst




## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 

## 3.1) Boxplot ---------------------------
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text.y = element_text(size = 24),
                    axis.text.x = element_text(size = 20),
                    title = element_text(size = 20),
                    legend.position = "None") 


# Specify order that we want the x-axis variables in for the plot
x_order <- c("Healthy", "Icp_TB_bl", "Icp_TB_fu", "Sub_TB_bl", "Sub_TB_fu", "Active_TB")


# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot <- boxplot_all

outcome = "TB" #no outcome data available so just use 'TB' as placeholder
boxplotfig <- boxplot_func(outcome = outcome) 

boxplotfig <- boxplotfig + scale_x_discrete(labels= c("Healthy" = "Healthy",
                             "Active_TB" = "Active\n TB", 
                             "Icp_TB_bl" = "Incipient TB \n baseline",
                             "Icp_TB_fu" = "Incipient TB \n followup",
                             "Sub_TB_bl" = "Subclinical TB \n baseline",
                             "Sub_TB_fu" = "Subclinical TB \n followup")) +
       labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Raw counts were VST normalised (DESeq2 1.42.1)\n",
       "P values from Mann-Whitney U test shown")) 
  


ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
       units = "px" )

 
## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Healthy", "Active_TB"),
  c("Healthy", "Icp_TB_bl"),
  c("Healthy", "Icp_TB_fu"),
  c("Healthy", "Sub_TB_bl"),
  c("Healthy", "Sub_TB_fu"),
  c("Active_TB", "Icp_TB_bl"),
  c("Active_TB", "Icp_TB_fu"),
  c("Active_TB", "Sub_TB_bl"),
  c("Active_TB", "Sub_TB_fu")
)

comparison_levels <- c(
  "Healthy vs Active_TB",
  "Healthy vs Icp_TB_bl",
  "Healthy vs Icp_TB_fu",
  "Healthy vs Sub_TB_bl",
  "Healthy vs Sub_TB_fu",
  "Active_TB vs Icp_TB_bl",
  "Active_TB vs Icp_TB_fu",
  "Active_TB vs Sub_TB_bl",
  "Active_TB vs Sub_TB_fu"
)

comparison_plotlabel_levels <- c(
  "Healthy vs Active TB",
  "Healthy vs Incipient TB baseline",
  "Healthy vs Incipient TB followup",
  "Healthy vs Subclinical TB baseline",
  "Healthy vs Subclinical TB followup",
  "Active TB vs Incipient TB baseline",
  "Active TB vs Incipient TB followup",
  "Active TB vs Subclinical TB baseline",
  "Active TB vs Subclinical TB followup"
)




# 5) ROC Curves & Forest plot-----------------------------------------------------------


# Run roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

clinical_treat <- clinical

roc_func_res <- roc_func(outcome = outcome)
res_table <- roc_func_res$res_table
forestplot_res_table <- roc_func_res$forestplot_res_table

write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))


  
  
  
disease_roc_subset <- c("Healthy vs Active TB",
  "Healthy vs Incipient TB baseline",
  "Healthy vs Incipient TB followup",
  "Healthy vs Subclinical TB baseline",
  "Healthy vs Subclinical TB followup")


timepoint_roc_subset <- c("Active TB vs Incipient TB baseline",
  "Active TB vs Incipient TB followup",
  "Active TB vs Subclinical TB baseline",
  "Active TB vs Subclinical TB followup")


disease_legend_nrow = 3
timepoint_legend_nrow = 2

roc_plot_width = 3200
roc_plot_height = 3500
  


# disease plot
disease_roc <- disease_roc_plot_func(legend_nrow = 2)
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

#timepoint plot
timepoint_roc <- timepoint_roc_plot_func(legend_nrow = 2)
ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome, ".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

#forestplot
panel_forest <- forestplot_func()
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1CP"), 
                       size = 12, hjust = 0, x = 0)
  )
ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  
  
  
  
  
# GSE79362 -----------------------------------------------------------------------------------------------------------------------------------
#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c("gene_annot",lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE79362"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))


this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

if(!exists(file.path(this.accession.res.dir, "figures"))) dir.create(file.path(this.accession.res.dir, "figures"))

## 1) Load public data ------------------------------------------------------
#Ran in HPC
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE79362", "file=GSE79362_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
# write.table(tbl,"/shared/homes/165861/RBMB/TB/data/public/GSE79362/raw_counts.txt")

gse=getGEO(filename="GSE79362_series_matrix.txt.gz", getGPL = FALSE)

raw_counts <- read.delim("raw_counts.txt", sep = " ")
raw_metadata <- gse@phenoData@data

colnames(raw_metadata) <- make.names(colnames(raw_metadata))

colnames(raw_metadata) 
raw_clinical <- raw_metadata[,c("age.ch1",
                                "gender.ch1",
                                "group.ch1",
                                "previous.diagnosis.of.tb.ch1", #collecting blood and culturing the Mycobacterium tuberculosis (MTB) bacteria
                                "qft.ch1")] #blood test used to diagnose tuberculosis (TB) by measuring the amount of interferon-gamma (IFN-γ) released by a person's immune cells when exposed to specific TB antigens

colnames(raw_clinical) <- c("age",
                            "sex",
                            "disease",
                            "previous_tb_diagnosis",
                            "qft")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))
table(clinical$disease)


## 2) Get groups to be compared --------------------------------------------

clinical[which(clinical$disease == "control (non-progressor)"),"group"] <- "Healthy" # infected but stayed healthy
clinical[which(clinical$disease == "case (TB progressor)"),"group"] <- "Active TB" # infected and developed intrathoracic TB

if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}

table(clinical$group)

counts_vst <- vst(as.matrix(raw_counts))
counts_vst <- counts_vst[,row.names(clinical)]

write.table(counts_vst, file.path("counts_vst.txt"))
write.table(clinical, file.path("clinical.txt"))


counts_norm <- counts_vst

## 3) Mean of z-scored expression ---------------------------

gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 

## 3.1) Boxplot ---------------------------
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text.y = element_text(size = 24),
                    axis.text.x = element_text(size = 20),
                    title = element_text(size = 20),
                    legend.position = "None") 


# Specify order that we want the x-axis variables in for the plot
x_order <- c("Healthy", "Active TB")


# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot <- boxplot_all

outcome = "TB" #no outcome data available so just use 'TB' as placeholder
boxplotfig <- boxplot_func(outcome = outcome) 

boxplotfig <- boxplotfig + labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Raw counts were VST normalised (DESeq2 1.42.1)\n",
       "P values from Mann-Whitney U test shown")) 


ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3000,
       height = 2700,
       units = "px" )



## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
 c("Healthy", "Active TB")
)


# Define for plotting
comparison_levels <- c(
"Healthy vs Active TB")


comparison_plotlabel_levels <- comparison_levels

# 5) ROC Curves & Forest plot-----------------------------------------------------------
# Run roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

clinical_treat <- clinical

roc_func_res <- roc_func(outcome = outcome)
res_table <- roc_func_res$res_table
forestplot_res_table <- roc_func_res$forestplot_res_table

write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))


disease_roc_subset <- c(
  "Healthy vs Active TB")

roc_plot_width = 3000
roc_plot_height = 3200


# disease plot
disease_roc <- disease_roc_plot_func(legend_nrow = 1)
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

# no timepoint data 

#forestplot
panel_forest <- forestplot_func()
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1CP"), 
                       size = 12, hjust = 0, x = 0)
  )
ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  






# =========================================================================================================== #
# ========================= 2. TB MULTIDRUG RESISTANT AMD DRUG SUSCEPTIBLE DATASETS ============================
# =========================================================================================================== #


# GSE147690 (German Identification cohort) ----------------------------------------------------------------------------------------------------------------------------------------
#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE147690"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)


  

#1) 50 DS-GIC and 30 MDR-GIC (2013-2016)
#2) 38 DS-GVC and 32 MDR-GVC (2015 - 2018)
#2) 52 MDR-RVC (2015 - 2017)


# 329 Samples of 80 patients with non MDR and MDR TB, repeated measurements after therapy start:
# (14 days, culture conversion, sputum conversion, 6 months, 10 months, 15 months and therapy end); 
# one time measurements of 10 healthy controls

# ## 1) Load public data -------------------------------------------------------------
# #One colour Agilent Intensity data
# data_dir <-file.path(my_directory,"TB", "data", "public", this.accession.no,"GSE147690_RAW")
# 
# #Get list of raw files
# files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
# 
# 
# #Get GSE/sample names
# targets <- data.frame(
#   FileName = basename(files),
#   SampleName = sub("_.*", "", basename(files)) , # Extract GSM IDs
#   SampleName_sg = sub("^GSM\\d+_(.*)\\.txt\\.gz$", "\\1", basename(files))  # Extract GSM IDs
# 
# )
# 
# 
# 
# 
# #The green.only argument tells read.maimages() to not expect a red colour (because its one-colour Agilent) output an EList object instead an RGList. The raw intensities will be stored in the E component of the data object.
# RG <- read.maimages(files, source = "agilent", green.only = TRUE)
# 
# #Make count matrix
# exprs_matrix <- RG$E #44495
# dim(exprs_matrix)
# rownames(exprs_matrix) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(exprs_matrix) <- targets$SampleName
# 
# View(RG$genes)
# #Note there are two different types of controls
# controls_to_remove <- RG$genes[which(RG$genes$ControlType != 0), "SystematicName"]#1377 controls to remove
# raw_counts <- exprs_matrix[-which(row.names(exprs_matrix) %in% controls_to_remove),] #43118 genes left
# dim(raw_counts)
# 
# gset <- getGEO("GSE147690", GSEMatrix =TRUE, getGPL=FALSE)
# gse=gset$GSE147690_series_matrix.txt.gz
# 
# raw_metadata <- gse@phenoData@data #329 samples
# 
# 
# #Loading the normalised data
# normalised_expr <- read.csv(file.path(my_directory,"TB", "data", "public", this.accession.no,
#                                       "GSE147690_Kopie_von_GA_Agilent_one_color_matrix_Heyckendorf_GIC_Normalized_data.csv"))
# normalised_expr <- as.matrix(normalised_expr[,-1])
# 
# raw_metadata$title <- make.names(raw_metadata$title)
# dim(normalised_expr)
# colnames(normalised_expr) == raw_metadata$title
# 
# row.names(normalised_expr) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(normalised_expr)  == targets$SampleName_sg #same order between normalised counts and target file, use GSM ids intsead of SG
# colnames(normalised_expr)  == raw_metadata$title
# 
# 
# colnames(normalised_expr) <- raw_metadata$geo_accession
# 
# 
# normalised_counts <- normalised_expr[-which(row.names(normalised_expr) %in% controls_to_remove),] #43118 genes left
# dim(normalised_counts)
# 
# 
# # Save files
# write.csv(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.csv")))
# write.csv(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.csv")))
# write.csv(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.csv")))
# saveRDS(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.rds" )))
# saveRDS(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.rds" )))
# saveRDS(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.rds" )))

## 1) Load public data -------------------------------------------------------------
raw_metadata <- readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_metadata.rds")))
raw_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_counts.rds")))
normalised_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_normalised_counts.rds")))

raw_metadata$data_processing
#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"


#Gene annotations
gpl <- getGEO("GPL13497", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("REFSEQ", "ENSEMBL_ID", "GENE_SYMBOL")]


colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "contact_country",
                                "disease.state.ch1",
                                "source_name_ch1")]

patient <- str_extract(raw_clinical$source_name_ch1, "^[A-Za-z]+\\d+")
patient_id <-str_extract(patient, "\\d+")
patient <- str_extract(raw_clinical$source_name_ch1, "^[A-Za-z]+")
day <- str_extract(raw_clinical$source_name_ch1, "(?<=Day ?)(\\d+)")

raw_clinical<-cbind(raw_clinical,patient,patient_id,day)
raw_clinical$patient_id <- as.numeric(raw_clinical$patient_id)
raw_clinical$day <- as.numeric(raw_clinical$day)
View(raw_clinical[order(raw_clinical$patient, raw_clinical$patient_id),])

#number of non MDR patients = 50
length(unique(raw_clinical[which(raw_clinical$disease.state.ch1 == "non MDR"), "patient_id"]))

#number of MDR patients=30
length(unique(raw_clinical[which(raw_clinical$disease.state.ch1 == "MDR"), "patient_id"]))

#Note patients 40, 68 and 33 are missing day info


#1) 50 DS-GIC and 30 MDR-GIC (2013-2016)
#2) 28 DS-GVC and 32 MDR-GVC (2015 - 2018)
#2) 52 MDR-RVC (2015 - 2017)

# Study visits were performed at 
## (ideally) before treatment initiation, 
## at 14 days of therapy, 
## at the times of smear conversion and following culture conversion (not available in the MDR-RVC), 
## at 6 months and/or therapy end in patients with DS-TB, 
## and additionally at 10, 15 and 20 months of therapy in patients with MDR-TB. 
# After completion of 4 weeks of therapy, an additional study visit was performed in patients from the MDR-RVC. 
# All patients completed 12 months of evaluation following the end of therapy to capture disease recurrence. 



raw_clinical$months <- raw_clinical$day/30


colnames(raw_clinical)[1:4] <- c("sample_id",
                                 "contact_country",
                                 "disease",
                                 "full_id")

clinical <- raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))


## 2) Get groups to be compared
#MDR = Multi-drug resistant TB
#DS-TB = drug-susceptible TB
table(clinical$disease)
table(clinical$group)

clinical[which(clinical$disease == "Healthy Controls"), "disease"] <- "Healthy"
clinical[which(clinical$disease == "MDR"), "disease"] <- "DS_TB"
clinical[which(clinical$disease == "non MDR"), "disease"] <- "MDR_TB"

clinical$group <- clinical$disease

if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)

#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"
any(is.na(raw_counts))

# counts_voom <- voom(as.matrix(raw_counts))
# counts_norm <- counts_voom

counts_norm <- as.matrix(normalised_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

# library(reshape2)
# ggplot(melt(raw_counts), aes(x = value)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Density Plot of Microarray Data",
#        x = "Expression Value",
#        y = "Density") +
#   theme_minimal()

hist(counts_norm)

GSE147690_clinical <- clinical
GSE147690_counts_norm <- counts_norm

write.table(GSE147690_counts_norm, file.path("counts_norm.txt"))
write.table(GSE147690_clinical, file.path("clinical.txt"))


## 3) Mean of z-scored expression ---------------------------
  gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

# For each sample, get the mean standardized expression of genes in the 4-gene signature
mean_sig_zscore <- mean_zscore_func()

## 3.1) Boxplot ---------------------------
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    disease = clinical$disease,
                                    day = clinical$day,
                                    months = clinical$months))

boxplot_all[which(boxplot_all$disease == "Healthy"), "months"] <- 0


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.text = element_text(size = 18)) 


outcome = "TB"

#the bocplot gunction is for boxplots but this is geom_smooth. just make it manually
#geom_point, split by disease
boxplotfig <- ggplot(boxplot_all, aes(
  x = as.numeric(months),
  y = as.numeric(boxplot_all[,1]),
  group = as.factor(disease))) +
  
  theme_bw()+
  boxplot_theme +
  geom_point(aes(colour=disease)) +
  geom_smooth(method = "loess", se = FALSE, aes(colour=disease))+
  
theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       color = "Disease", #legend title
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy", "\n",
                        "n=", nrow(boxplot_all),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from the original study: 'Data has been log-transformed and\nnormalized (percentile shift) within the Agilent GeneSpring V13 software'","\n",
       "Curve generated using loess method (non-parametric local regression)"
)) +
  scale_color_manual(values = c("Healthy" = "#7CAE00",
                                "MDR_TB" = "#00BFC4",
                                "DS_TB" = "#F8766D"))+
  ylab (label = "Signature Score") +
  xlab (label = "Months of Treatment")



ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600, 
       height = 3250, 
       units = "px" )


## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Healthy", "MDR_TB"),
  c('Healthy', "DS_TB"),
  c('DS_TB', "MDR_TB")
  
)


# Define for plotting
comparison_levels <- c(
  "Healthy vs MDR_TB",
  "Healthy vs DS_TB",
  "DS_TB vs MDR_TB"
)


comparison_plotlabel_levels <- comparison_levels

#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

#Only get T0 and Healthy
clinical_treat <- clinical[which(clinical$months == 0| clinical$group == "Healthy"),]
clinical_treat[which(clinical_treat$group == "Healthy"), "months"] <- 0

roc_func_res <- roc_func(outcome = outcome)
res_table <- roc_func_res$res_table
forestplot_res_table <- roc_func_res$forestplot_res_table

write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))






disease_roc_subset <- comparison_plotlabel_levels

roc_plot_width = 3200
roc_plot_height = 3500



# disease plot
disease_roc <- disease_roc_plot_func(legend_nrow = 2)

disease_roc <- disease_roc +   labs(
    title = paste0("ROC - Control vs TB", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))
  
  
    
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")


#forestplot
panel_forest <- forestplot_func()
panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was \n used as a proxy"),
                       size = 12, hjust = 0, x = 0))
    
ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  
  
  
# 
# 
# # Loop through each pairwise comparison
# for (pair in pairwise_comparisons) {
#   
#   group1 <- pair[1]
#   group2 <- pair[2]
#   
#   subset_clinical <- clinical[which(clinical$months == 0| clinical$group == "Healthy"),]
#   subset_clinical[which(subset_clinical$group == "Healthy"), "months"] <- 0
#   
#   # Subset data to omly include the 2 rgroups of interest
#   subset_clinical <- subset_clinical[subset_clinical$group %in% c(group1,group2),]
#   subset_counts <- counts_norm[, row.names(subset_clinical)]
#   
#   subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
#   
#   # Mean-zscore
# # Get the gene IDs instead of HGNCs
# gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))
# 
# signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
# gene_set_list <- c(signature_geneid)
# 
# 
# if(length(signature_geneid) < 4){ 
#   print("Missing gene in signature after genone_annot conversion")
#   stop() }
# 
# # For each sample, get the mean standardized expression of genes in the 4-gene signature
# 
# gene_set <- subset_counts[gene_set_list,]
# gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
#       # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
#       # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
# gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
#       #take the mean of the scaled+centred data for each gene and use this as each sample's score
#       #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
# mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
# colnames(mean_sig_zscore) <- "score"
# 
# all(row.names(mean_sig_zscore) == row.names(subset_clinical))
# 
#   
#   glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
#   
#   table(glm_data$Group)
#   
#   
#   glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
#   glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
#   
#   test_probs <- predict(glm_model, type = "response")
#   
#   roc_obj <- roc(glm_data$Group, test_probs)
#   
#   plot(roc_obj)
#   auc(roc_obj)
#   auc_ci <- ci.auc(roc_obj)
#   
#   #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
#   optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
#   
#   
#   if(nrow(optimal_threshold_coords) > 1) {
#     optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
#   }
#   
#   
# 
#   # Sensitivity confidence interval (at optimal specificity)
#   ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
#   
#   # Specificity confidence interval (at optimal sensitivity)
#   ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
#   
#   res_current <-cbind(
#     comparison = paste0(group1,"vs",group2),
#     samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
#     samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
#     auc = auc(roc_obj),
#     ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
#     sensitivity = optimal_threshold_coords$sensitivity, 
#     specificity = optimal_threshold_coords$specificity
#     
#   )
#   
#   res_table <- rbind(res_table, res_current)
#   
#   forestplot_res_table <- rbind(forestplot_res_table, 
#                                 cbind(comparison = paste0(group1,"vs",group2),
#                                       auc = auc(roc_obj),
#                                       auc_ci_low = as.numeric(auc_ci[1]),
#                                       auc_ci_high = as.numeric(auc_ci[3]),
#                                       
#                                       sensitivity = optimal_threshold_coords$sensitivity, 
#                                       sensitivity_ci_low = ci_sens[, "2.5%"],
#                                       sensitivity_ci_high = ci_sens[, "97.5%"],
#                                       
#                                       specificity = optimal_threshold_coords$specificity,
#                                       specificity_ci_low = ci_spec[, "2.5%"],
#                                       specificity_ci_high = ci_spec[, "97.5%"])
#   )
#   
#   roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
#   
# } # close pair
#   
# write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
# write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))
# 
# 
# 
# # 5) ROC Curves & Forest plot-----------------------------------------------------------
# 
# # Convert ROC data to a format suitable for ggplot
# roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
#   data.frame(
#     TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
#     FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
#     Comparison = comparison,
#     auc = rev(roc_objects[[comparison]]$auc)
#   )
# }))
# 
# roc_data$Comparison <- factor(roc_data$Comparison, levels = comparison_levels)
# 
# roc_data$Comparison_plotlabel <- roc_data$Comparison
# levels(roc_data$Comparison_plotlabel) <-  comparison_plotlabel_levels
# 
# roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]
# 
# roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
#                           round(roc_data$auc, 2), " (", roc_data$ci, ")")
# 
# 
# 
# # Disease plot
# disease_roc_data <- roc_data
# 
# disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
#   geom_line(size = 1.2) +
#   theme_bw() +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
#   guides(colour = guide_legend(nrow = disease_legend_nrow)) +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         axis.title = element_text(size = 24),
#         axis.text = element_text(size = 24),
#         legend.text = element_text(size = 16),
#         title = element_text(size = 20)) +
# 
#   labs(
#     title = paste0("ROC - Control vs TB", " (", outcome, ")"),
#     x = "FPR (1 - Specificity)",
#     y = "TPR (Sensitivity)",
#     color = "Comparison",
#        caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
#                         "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))
#     
#     
# ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
#        width = roc_plot_width, 
#        height = roc_plot_height, 
#        units = "px")
# 
# ## Forest plot -------
# res_table <- forestplot_res_table
#   
# res_table$comparison <- factor(res_table$comparison, levels = comparison_levels)
#   
# res_table[,-1] <- lapply(res_table[,-1], as.numeric)
# 
# 
# forestplot_theme <- theme(axis.title = element_text(size =  16),
#                             axis.text = element_text(size = 14),
#                             legend.position = "None") 
# 
# 
# auc_plot <- res_table %>% 
#     ggplot(aes(y = comparison)) + 
#     theme_bw() +
#     geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
#     geom_linerange(aes(xmin=auc_ci_low, 
#                        xmax=auc_ci_high, 
#                        color = comparison),
#                    size = 1) +
#           scale_y_discrete(labels = comparison_plotlabel_levels)+
#     forestplot_theme +
#     ylab(NULL)+
#     xlab("AUC")+
#     coord_cartesian( xlim=c(0.1, 1))
# 
# 
#   sens_plot <- res_table %>% 
#     ggplot(aes(y = comparison)) + 
#     theme_bw() +
#     geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
#     geom_linerange(aes(xmin=sensitivity_ci_low, 
#                        xmax=sensitivity_ci_high, 
#                        color = comparison),
#                    size = 1) +
#         scale_y_discrete(labels = comparison_plotlabel_levels)+
#     forestplot_theme +
#     ylab(NULL)+
#     xlab("Sensitivity")+
#     coord_cartesian( xlim=c(0.1, 1))
#   
#   spec_plot <- res_table %>% 
#     ggplot(aes(y = comparison)) + 
#     theme_bw() +
#     geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
#     geom_linerange(aes(xmin=specificity_ci_low, 
#                        xmax=specificity_ci_high, 
#                        color = comparison),
#                    size = 1) +
#         scale_y_discrete(labels = comparison_plotlabel_levels)+
#     forestplot_theme +
#     ylab(NULL)+
#     xlab("Specificity")+
#     coord_cartesian( xlim=c(0.1, 1))
#   
#   
#   panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
#                             ncol = 1,
#                             nrow = 3)
#   
#   panel_forest <- annotate_figure(
#     panel_forest,
#     top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
#     bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
#                               "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
#                         "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"),
#                        size = 12, hjust = 0, x = 0)
#   )
#   
#     
#   ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
#          width = 15, height = 20, units = "cm",   bg = "white"  )
# 
#   
#   
#   
  
# GSE147689 (German validation cohort) ----------------------------------------------------------------------------------------------------------------------------------------
#free up space and remove objects from previous gse (except for the functions I made and gene_annot)

this.accession.no <- "GSE147689"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)


# 335 Samples of 60 patients with non MDR and MDR TB, 
#repeated measurements after therapy start
## (14 days, culture conversion, sputum conversion, 6 months, 10 months, 15 months and therapy end)


# ## 1) Load public data -------------------------------------------------------------
# 
# #One colour Agilent Intensity data
# data_dir <-file.path(my_directory,"TB", "data", "public", this.accession.no,"GSE147689_RAW")
# 
# #Get list of raw files
# files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
# 
# 
# #Get GSE/sample names
# targets <- data.frame(
#   FileName = basename(files),
#   SampleName = sub("_.*", "", basename(files)) , # Extract GSM IDs
#   SampleName_sg =sub("^GSM\\d+_(.*)\\.txt\\.gz$", "\\1", basename(files))  # Extract GSM IDs
#   
# )
# 
# 
# 
# 
# #The green.only argument tells read.maimages() to not expect a red colour (because its one-colour Agilent) output an EList object instead an RGList. The raw intensities will be stored in the E component of the data object.
# RG <- read.maimages(files, source = "agilent", green.only = TRUE)
# 
# #Make count matrix
# exprs_matrix <- RG$E #44495
# rownames(exprs_matrix) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(exprs_matrix) <- targets$SampleName
# 
# View(RG$genes)
# #Note there are two different types of controls
# controls_to_remove <- RG$genes[which(RG$genes$ControlType != 0), "SystematicName"]#1377 controls to remove
# raw_counts <- exprs_matrix[-which(row.names(exprs_matrix) %in% controls_to_remove),] #43118 genes left
# 
# gset <- getGEO("GSE147689", GSEMatrix =TRUE, getGPL=FALSE)
# gse=gset$GSE147689_series_matrix.txt.gz
# 
# raw_metadata <- gse@phenoData@data #335 samples
# 
# 
# #Loading the normalised data
# normalised_expr <- read.csv(file.path(my_directory,"TB", "data", "public", this.accession.no,
#                                       "GSE147689_Kopie_von_GA_Agilent_one_color_matrix_Heyckendorf_GVC_Normalized_data.csv"))
# 
# dim(normalised_expr)
# normalised_expr <- as.matrix(normalised_expr[,-1])
# row.names(normalised_expr) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(normalised_expr)  == targets$SampleName_sg #same order between normalised counts and target file, use GSM ids intsead of SG
# colnames(normalised_expr)  == raw_metadata$title
# 
# colnames(normalised_expr) <- raw_metadata$geo_accession
# 
# normalised_counts <- normalised_expr[-which(row.names(normalised_expr) %in% controls_to_remove),] #43118 genes left
# 
# Save files
# write.csv(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.csv")))
# write.csv(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.csv")))
# write.csv(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.csv")))
# saveRDS(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.rds" )))
# saveRDS(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.rds" )))
# saveRDS(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.rds" )))

# ## 1) Load public data -------------------------------------------------------------

raw_metadata <-readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_metadata.rds")))
raw_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_counts.rds")))
normalised_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_normalised_counts.rds")))

raw_metadata$data_processing
#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"


#Gene annotations
gpl <- getGEO("GPL13497", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("REFSEQ", "ENSEMBL_ID", "GENE_SYMBOL")]


colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "contact_country",
                                "disease.state.ch1",
                                "source_name_ch1")]

patient_id <-str_extract(raw_clinical$source_name_ch1, "(?<=Pat_Valid_)\\d+")
day <- str_extract(raw_clinical$source_name_ch1, "(?<=Day ?)(\\d+)")

raw_clinical<-cbind(raw_clinical,patient_id,day)
raw_clinical$patient_id <- as.numeric(raw_clinical$patient_id)
raw_clinical$day <- as.numeric(raw_clinical$day)
View(raw_clinical[order(raw_clinical$disease.state.ch1, raw_clinical$patient_id),])

#number of non MDR patients = 28
length(unique(raw_clinical[which(raw_clinical$disease.state.ch1 == "Non MDR"), "patient_id"]))

#number of MDR patients= 34
length(unique(raw_clinical[which(raw_clinical$disease.state.ch1 == "MDR"), "patient_id"]))

#1) 50 DS-GIC and 30 MDR-GIC (2013-2016)
#2) 28 DS-GVC and 32 MDR-GVC (2015 - 2018)
#2) 52 MDR-RVC (2015 - 2017)

# Study visits were performed at 
## (ideally) before treatment initiation, 
## at 14 days of therapy, 
## at the times of smear conversion and following culture conversion (not available in the MDR-RVC), 
## at 6 months and/or therapy end in patients with DS-TB, 
## and additionally at 10, 15 and 20 months of therapy in patients with MDR-TB. 
# After completion of 4 weeks of therapy, an additional study visit was performed in patients from the MDR-RVC. 
# All patients completed 12 months of evaluation following the end of therapy to capture disease recurrence. 



raw_clinical$months <- raw_clinical$day/30


colnames(raw_clinical)[1:4] <- c("sample_id",
                                 "contact_country",
                                 "disease",
                                 "full_id")

clinical <- raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))



## 2) Get groups to be compared --------------------------------------------
#MDR = Multi-drug resistant TB
#DS-TB = drug-susceptible TB
table(clinical$disease)
table(clinical$group)

clinical[which(clinical$disease == "MDR"), "disease"] <- "DS_TB"
clinical[which(clinical$disease == "Non MDR"), "disease"] <- "MDR_TB"

clinical$group <- clinical$disease

if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)

#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"

# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(normalised_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

# library(reshape2)
# ggplot(melt(raw_counts), aes(x = value)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Density Plot of Microarray Data",
#        x = "Expression Value",
#        y = "Density") +
#   theme_minimal()

# hist(raw_counts)

GSE147689_clinical <- clinical
GSE147689_counts_norm <- counts_norm

write.table(GSE147689_counts_norm, file.path("counts_norm.txt"))
write.table(GSE147689_clinical, file.path("clinical.txt"))



## 3) Mean of z-scored expression ---------------------------

  gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- counts_norm[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

if(all(row.names(mean_sig_zscore) != row.names(clinical))){ stop("Row names do not match", this.accession.no)}



## 3.1) Boxplot ---------------------------
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    disease = clinical$disease,
                                    day = clinical$day,
                                    months = clinical$months))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.text = element_text(size = 18)) 


outcome = "TB"
#geom_point, split by disease
boxplotfig <- ggplot(boxplot_all, aes(
  x = as.numeric(months),
  y = as.numeric(boxplot_all[,1]),
  group = as.factor(disease))) +
  
  theme_bw()+
  boxplot_theme +
  geom_point(aes(colour=disease)) +
  geom_smooth(method = "loess", se = FALSE, aes(colour=disease))+
  
theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       color = "Disease", #legend title
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy", "\n",
                        "n=", nrow(boxplot_all),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from the original study: 'Data has been log-transformed and\nnormalized (percentile shift) within the Agilent GeneSpring V13 software'","\n",
       "Curve generated using loess method (non-parametric local regression)"
)) +
  scale_color_manual(values = c("Healthy" = "#7CAE00",
                                "MDR_TB" = "#00BFC4",
                                "DS_TB" = "#F8766D"))+
  ylab (label = "Signature Score") +
  xlab (label = "Months of Treatment")



ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600, 
       height = 3250, 
       units = "px" )




## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("DS_TB", 'MDR_TB')
)


# Define for plotting
comparison_levels <- c(
  "DS_TBvsMDR_TB"
)


comparison_plotlabel_levels <- c(
  "DS_TB vs MDR_TB"
)


#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

disease_legend_nrow = 2
roc_plot_width = 3200
roc_plot_height = 3500


# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  subset_clinical <- clinical[which(clinical$months == 0),]

  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- subset_clinical[subset_clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # Mean-zscore
# Get the gene IDs instead of HGNCs
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
gene_set_list <- c(signature_geneid)

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  

  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1,"vs",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
} # close pair
  
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))



# 5) ROC Curves & Forest plot-----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison, levels = comparison_levels)

roc_data$Comparison_plotlabel <- roc_data$Comparison
levels(roc_data$Comparison_plotlabel) <-  comparison_plotlabel_levels

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")



# Disease plot
disease_roc_data <- roc_data

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = disease_legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +

  labs(
    title = paste0("ROC - Control vs TB", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))
    
    
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

## Forest plot -------
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = comparison_levels)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 


auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
          scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
        scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
        scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"),
                       size = 12, hjust = 0, x = 0)
  )
  
    
  ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  
  

# COMBINE GERMAN IDENTIFICATION + VALIDATIONS COHORTS) ----------------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE147690_GSE147689"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)



german_clinical <-rbind(GSE147689_clinical, GSE147690_clinical[,-5])
german_counts <- cbind(GSE147689_counts_norm, GSE147690_counts_norm)

colnames(german_counts) == row.names(german_clinical)

table(german_clinical$group)
table(german_clinical$months)

german_clinical[which(german_clinical$disease == "Healthy"), "months"] <- 0


## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(german_counts)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }


# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- german_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

if(all(row.names(mean_sig_zscore) != row.names(clinical))){ stop("Row names do not match", this.accession.no)}


## 3.1) Boxplot ---------------------------
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    day = german_clinical$day,
                                    months = german_clinical$months,
                                    group = german_clinical$group))



boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.text = element_text(size = 18)) 


x_order <- c("Healthy", "DS_TB", "MDR_TB")


my_comparisons <- combn(unique(german_clinical$group), 2, simplify = FALSE)


boxplot_all$group <- factor(boxplot_all$group, levels = x_order)
boxplot_all$score <- as.numeric(boxplot_all$score)


outcome = "TB"
#geom_point, split by disease
boxplotfig<- ggplot(boxplot_all, aes(
  x = as.numeric(months),
  y = as.numeric(boxplot_all[,1]),
  group = as.factor(group))) +
  
  theme_bw()+
  boxplot_theme +
  geom_point(aes(colour=group)) +
  geom_smooth(method = "loess", se = FALSE, aes(colour=group))+
  scale_color_manual(values = c("Healthy" = "#7CAE00",
                                "MDR_TB" = "#00BFC4",
                                "DS_TB" = "#F8766D"))+
  # stat_pvalue_manual(stat.table.gsva,
  #                    label = "p",
  #                    tip.length = 0.01,
  #                    size = 7)+
  # 
  # # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  # stat_summary(fun.y = mean, fill = "red",
  #              geom = "point", shape = 21, size =4,
  #              show.legend = TRUE) +
  # 


theme(axis.text.x = element_text(size = 18))+
  labs(title = "Signature Analysis: GSE147689 & GSE147690",
       color = "Disease", #legend title
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy", "\n",
                        "n=", nrow(boxplot_all),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from the original study: 'Data has been log-transformed and\nnormalized (percentile shift) within the Agilent GeneSpring V13 software'","\n",
       "Curve generated using loess method (non-parametric local regression)"
)) +
  ylab (label = "Signature Score") +
  xlab (label = "Months of Treatment")


ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3700, 
       height = 3250, 
       units = "px" )




## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Healthy", "MDR_TB"),
  c('Healthy', "DS_TB"),
  c('DS_TB', "MDR_TB")
  
)


# Define for plotting
comparison_levels <- c(
  "HealthyvsMDR_TB",
  "HealthyvsDS_TB",
  "DS_TBvsMDR_TB"
)


comparison_plotlabel_levels <- c(
  "Healthy vs MDR_TB",
  "Healthy vs DS_TB",
  "DS_TB vs MDR_TB"
)


#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

disease_legend_nrow = 2
roc_plot_width = 3200
roc_plot_height = 3500

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  subset_clinical <- german_clinical[which(german_clinical$months == 0),]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- subset_clinical[subset_clinical$group %in% c(group1,group2),]
  subset_counts <- german_counts[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(subset_counts)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }


# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"



all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  

  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1,"vs",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
} # close pair
  
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))



# 5) ROC Curves & Forest plot-----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison, levels = comparison_levels)

roc_data$Comparison_plotlabel <- roc_data$Comparison

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease plot
disease_roc_data <- roc_data

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = disease_legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +

  labs(
    title = paste0("ROC - Control vs TB", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))
    
    
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")

## Forest plot -------
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = comparison_levels)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 


auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
          scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
        scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
        scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"),
                       size = 12, hjust = 0, x = 0)
  )
  
    
  ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  
  
  
  
  
  
  
  
## 4) ROUNDED Validation  ------------------------------------------------------
# CREATE TIMEPOINTS

round_half <- function(x) {
  round(x * 2) / 2
}
table(round_half(german_clinical$months))

table(round(german_clinical$months))
# Can see from this that there are the below no. of samples available
# T0: 105
# 14 days: 113
# ~ 6 months: 5 months (15), 6 months(64), 7 months(24)
# ~ 10 months: 9 motnhs(9), 10 months (28), 11 months (16)
# ~ 15 months: 14 months (2), 15 months (22), 16 months (14)
# ~ 20 months: 19 months (3), 20 months(17), 21 months (15)

german_clinical$months_round <- round(german_clinical$months)
german_clinical$months_round_half <- round_half(german_clinical$months)

german_clinical[which(german_clinical$months == 0), "group_round"] <- "0_Days"

german_clinical[which(german_clinical$months_round_half == 0.5), "group_round"] <- "14_Days"

german_clinical[which(german_clinical$months_round_half == 5.5|
                        german_clinical$months_round_half == 6|
                        german_clinical$months_round_half == 6.5), "group_round"] <- "6_Mths"

german_clinical[which(german_clinical$months_round_half == 9.5|
                        german_clinical$months_round_half == 10|
                        german_clinical$months_round_half == 10.5), "group_round"] <- "10_Mths"


german_clinical[which(german_clinical$months_round_half == 14.5|
                        german_clinical$months_round_half == 15|
                        german_clinical$months_round_half == 15.5), "group_round"] <- "15_Mths"


german_clinical[which(german_clinical$months_round_half == 19.5|
                        german_clinical$months_round_half == 20|
                        german_clinical$months_round_half == 20.5), "group_round"] <- "20_Mths"

german_clinical2 <- german_clinical[!is.na(german_clinical$group_round),]
german_clinical2$group_new <- paste0(german_clinical2$disease,"_", german_clinical2$group_round)

german_counts2 <- german_counts[,row.names(german_clinical2)]

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest

table(german_clinical2$group_new)  #10 and 15 months dont have enough samples
pairwise_comparisons <- list(
  c("Healthy_0_Days", "MDR_TB_0_Days"),
  c('Healthy_0_Days', "DS_TB_0_Days"),
  c('DS_TB_0_Days', "MDR_TB_0_Days"),
  
  c('MDR_TB_0_Days', "MDR_TB_14_Days"),
  c('MDR_TB_0_Days', "MDR_TB_6_Mths"), 

  c('DS_TB_0_Days', "DS_TB_14_Days"),
  c('DS_TB_0_Days', "DS_TB_6_Mths"),
  c('DS_TB_0_Days', "DS_TB_10_Mths"),
  c('DS_TB_0_Days', "DS_TB_20_Mths")
)


# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- german_clinical2[german_clinical2$group_new %in% c(group1,group2),]
  subset_counts <- german_counts2[, row.names(subset_clinical)]
  
  subset_clinical$group_new <- factor(subset_clinical$group_new, levels = c(group1, group2))
  
  
## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(subset_counts)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }


# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"



all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], group_new = subset_clinical$group_new)
  
  table(glm_data$group_new)
  
  
  glm_data$group_new <- factor(glm_data$group_new, levels = c(group1, group2))
  glm_model <- glm(group_new ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$group_new, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  

  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1," vs ",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$group_new == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$group_new == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1," vs ",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
} # close pair
  
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))



## 5) ROUNDED ROC Curves -----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))


roc_data$Comparison <- factor(roc_data$Comparison)


roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data[which(roc_data$Comparison == "DS_TB_0_Days vs MDR_TB_0_Days" | 
                                     roc_data$Comparison == "Healthy_0_Days vs DS_TB_0_Days" |
                                     roc_data$Comparison == "Healthy_0_Days vs MDR_TB_0_Days" ), ]

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - TB vs Healthy at T0",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))
    
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = 3200, 
       height = 3500, 
       units = "px")


# Timepoint ROC
timepoint_roc_data <- roc_data[which(roc_data$Comparison == "MDR_TB_0_Days vs MDR_TB_14_Days" | 
                                       roc_data$Comparison == "MDR_TB_0_Days vs MDR_TB_6_Mths" |
                                       
                                       roc_data$Comparison == "DS_TB_0_Days vs DS_TB_14_Days" |
                                       roc_data$Comparison == "DS_TB_0_Days vs DS_TB_6_Mths" |
                                       roc_data$Comparison == "DS_TB_0_Days vs DS_TB_10_Mths" |
                                       roc_data$Comparison == "DS_TB_0_Days vs DS_TB_20_Mths"), ]

timepoint_roc <- ggplot(timepoint_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom",
          legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - TB treatment timepoints",
    x = "FPR (1 - Specificity)",
    y = "TPR(Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))


ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome,".png")), 
       width = 3400, 
       height = 3200, 
       units = "px")



## Forest plot -------
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = roc_data$Comparison)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 


auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
          # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
        # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
        # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"),
                       size = 12, hjust = 0, x = 0)
  )
  
    
  ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 20, height = 20, units = "cm",   bg = "white"  )

# Study visits were performed at 
## (ideally) before treatment initiation, 
## at 14 days of therapy, 
## at the times of smear conversion and following culture conversion (not available in the MDR-RVC), 
## at 6 months and/or therapy end in patients with DS-TB, 
## and additionally at 10, 15 and 20 months of therapy in patients with MDR-TB. 
# After completion of 4 weeks of therapy, an additional study visit was performed in patients from the MDR-RVC. 
# All patients completed 12 months of evaluation following the end of therapy to capture disease recurrence. 


  
  
  
  
# GSE147691 (Romanian Validation cohort) ----------------------------------------------------------------------------------------------------------------------------------------

#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE147691"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

#1) 50 DS-GIC and 30 MDR-GIC (2013-2016)
#2) 38 DS-GVC and 32 MDR-GVC (2015 - 2018)
#2) 52 MDR-RVC (2015 - 2017)

# 262 samples of 52 at different measurement times during TB-treatment
# (0 days, 14 days, 28 days, 6 months, 10 months, 15 months and 20 months)

# ## 1) Load public data -------------------------------------------------------------
#One colour Agilent Intensity data
# data_dir <-file.path(my_directory,"TB", "data", "public", this.accession.no,"GSE147691_RAW")
# 
# #Get list of raw files
# files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
# 
# 
# #Get GSE/sample names
# targets <- data.frame(
#   FileName = basename(files),
#   SampleName = sub("_.*", "", basename(files)) , # Extract GSM IDs
#   SampleName_sg = sub("^GSM\\d+_(.*)\\.txt\\.gz$", "\\1", basename(files))  # Extract GSM IDs
# 
# )
# 
# 
# 
# 
# #The green.only argument tells read.maimages() to not expect a red colour (because its one-colour Agilent) output an EList object instead an RGList. The raw intensities will be stored in the E component of the data object.
# RG <- read.maimages(files, source = "agilent", green.only = TRUE)
# 
# #Make count matrix
# exprs_matrix <- RG$E #44495
# dim(exprs_matrix)
# rownames(exprs_matrix) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(exprs_matrix) <- targets$SampleName
# 
# View(RG$genes)
# #Note there are two different types of controls
# controls_to_remove <- RG$genes[which(RG$genes$ControlType != 0), "SystematicName"]#1377 controls to remove
# raw_counts <- exprs_matrix[-which(row.names(exprs_matrix) %in% controls_to_remove),] #43118 genes left
# dim(raw_counts)
# 
# gset <- getGEO("GSE147691", GSEMatrix =TRUE, getGPL=FALSE)
# gse=gset$GSE147691_series_matrix.txt.gz
# 
# raw_metadata <- gse@phenoData@data #262samples
# 
# #Loading the normalised data
# normalised_expr <- read.csv(file.path(my_directory,"TB", "data", "public", this.accession.no,
#                                       "GSE147691_Kopie_von_GA_Agilent_one_color_matrix_Heyckendorf_RVC_Normalized_data.csv"))
# normalised_expr <- as.matrix(normalised_expr[,-1])
# 
# raw_metadata$title <- make.names(raw_metadata$title)
# dim(normalised_expr)
# colnames(normalised_expr) == raw_metadata$title
# 
# row.names(normalised_expr) <- RG$genes$SystematicName  #see View(RG$genes)
# colnames(normalised_expr)  == targets$SampleName_sg #same order between normalised counts and target file, use GSM ids intsead of SG
# colnames(normalised_expr)  == raw_metadata$title
# 
# 
# colnames(normalised_expr) <- raw_metadata$geo_accession
# 
# 
# normalised_counts <- normalised_expr[-which(row.names(normalised_expr) %in% controls_to_remove),] #43118 genes left
# dim(normalised_counts)
# 
# 
# # Save files
# write.csv(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.csv")))
# write.csv(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.csv")))
# write.csv(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.csv")))
# saveRDS(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.rds" )))
# saveRDS(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.rds" )))
# saveRDS(normalised_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_normalised_counts.rds" )))

  
# ## 1) Load public data -------------------------------------------------------------
raw_metadata <-readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_metadata.rds")))
raw_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_raw_counts.rds")))
normalised_counts <-  readRDS(file.path(getwd(), paste0(this.accession.no,"_normalised_counts.rds")))

raw_metadata$data_processing
#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"


#Gene annotations
gpl <- getGEO("GPL13497", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("REFSEQ", "ENSEMBL_ID", "GENE_SYMBOL")]


colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "disease.state.ch1",
                                "source_name_ch1")]

patient_id <-str_extract(raw_clinical$source_name_ch1, "(?<=Pat_RVC_)\\d+")
day <- str_extract(raw_clinical$source_name_ch1, "(?<=Day ?)(\\d+)")

raw_clinical<-cbind(raw_clinical,patient_id,day)
raw_clinical$patient_id <- as.numeric(raw_clinical$patient_id)
raw_clinical$day <- as.numeric(raw_clinical$day)
View(raw_clinical[order(raw_clinical$disease.state.ch1, raw_clinical$patient_id),])

#number of MDR patients=52
length(unique(raw_clinical$patient_id))

#Note patients 40, 68 and 33 are missing day info


unique(raw_clinical$day)
raw_clinical$time <- raw_clinical$day
raw_clinical[which(raw_clinical$day == "14"), "time"] <- "14 Days"
raw_clinical[which(raw_clinical$day == "28"), "time"] <- "28 Days"
raw_clinical[which(raw_clinical$day == "180"), "time"] <- "6 Mths"
raw_clinical[which(raw_clinical$day == "300"), "time"] <- "10 Mths"
raw_clinical[which(raw_clinical$day == "450"), "time"] <- "15 Mths"
raw_clinical[which(raw_clinical$day == "600"), "time"] <- "20 Mths"

colnames(raw_clinical)[1:3] <- c("sample_id",
                                 "disease",
                                 "full_id")

clinical <- raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))


## 2) Get groups to be compared --------------------------------------------
#MDR = Multi-drug resistant TB
#DS-TB = drug-susceptible TB
table(clinical$time)

clinical$group <- clinical$time

if(any(is.na(clinical$group))){
  clinical <- clinical[-which(is.na(clinical$group)),]
}
length(unique(clinical$patient_id)) #no patients lost from NAs


#[331] "Data has been log-transformed and normalized (percentile shift) within the Agilent GeneSpring V13 software"

# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(normalised_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

library(reshape2)
# ggplot(melt(raw_counts), aes(x = value)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Density Plot of Microarray Data",
#        x = "Expression Value",
#        y = "Density") +
#   theme_minimal()

hist(raw_counts)


## 3) GSVA and boxplot to see comparisons (genesig_D_7---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- counts_norm[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

if(all(row.names(mean_sig_zscore) != row.names(clinical))){ stop("Row names do not match", this.accession.no)}

## 3.1) Boxplot ---------------------------
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))

boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("0", "14 Days", "28 Days", "6 Mths",  "10 Mths",  "15 Mths",  "20 Mths")

boxplot_all$group <- factor(boxplot_all$group, levels = x_order)
boxplot_all$score <- as.numeric(boxplot_all$score)

# calculate p values
stat.table<- boxplot_all  %>%
  wilcox_test(score ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]
lowest_bracket <- max(boxplot_all$score) + 0.05*(max(boxplot_all$score))
stat.table$y.position <- seq(lowest_bracket, by= 0.25, length.out = nrow(stat.table))

outcome = "TB"
boxplotfig <- ggplot(boxplot_all, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_all[,1]),
  group = group)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  {if(nrow(stat.table) >0 )
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       vjust = 0.35,
                       size = 5)
  } + 
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       color = "Disease", #legend title
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy", "\n",
                        "n=", nrow(boxplot_all),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from the original study: 'Data has been log-transformed and\nnormalized (percentile shift) within the Agilent GeneSpring V13 software'","\n",
       "Curve generated using loess method (non-parametric local regression)"
)) +
  ylab (label = "Signature Score") +
  xlab (label = "Disease")

ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600, 
       height = 3600, 
       units = "px" )






## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("0", "14 Days"),
  c("0", "28 Days"),
  c("0", "6 Mths"),
  c("0", "10 Mths"),
  c("0", "15 Mths"),
  c("0", "20 Mths"),
  c("14 Days", "28 Days"),
  c("28 Days", "6 Mths")
)

#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  
  
## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1B"))

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(subset_counts)
gene_set_list <- c(signature_geneid)
#CANNOT FIND FCGR1C OR ITS ALIASES

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }


# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"



all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  

  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  

  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  

  res_current <-cbind(
    comparison = paste0(group1," vs ",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1," vs ",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
} # close pair
  
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))


# 5) ROC Curves & Forest plot-----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison)

roc_data$Comparison_plotlabel <- roc_data$Comparison

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")



# Timepoint ROC
timepoint_roc_data <- roc_data

timepoint_roc <- ggplot(timepoint_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - TB treatment timepoints",
    x = "FPR (1 - Specificity)",
    y = "TPR(Sensitivity)",
    color = "Comparison",
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"))


ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome,".png")), 
       width = 3200, 
       height = 3500, 
       units = "px")



## Forest plot -------
res_table <- forestplot_res_tablexx
  
res_table$comparison <- factor(res_table$comparison)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 


auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
          # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
        # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
        # scale_y_discrete(labels = comparison_plotlabel_levels)+
    forestplot_theme +
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1B", "\n", 
                        "Note: FCGR1CP was not annotated in GPL13497; FCGR1B was used as a proxy"),
                       size = 12, hjust = 0, x = 0)
  )
  
    
  ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

  
  
  
  
  
  
# ================================================================================== #
# 3.  OTHER DISEASES DATASETS ======================================================
# ================================================================================== #

# GSE42834  -----------------------------------------------------------------------------------------------------------------------------------

#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE42834"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

if(!exists(file.path(this.accession.res.dir, "figures"))) dir.create(file.path(this.accession.res.dir, "figures"))


## 1) Load public data -------------------------------------------------------------
#Ran in HPC
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE193777", "file=GSE193777_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

gse=getGEO(filename="GSE42834_series_matrix.txt.gz", getGPL = FALSE)

raw_metadata <- gse@phenoData@data
raw_counts <- exprs(gse)
ex <- exprs(gse)

colnames(raw_metadata) <- make.names(colnames(raw_metadata))



raw_clinical <- raw_metadata[,c("geo_accession",
                                "disease.state.ch1",
                                "gender.ch1",
                                "title")] #blood test used to diagnose tuberculosis (TB) by measuring the amount of interferon-gamma (IFN-γ) released by a person's immune cells when exposed to specific TB antigens

colnames(raw_clinical) <- c("sample_id",
                            "disease",
                            "sex",
                            "title")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))

gpl <- getGEO("GPL10558", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("ID", "Symbol")]


## 2) Get groups to be compared --------------------------------------------
table(clinical$disease)
#NOTE: all patients are household contacts, none are community controls
clinical[which(clinical$disease == "Active Sarcoid"),"group"] <- "Sarcoidosis" 
clinical[which(clinical$disease == "Active sarcoidosis"),"group"] <- "Sarcoidosis" 
clinical[which(clinical$disease == "Non-active sarcoidosis"),"group"] <- "Sarcoidosis" 
clinical[which(clinical$disease == "Sarcoid"),"group"] <- "Sarcoidosis" 

clinical[which(clinical$disease == "lung cancer"),"group"] <- "Lung cancer"

clinical[which(clinical$disease == "Pneumonia"),"group"] <- "Pneumonia"
clinical[which(clinical$disease == "pneumonia"),"group"] <- "Pneumonia"
clinical[which(clinical$disease == "Baseline"),"group"] <- "Pneumonia"

clinical[which(clinical$disease == "TB"),"group"] <- "TB"

clinical[which(clinical$disease == "Control"),"group"] <- "Control"

if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)

# Data already normalised
# background corrected, log2-transformed and 75th percentile normalised using GeneSpring 11.5):
# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(raw_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))


write.table(counts_norm, file.path("counts_norm.txt"))
write.table(clinical, file.path("clinical.txt"))



## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1C")) #Used alias FCGR1C instead of FCGR1CP

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 


## 3.1) Boxplot ---------------------------
outcome = "Lung Disease"
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

# Specify order that we want the x-axis variables in for the plot
x_order <- c("Control", "TB", "Sarcoidosis", "Pneumonia", "Lung cancer")

# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot <- boxplot_all
boxplotfig <- boxplot_func(outcome = "TB") 


boxplotfig <- boxplotfig +
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from original study: described as background \n corrected, log2-transformed and 75th percentile normalised using GeneSpring 11.5\n",
       "P values from Mann-Whitney U test shown")) 



ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
       units = "px" )




## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Control", "TB"),
  c("TB", "Lung cancer"),
  c("TB", "Pneumonia"),
  c("TB", "Sarcoidosis")
)


# Define for plotting
comparison_levels <- c(
  "Control vs TB",
  "TB vs Lung cancer",
  "TB vs Pneumonia",
  "TB vs Sarcoidosis"
)
# 
# 
# comparison_plotlabel_levels <- c(
#   "Control vs TB",
#   "TB vs Lung cancer",
#   "TB vs Pneumonia",
#   "TB vs Sarcoidosis"
# )


disease_roc_subset <- comparison_levels
timepoint_roc_subset <- NULL


disease_legend_nrow = 2
timepoint_legend_nrow = NULL

roc_plot_width = 3200
roc_plot_height = 3500



# Run roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)


roc_func2 <- function(outcome){

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical_treat[clinical_treat$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)  #default 95% CI is computed with 2000 stratified bootstrap replicates.
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  
  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1," vs ",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1," vs ",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
} # close pair
  


# saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no,"_roc_objects_", outcome, ".rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))



# 5) ROC Curves & Forest plot-----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison, levels = comparison_levels)

roc_data$Comparison_plotlabel <- roc_data$Comparison

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")


# Disease plot
disease_roc_data <- roc_data[which(roc_data$Comparison_plotlabel %in% disease_roc_subset),]

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = disease_legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = paste0("TB vs Other lung diseases", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature:  TAP1, GBP5, GBP2, FCGR1CP") 

if(nrow(disease_roc_data) >0 ){
ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")
}
#Timepoint plot
timepoint_roc_data <- roc_data[which(roc_data$Comparison_plotlabel %in% timepoint_roc_subset),]

timepoint_roc <- ggplot(timepoint_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = timepoint_legend_nrow)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = paste0("TB vs Other lung diseases", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR(Sensitivity)",
    color = "Comparison",
    caption = "Signature:TAP1, GBP5, GBP2, FCGR1CP")

if(nrow(timepoint_roc_data) > 0){
ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome, ".png")), 
       width = roc_plot_width, 
       height = roc_plot_height, 
       units = "px")
}
## Forest plot -------
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = comparison_levels)
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)


forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 
auc_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("AUC")+
    coord_cartesian( xlim=c(0.1, 1))


  sens_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=sensitivity_ci_low, 
                       xmax=sensitivity_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("Sensitivity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  spec_plot <- res_table %>% 
    ggplot(aes(y = comparison)) + 
    theme_bw() +
    geom_point(aes(x=specificity, color = comparison), shape=15, size=3) +
    geom_linerange(aes(xmin=specificity_ci_low, 
                       xmax=specificity_ci_high, 
                       color = comparison),
                   size = 1) +
    forestplot_theme +
    scale_y_discrete(labels = comparison_plotlabel_levels)+
    ylab(NULL)+
    xlab("Specificity")+
    coord_cartesian( xlim=c(0.1, 1))
  
  
  panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                            ncol = 1,
                            nrow = 3)
  
  panel_forest <- annotate_figure(
    panel_forest,
    top = text_grob(paste0(this.accession.no, " (", outcome, ")"), size = 14, hjust = 0, x = 0),
    bottom = text_grob(paste0("Senstivity and specificity calculated at Youden threshold \n",
                              "Signature:TAP1, GBP5, GBP2, FCGR1CP"), 
                       size = 12, hjust = 0, x = 0)
  )
  
  


} #close roc2 function

clinical_treat <- clinical
roc_func2(outcome = "Lung Disease")



  
# GSE19491  -----------------------------------------------------------------------------------------------------------------------------------

#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE19491"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

if(!exists(file.path(this.accession.res.dir, "figures"))) dir.create(file.path(this.accession.res.dir, "figures"))



## 1) Load public data -------------------------------------------------------------
#Ran in HPC
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE193777", "file=GSE193777_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
# gse=getGEO(filename="GSE19491_series_matrix.txt.gz", getGPL = FALSE)


gset <- getGEO("GSE19491", GSEMatrix =TRUE, getGPL=FALSE)
gse=gset$GSE19491_series_matrix.txt.gz

raw_metadata <- gse@phenoData@data
raw_counts_prelog2 <- exprs(gse)

# "The data were analyzed using BeadStudio version 1.5.1.3 using the average chip normalization method."
qx <- as.numeric(quantile(raw_counts_prelog2, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { raw_counts_prelog2[which(raw_counts_prelog2 <= 0)] <- NaN
raw_counts <- log2(raw_counts_prelog2) }


# boxplot(raw_counts_prelog2, main = "Boxplot of Average Chip Normalized Values",
#         las = 2, col = "lightblue", outline = FALSE)
# 
# 
# boxplot(raw_counts, main = "Boxplot of Average Chip Normalized Values",
#         las = 2, col = "lightblue", outline = FALSE)
# 
# hist(raw_counts)
counts_norm <- raw_counts



colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "source_name_ch1",
                                "characteristics_ch1.3",
                                "characteristics_ch1.4",
                                "title")] #blood test used to diagnose tuberculosis (TB) by measuring the amount of interferon-gamma (IFN-γ) released by a person's immune cells when exposed to specific TB antigens

colnames(raw_clinical) <- c("sample_id",
                            "classification",
                            "disease",
                            "time",
                            "title")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))

#annotations
# library("illuminaHumanv4.db") 
library("AnnotationDbi")

gpl <- getGEO("GPL6947", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("ID", "Symbol")]


## 2) Get groups to be compared --------------------------------------------
table(clinical$disease) 
#THERE ARE SEVERAL DATASETS IN HERE
#1) LONDON (USED AS TRAINING AND TESTING)
#2) SOUTH AFRICA (USED AS VALIDATION)
#3) 4 OTHER (Still, Staph, ASLE and pSLE)
#ASLE (adult systemic lupus erythematosus)
#pSLE (pediatric SLE) #

# clinical[which(clinical$disease == "healthy control: pSLE "),"group"] <- "healthy_pSLE" 

table(clinical$time)
#clean up values
clinical$time <- gsub("time post initiation of treatment: ", "", clinical$time)
clinical$disease <- gsub("illness: ", "", clinical$disease)

#get study
london_indexes <- which(grepl("LON", clinical$title))
clinical[london_indexes, "study"] <- "LON"

sa_indexes <- which(grepl("SA", clinical$title))
clinical[sa_indexes, "study"] <- "SA"

clinical$group <- ifelse(!is.na(clinical$study),
                         paste0(clinical$study, "_", clinical$disease),
                         clinical$disease)


timepoint_indexes <- which(grepl("months", clinical$time))
clinical[timepoint_indexes, "group" ] <- paste0(clinical[timepoint_indexes, "group"], 
                                                "_",
                                                clinical[timepoint_indexes, "time"])


if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)

# clinical[which(clinical$disease == "Latent"),"group"] <- "Latent TB" #from london data
# clinical[which(clinical$disease == "LATENT TB"),"group"] <- "Latent TB" #from South africa data

# 
# # Split into the different studies
# london_indexes <- which(grepl("LON", clinical$title))
# clinical_london <- clinical[london_indexes,]
# 
# #Split into the different studies
# sa_indexes <- which(grepl("SA", clinical$title))
# clinical_sa <- clinical[sa_indexes,]
# 
# #Split into the different studies
# otherdisease_indexes <- which(grepl("WB", clinical$title))
# clinical_otherdisease <- clinical[otherdisease_indexes,]


# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(raw_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))



## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP")) 

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])

#https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=503
#according to platform GPL6947, FCGR1CP is both of those probs?

"ILMN_2261600" %in% row.names(counts_norm)
"ILMN_2176063" %in%  row.names(counts_norm)

signature_geneid[4] <- "ILMN_2261600"
# signature_geneid[4] <- "ILMN_2176063"

gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 



## 3.1) Boxplot ---------------------------
outcome = "Lung Disease"
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


# Specify order that we want the x-axis variables in for the plot
x_order <- c("LON_Control (BCG+)",
             "LON_Control (BCG-)",
             "LON_Latent",
             "LON_PTB",
             "LON_Control_0_months", 
             "LON_PTB_0_months", 
             "LON_PTB_2_months", 
             "LON_PTB_12_months", 
             "SA_LATENT TB",
             "SA_PTB",
             "healthy control: Still",
             "Still",
             "healthy control: Strep and Staph",
             "Strep",
             "Staph",
             "healthy control: ASLE",
             "ASLE",
             "healthy control: pSLE",
             "PSLE")


boxplot <- boxplot_all

#get all comparisons
my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)

boxplot$group <- factor(boxplot$group, levels = x_order)
boxplot$score <- as.numeric(boxplot$score)

# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)



boxplotfig <- ggplot(boxplot, aes(
  x = factor(group, level = x_order),
  y = score,
  group = group)) +
  
  theme_bw()+
  
  boxplot_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  # stat_pvalue_manual(stat.table.gsva,
  #                    label = "p",
  #                    tip.length = 0.01,
  #                    size = 7)+
  
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18,
                                   angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))+
  
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised data was obtained from the original study (BeadStudio average chip normalisation)and log2-transformed prior to analysis")) +
  ylab (label = "Signature Score") +
  xlab (label = "Disease")


ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 5000,
       height = 3500,
       units = "px" )


LON_groups <-  c("LON_Control (BCG+)",
                 "LON_Control (BCG-)",
                 "LON_Latent",
                 "LON_PTB",
                 "LON_Control_0_months", 
                 "LON_PTB_0_months", 
                 "LON_PTB_2_months", 
                 "LON_PTB_12_months")

SA_groups <-  c("SA_LATENT TB",
                "SA_PTB")

Other_groups <- c("healthy control: Still",
                  "Still",
                  "healthy control: Strep and Staph",
                  "Strep",
                  "Staph",
                  "healthy control: ASLE",
                  "ASLE",
                  "healthy control: pSLE",
                  "PSLE")

listof_sep_plots <- list(LON = LON_groups, SA = SA_groups, Other_dx = Other_groups)

i = "SA"
## Seperate plots only ------
for (i in names(listof_sep_plots)){
  boxplot_sep_subset <-  boxplot[boxplot$group %in% listof_sep_plots[[i]],]
  
  boxplot_sep <- ggplot(
    boxplot_sep_subset, 
    aes(x = factor(group, level = listof_sep_plots[[i]]),
        y = as.numeric(boxplot_sep_subset[,1]),
        group = group)) +
    
    theme_bw()+
    
    boxplot_theme +
    
    geom_boxplot(position = position_dodge(1)) +
    
    geom_jitter(aes(color = group),
                alpha = 0.5,
                size = 2.5, 
                width = 0.3) +
    
    
    stat_summary(fun.y = mean, fill = "red",
                 geom = "point", shape = 21, size =4,
                 show.legend = TRUE) +
    

    theme(axis.text.x = element_text(size = 22
                                     # ,
                                     # angle = 90,
                                     # vjust = 0.5,
                                     # hjust=1
                                     )) +
    
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised data was obtained from the original study (BeadStudio average chip normalisation)\nand log2-transformed prior to analysis")) +
  ylab (label = "Signature Score") +
  xlab (label = "Disease")
  
  ggsave(boxplot_sep, filename = file.path(this.figure.dir, paste0("meanzscore_plot_",  i, "_",this.accession.no, ".png")),
       width = 2500,
       height = 3500,
       units = "px" )
  
}






## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  #SA TB VS Healthy
  c("SA_LATENT TB", "SA_PTB"),
  
  #SA TB VS Other Dx
  c("SA_PTB", "Still"),
  c("SA_PTB", "ASLE"),
  c("SA_PTB", "PSLE"),
  c("SA_PTB", "Strep"),
  c("SA_PTB", "Staph"),
  
  #LONDON TB VS HEALTHY
  c("LON_Control (BCG+)", "LON_PTB"),
  c("LON_Control (BCG-)", "LON_PTB"),
  c("LON_Control (BCG+)", "LON_Latent"),
  c("LON_Control (BCG-)", "LON_Latent"),
  c("LON_Latent", "LON_PTB"),
  c("LON_Control_0_months", "LON_PTB_2_months"),
  c("LON_Control_0_months", "LON_PTB_12_months"),
  
  #LONDON TB VS Other Dx
  c("LON_PTB", "Still"),
  c("LON_PTB", "ASLE"),
  c("LON_PTB", "PSLE"),
  c("LON_PTB", "Strep"),
  c("LON_PTB", "Staph"),
  
  #DISEASE
  c("healthy control: Still", "Still"),
  c("healthy control: ASLE", "ASLE"),
  c("healthy control: pSLE", "PSLE"),
  c("healthy control: Strep and Staph", "Strep"),
  c("healthy control: Strep and Staph", "Staph"))


comparison_labels <- sapply(pairwise_comparisons, function(x) {
  paste(x[1], x[2], sep = " vs ")
})

#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()




clinical_treat <- clinical
# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical_treat[clinical_treat$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
gene_set <- subset_counts[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

all(row.names(mean_sig_zscore) == row.names(subset_clinical))

  
  glm_data <- data.frame(Score = mean_sig_zscore[,"score"], Group = subset_clinical$group)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(group1, group2))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)  #default 95% CI is computed with 2000 stratified bootstrap replicates.
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  
  
  # Sensitivity confidence interval (at optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"])) 
  
  # Specificity confidence interval (at optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  
  forestplot_res_table <- rbind(forestplot_res_table, 
                                cbind(comparison = paste0(group1,"vs",group2),
                                      auc = auc(roc_obj),
                                      auc_ci_low = as.numeric(auc_ci[1]),
                                      auc_ci_high = as.numeric(auc_ci[3]),
                                      
                                      sensitivity = optimal_threshold_coords$sensitivity, 
                                      sensitivity_ci_low = ci_sens[, "2.5%"],
                                      sensitivity_ci_high = ci_sens[, "97.5%"],
                                      
                                      specificity = optimal_threshold_coords$specificity,
                                      specificity_ci_low = ci_spec[, "2.5%"],
                                      specificity_ci_high = ci_spec[, "97.5%"])
  )
  
  roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
} # close pair


# saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no,"_roc_objects_", outcome, ".rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table_", outcome, ".csv")))
write.csv(forestplot_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_mean_ztransformed_scores_forestplot_res_table_", outcome, ".csv")))





# GSE42826 -----------------------------------------------------------------------------------------------------------------------------------

#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE42826"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

if(!exists(file.path(this.accession.res.dir, "figures"))) dir.create(file.path(this.accession.res.dir, "figures"))


## 1) Load public data -------------------------------------------------------------

gset <- getGEO("GSE42826", GSEMatrix =TRUE, getGPL=FALSE)
gse=gset$GSE42826_series_matrix.txt.gz

raw_metadata <- gse@phenoData@data
raw_counts <- exprs(gse)

#Gene annotations
gpl <- getGEO("GPL10558", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("ID", "Symbol")]

# boxplot(raw_counts_prelog2, main = "Boxplot of Average Chip Normalized Values",
#         las = 2, col = "lightblue", outline = FALSE)
#
#
# boxplot(raw_counts, main = "Boxplot of Average Chip Normalized Values",
#         las = 2, col = "lightblue", outline = FALSE)
#
# hist(raw_counts)

colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "title",
                                "disease.state.ch1")]

colnames(raw_clinical) <- c("sample_id",
                            "title",
                            "disease")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))


## 2) Get groups to be compared --------------------------------------------
table(clinical$disease)
table(clinical$group)

clinical$group <- clinical$disease
clinical[which(clinical$disease == "Active Sarcoid"), "group"] <- "Active Sarcoidosis"
clinical[which(clinical$disease == "lung cancer"), "group"] <- "Lung Cancer"
clinical[which(clinical$disease == "Non-active sarcoidosis"), "group"] <- "Non-active Sarcoidosis"


if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)


#[102] "Background subtracted, normalised and log2 transformed for normalised values and non-normalised values have only had background subtraction"
# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(raw_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

library(reshape2)
ggplot(melt(raw_counts), aes(x = value)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Microarray Data",
       x = "Expression Value",
       y = "Density") +
  theme_minimal()

hist(raw_counts)




## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1C")) #Used alias FCGR1C instead of FCGR1CP

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 




## 3.1) Boxplot ---------------------------
outcome = "Lung Disease"
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

# Specify order that we want the x-axis variables in for the plot
x_order <- c("Control", "TB", "Non-active Sarcoidosis", "Active Sarcoidosis", "Pneumonia", "Lung Cancer")

# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot <- boxplot_all
boxplotfig <- boxplot_func(outcome = "Lung Disease") 


boxplotfig <- boxplotfig +
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1C", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from original study: described as background \n corrected, log2-transformed and 75th percentile normalised using GeneSpring 11.5\n",
       "P values from Mann-Whitney U test shown"))+
    scale_x_discrete(labels = function(x) ifelse(x == "Non-active Sarcoidosis", "Non-active \n Sarcoidosis", 
                                               ifelse(x == "Active Sarcoidosis", "Active \n Sarcoidosis",
                                                      x)))



ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
       units = "px" )





## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Control", "TB"),
  # c("Control", "Non-active Sarcoidosis"),
  # c("Control", "Active Sarcoidosis"),
  # c("Control", "Pneumonia"),
  # c("Control", "Lung Cancer"),
  c("TB", "Non-active Sarcoidosis"),
  c("TB", "Active Sarcoidosis"),
  c("TB", "Pneumonia"),
  c("TB", "Lung Cancer")
)



comparison_levels <- sapply(pairwise_comparisons, function(x) {
  paste(x[1], x[2], sep = " vs ")
})

comparison_plotlabel_levels <- comparison_levels


#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()
forestplot_res_table <- data.frame()




clinical_treat <- clinical

disease_roc_subset <- comparison_levels
timepoint_roc_subset <- NULL


disease_legend_nrow = 2
timepoint_legend_nrow = NULL

roc_plot_width = 3200
roc_plot_height = 3500



# Run roc function

clinical_treat <- clinical
roc_func2(outcome = "Lung Disease")



# GSE83456  ----------------------------------------------------------------------------------------------------------------------------------------

#free up space and remove objects from previous gse (except for the functions I made and gene_annot)
rm(list = setdiff(ls(), c(lsf.str())))

gc()

my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
data.dir <- file.path(main.dir, "data")
processed.dir <- file.path(data.dir, "processed")
output.dir <- file.path(main.dir, "output_v2", "public_validation")

this.accession.no <- "GSE83456"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

if(!exists(file.path(this.accession.res.dir, "figures"))) dir.create(file.path(this.accession.res.dir, "figures"))


## 1) Load public data -------------------------------------------------------------

gset <- getGEO("GSE83456", GSEMatrix =TRUE, getGPL=FALSE)
gse=gset$GSE83456_series_matrix.txt.gz

raw_metadata <- gse@phenoData@data
raw_counts <- exprs(gse)

write.csv(raw_metadata, file.path(main.dir, "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.csv")))
write.csv(raw_counts, file.path(my_main.dir, "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.csv")))

#Gene annotations
gpl <- getGEO("GPL10558", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("ID", "Symbol")]


colnames(raw_metadata) <- make.names(colnames(raw_metadata))


raw_clinical <- raw_metadata[,c("geo_accession",
                                "title",
                                "disease.state.ch1")]

colnames(raw_clinical) <- c("sample_id",
                            "title",
                            "disease")

clinical <-raw_clinical
all(row.names(raw_clinical) == colnames(raw_counts))



## 2) Get groups to be compared --------------------------------------------
table(clinical$disease)
table(clinical$group)

clinical$group <- clinical$disease
clinical[which(clinical$disease == "Sarcoid"), "group"] <- "Sarcoidosis"



if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}
table(clinical$group)


#[114] "quantile normalization, background subtraction using Genomestudio"

# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(raw_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))

library(reshape2)
ggplot(melt(raw_counts), aes(x = value)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Microarray Data",
       x = "Expression Value",
       y = "Density") +
  theme_minimal()

hist(raw_counts)



## 3) Mean of z-scored expression ---------------------------
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1C")) #Used alias FCGR1C instead of FCGR1CP

# Get the gene IDs instead of HGNCs
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
gene_set_list <- c(signature_geneid)

if(length(signature_geneid) < 4){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

mean_sig_zscore <- mean_zscore_func() #Function already defined for previous gse, same code 




## 3.1) Boxplot ---------------------------
outcome = "Lung Disease"
boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                    group = clinical$group))


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

# Specify order that we want the x-axis variables in for the plot
my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("HC", "PTB", "EPTB", "Sarcoidosis")



my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("HC", "PTB", "EPTB", "Sarcoidosis")


# Run boxplot function for all (no treatment outcome data for this gse)
# Function defined previously
this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)

boxplot <- boxplot_all
boxplotfig <- boxplot_func(outcome = "Lung Disease") 


boxplotfig <- boxplotfig +
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes\n",
       "Normalised expression data was obtained from original study: quantile normalisation and background subtracted using GenomeStudio"
       ))



ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
       units = "px" )



## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("HC", "PTB"),
  c("HC", "EPTB"),
  c("HC", "Sarcoidosis"),
  c("PTB", "EPTB"),
  c("PTB", "Sarcoidosis"),
  c("EPTB", "Sarcoidosis")
)

comparison_levels <- list(
  c("HC vs PTB"),
  c("HC vs EPTB"),
  c("HC vs Sarcoidosis"),
  c("PTB vs EPTB"),
  c("PTB vs Sarcoidosis"),
  c("EPTB vs Sarcoidosis")
)

comparison_plotlabel_levels <- comparison_levels

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

#Make roc function
this.figure.dir <- file.path(this.accession.res.dir, "figures", "roc")
if(!exists(this.figure.dir)) dir.create(this.figure.dir, recursive = TRUE)


clinical_treat <- clinical

disease_roc_subset <- comparison_levels
timepoint_roc_subset <- NULL


disease_legend_nrow = 2
timepoint_legend_nrow = NULL

roc_plot_width = 3200
roc_plot_height = 3500



# Run roc function

clinical_treat <- clinical
roc_func2(outcome = "Lung Disease")




## Differential expression TB vs Sarcoid Microarray ------------------------------------------------------------------------------------------
library(ggrepel)
diffexp.dir <- file.path(this.accession.res.dir, "diffexp")
if(!exists(diffexp.dir)) dir.create(diffexp.dir)
results.dir <- file.path(diffexp.dir, "results")
figures.dir <- file.path(diffexp.dir, "figures")

if(!exists(diffexp.dir)) dir.create(diffexp.dir)
if(!exists(results.dir)) dir.create(results.dir)
if(!exists(figures.dir)) dir.create(figures.dir)


clinical_de <- raw_metadata[,c("geo_accession",
                               "title",
                               "disease.state.ch1",
                               "gender.ch1",
                               "age.ch1",
                               "ethnicity.ch1")]

colnames(clinical_de) <- c("sample_id",
                           "title",
                           "disease",
                           "sex",
                           "age",
                           "ethnicity")

clinical_de$group <- clinical_de$disease
clinical_de[which(clinical_de$disease == "Sarcoid"), "group"] <- "Sarcoidosis"



if(any(is.na(clinical_de$disease))){
  clinical_de <- clinical_de[-which(is.na(clinical_de$disease)),]
}
table(clinical_de$group)


all(row.names(clinical_de) == colnames(counts_norm))

clinical_de$age <- as.numeric(clinical_de$age)

#Design matrix
design <- model.matrix(
  ~0 + group + age + sex,
  data = clinical_de 
)

colnames(design)[1:4] <- levels(as.factor(clinical_de$group))

#fit a linear model to each gene
fit <- lmFit(counts_norm,
             design)
# ,
# block = clinical_de$study.id, 
# correlation = corfit$consensus)

##Specify columns to compare
cont.matrix <- makeContrasts(
  test1 = PTB - Sarcoidosis,
  test2 = EPTB - Sarcoidosis,
  test3 = PTB - HC,
  test4 = EPTB - HC,
  test5 = Sarcoidosis - HC,
  levels = design
  
) 


nameconvert <- as.data.frame(cbind(name = colnames(cont.matrix),
                                   contrast = c(
                                     "PTB - Sarcoidosis",
                                     "EPTB - Sarcoidosis",
                                     "PTB - HC",
                                     "EPTB - HC",
                                     "Sarcoidosis - HC"
                                   )))
#Rename certain genes to their aliases
gene_annot$Symbol_alias <- gene_annot$Symbol
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "FCGR1C")] <- "FCGR1CP"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "Septin 4")] <- "SEPTIN4"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "FAM26F")] <- "CALHM6"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "TRMT2A")] <- "TRAMT2A"


listoftT <- list()
listoftT2 <- list()
listofvolcano <- list()

i = "test1"
for (i in colnames(cont.matrix)){
  #contrasts.fit converts the coefficients and standard errors to reflect the contrasts rather than the original design matrix, but does not compute t-statistics or p-values. 
  fit2 <- contrasts.fit(fit, contrast=cont.matrix[,i])
  
  #eBayes computes t-statistics and p-values from the coefficients and standard errors.
  fit2 <- eBayes(fit2)
  
  #topTableF ranks genes on the basis of moderated F-statistics for one or more coefficients.
  tT <- topTable(fit2, adjust="BH", sort.by="P", number=nrow(fit2))
  
  tT$legend <-  ifelse(
    tT$adj.P.Val< 0.05 & tT$logFC > 1, "Upregulated",
    ifelse(
      tT$adj.P.Val < 0.05 & tT$logFC < -1, "Downregulated",
      "Not Sig"))
  
  tT$gene <- gene_annot[match(row.names(tT), gene_annot$ID), "ID"]
  
  selection <-which((tT$logFC>1|tT$logFC< -1)&tT$adj.P.Val<0.05)
  tT2 <- tT[selection,]
  
  listoftT[[i]] <- tT
  listoftT2[[i]] <- tT2
  
  write.csv(tT, file.path(results.dir, paste0("tT_", names(listoftT[i]), ".csv")))
  write.csv(tT2, file.path(results.dir, paste0("tT2_", names(listoftT[i]), ".csv")))
  
  
  volcano <- ggplot(tT, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = legend)) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not Sig" = "grey", "Upregulated" = "red"))+
    geom_hline(yintercept =-log10(max(as.data.frame(tT2[,"P.Value"]))),colour="black", linetype="dashed") +
    geom_vline(xintercept =-1,colour="black", linetype="dashed")+
    geom_vline(xintercept =1,colour="black", linetype="dashed")+
    geom_text_repel(data = subset(tT2[1:20,]),
                    aes(label= gene),
                    size = 3, 
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines") 
    ) +
    theme_bw(base_size = 12) +
    
     theme(legend.position = "bottom") +
    
    labs(title = nameconvert[which(nameconvert$name == i), "contrast"])
  
  listofvolcano[[i]] <- volcano
  
  ggsave(volcano, 
         filename = file.path(figures.dir, 
                              paste0("volcano_",
                                     nameconvert[which(nameconvert$name == i), "contrast"],
                                     ".png")),
         width = 16,
         height = 16,
         units = "cm")
}




