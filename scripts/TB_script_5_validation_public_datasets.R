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

setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))
# 
# shiny.data.dir <- file.path(my_directory,"TB", "shiny", "data", "public", this.accession.no)
# if(!exists(shiny.data.dir)) dir.create(shiny.data.dir, recursive = TRUE)
# write.table(counts_vst, file.path(shiny.data.dir, "counts_vst.txt"))
# write.table(clinical, file.path(shiny.data.dir, "clinical.txt"))

## 3) Mean of z-scored expression and boxplot to see comparisons ---------------------------


# Get the gene IDs instead of HGNCs
mean_zscore_func <- function(){
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- c(signature_geneid)

# For each sample, get the mean standardized expression of genes in the 4-gene signature

gene_set <- counts_vst[gene_set_list,]
gene_set <- t(gene_set) #genes are columns so we can z-score column-wise (centre = centre each gene around its own mean)
      # Gene-wise z-score = for each gene, how does this sample compare to all other samples for that gene, then average across our 7 genes
      # This tells you if a sample has collectively high expression of your gene set relative to the cohort      
gene_set_zscore<-scale(gene_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
      #take the mean of the scaled+centred data for each gene and use this as each sample's score
      #each row is a sample and each column is a gene. we are taking the average ACROSS the columns so every sample gets 1 score
mean_sig_zscore<- data.frame(rowMeans(gene_set_zscore))
colnames(mean_sig_zscore) <- "score"

if(all(row.names(mean_sig_zscore) != row.names(clinical))){ stop()}
return(mean_sig_zscore)


} #close function

mean_sig_zscore <- mean_zscore_func()

boxplot_all <- as.data.frame(cbind(mean_zscore = mean_sig_zscore,
                                        disease = clinical$disease,
                                    group = clinical$group,
                                    treatment_outcome = clinical$treatmentresult))



boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

# Make boxplot function
boxplot_func <- function(outcome){
    
  

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures", "boxplot")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

  
my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)


x_order <- c("Healthy", "Lungdx_ctrl", "MTP_ctrl", "TB_DX", "TB_day_7", "TB_week_4", "TB_week_24")


boxplot$group <- factor(boxplot$group, levels = x_order)
boxplot$score <- as.numeric(boxplot$score)


stat.table<- boxplot  %>%
  wilcox_test(score ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]
lowest_bracket <- max(boxplot$score) + 0.05*(max(boxplot$score))
stat.table$y.position <- seq(lowest_bracket, by= 0.2, length.out = nrow(stat.table))


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
  labs(title = paste0("Signature Analysis: ", this.accession.no, " (", outcome, ")"),
       caption = paste0("Signature:TAP1, GBP5, GBP2, FCGR1CP", "\n", "n=", nrow(boxplot),"\n",
       "Signature scores calculated as mean of z-scored expression of signature genes \n",
       "Raw counts were VST normalised (DESeq2 1.42.1)",
       "P values from Mann-Whitney U test shown")) +
  ylab (label = "Signature Score") +
  xlab (label = "Disease")


ggsave(boxplotfig, filename = file.path(this.figure.dir, paste0("meanzscore_plot_", this.accession.no, "_", outcome, ".png")),
       width = 3500,
       height = 3600,
       units = "px" )
} #end boxplot function

#Run boxplot function for all, cured and non-cured

for (outcome in c("all_outcomes", "not_cured", "cured")){
    if(outcome == "all_outcomes"){
  boxplot <- boxplot_all
  }
  
  if(outcome == "not_cured"){
  boxplot <- boxplot_all[which(boxplot_all$treatment_outcome == "Not Cured" | boxplot_all$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    boxplot <- boxplot_all[which(boxplot_all$treatment_outcome == "Definite Cure" | boxplot_all$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
}
boxplot_func(outcome = "all_outcomes")
boxplot_func(outcome = "not_cured")
boxplot_func(outcome = "cured")


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

roc_func <- function(outcome){
  
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
  subset_clinical <- clinical_treat[clinical_treat$group %in% c(group1,group2),]
  subset_counts <- counts_vst[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # Mean-zscore
# Get the gene IDs instead of HGNCs
gene_set_list <- list(c("TAP1","GBP5","GBP2","FCGR1CP"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
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



## 5) ROC Curves -----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(roc_objects), function(comparison) {
  data.frame(
    TPR = rev(roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison, levels = c(
  "HealthyvsTB_DX",
  "MTP_ctrlvsTB_DX",
  "TB_DXvsTB_day_7",
  "TB_DXvsTB_week_4",
  "TB_DXvsTB_week_24",
  "HealthyvsTB_week_24",
  "MTP_ctrlvsTB_week_24"
))

roc_data$Comparison_plotlabel <- roc_data$Comparison

levels(roc_data$Comparison_plotlabel) <-  c(
  "Healthy vs TB_T0",
  "MTP Controls vs TB_T0",
  "TB_T0 vs TB_Day7",
  "TB_T0 vs TB_Wk4",
  "TB_T0 vs TB_Wk24",
  "Healthy vs TB_Wk24",
  "MTP Controls vs TB_Wk24"
)  

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data[which(roc_data$Comparison == "HealthyvsTB_DX" | 
                                     roc_data$Comparison == "MTP_ctrlvsTB_DX" |
                                     roc_data$Comparison == "HealthyvsTB_week_24" |
                                     roc_data$Comparison == "MTP_ctrlvsTB_week_24" ),]

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = paste0("ROC - TB treatment timepoints", " (", outcome, ")"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature:  TAP1, GBP5, GBP2, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, paste0("disease_roc_", outcome,".png")), 
       width = 3000, 
       height = 3200, 
       units = "px")

# Timepoint ROC
timepoint_roc_data <- roc_data[which(roc_data$Comparison == "TB_DXvsTB_day_7" | 
                                       roc_data$Comparison == "TB_DXvsTB_week_4" |
                                       roc_data$Comparison == "TB_DXvsTB_week_24" ), ]

timepoint_roc <- ggplot(timepoint_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom",
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


ggsave(timepoint_roc, filename = file.path(this.figure.dir, paste0("timepoint_roc_", outcome, ".png")), 
       width = 3000, 
       height = 3200, 
       units = "px")


## Forest plot -------
res_table <- forestplot_res_table
  
res_table$comparison <- factor(res_table$comparison, levels = levels(roc_data$Comparison))
  
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
    # scale_colour_manual(values = colours)+
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
    # scale_colour_manual(values = colours)+
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
    # scale_colour_manual(values = colours)+
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
  
  
  ggsave(panel_forest, filename= file.path(this.figure.dir, paste0("forestplot_panel_", outcome,".png")),
         width = 15, height = 20, units = "cm",   bg = "white"  )

} #close roc function


for (outcome in c("all_outcomes", "non_cured","cured")){
    if(outcome == "all_outcomes"){
  clinical_treat <- clinical
  }
  
  if(outcome == "not_cured"){
  clinical_treat <- clinical[which(clinical$treatmentresult == "Not Cured" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
  
  if( outcome == "cured"){
    clinical_treat <- clinical[which(clinical$treatmentresult == "Definite Cure" | clinical$disease %in% c("Healthy", "Lungdx_ctrl", "MTP_ctrl")) ,]
  }
    roc_func(outcome = outcome)
    
 
}
