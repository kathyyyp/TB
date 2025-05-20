
# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

#Mac
# my_directory <- "/Users/kathyphung/Library/CloudStorage/OneDrive-UTS/Documents/RBMB"
# my_directory <- "/Volumes/One Touch/RBMB"

# Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")

# #Linux
# my_directory <- "/shared/homes/165861/RBMB"
# setwd(file.path(my_directory,"TB"))
# .libPaths("~/RBMB/rlibrary")

setwd(file.path(my_directory,"TB"))

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

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #

output.dir <- file.path(my_directory,"TB","output", "public_validation")
if(!exists(output.dir)) dir.create(output.dir, recursive = TRUE)

figures.dir <- file.path(output.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)


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
                                "characteristics_ch1.6")]

colnames(raw_clinical) <- c("sample_id",
                            "disease",
                            "subject",
                            "treatmentresult",
                            "time",
                            "timetonegativity")

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


## 3) GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_vst), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 
  



my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("Healthy", "Lungdx_ctrl", "MTP_ctrl", "TB_DX", "TB_day_7", "TB_week_4", "TB_week_24")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))



boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)+
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  theme(axis.text.x = element_text(size = 15))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")), 
       width = 3500, 
       height = 3200, 
       units = "px" )


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

# Create a list to store AUC values and roc objects
GSE89403_res_table <- data.frame()
GSE89403_roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_vst[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  GSE89403_res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  GSE89403_res_table <- rbind(GSE89403_res_table, GSE89403_res_current)
  
  GSE89403_roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj

}
  
 saveRDS(GSE89403_roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
 write.csv(GSE89403_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))
  


## 5) ROC Curves -----------------------------------------------------------

 # Convert ROC data to a format suitable for ggplot
 roc_data <- do.call(rbind, lapply(names(GSE89403_roc_objects), function(comparison) {
   data.frame(
     TPR = rev(GSE89403_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
     FPR = rev(1 - GSE89403_roc_objects[[comparison]]$specificities),  # False Positive Rate
     Comparison = comparison,
     auc = rev(GSE89403_roc_objects[[comparison]]$auc)
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
 
 roc_data$ci <- GSE89403_res_table[match(roc_data$Comparison, GSE89403_res_table$comparison), "ci"]

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
     title = "ROC - Control vs TB",
     x = "FPR (1 - Specificity)",
     y = "TPR (Sensitivity)",
     color = "Comparison",
     caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 
 
 ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc.png"), 
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
     title = "ROC - TB treatment timepoints",
     x = "FPR (1 - Specificity)",
     y = "TPR(Sensitivity)",
     color = "Comparison",
     caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP")
 
 
 ggsave(timepoint_roc, filename = file.path(this.figure.dir, "timepoint_roc.png"), 
        width = 3000, 
        height = 3200, 
        units = "px")
 
 
 
 

# GSE193777 -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE193777"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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




## 3) GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_vst), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)

x_order <- c("Healthy", "Icp_TB_bl", "Icp_TB_fu", "Sub_TB_bl", "Sub_TB_fu", "Active_TB")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 7)+
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  scale_x_discrete(labels= c("Healthy" = "Healthy",
                             "Active_TB" = "Active\n TB", 
                             "Icp_TB_bl" = "Incipient TB \n baseline",
                             "Icp_TB_fu" = "Incipient TB \n followup",
                             "Sub_TB_bl" = "Subclinical TB \n baseline",
                             "Sub_TB_fu" = "Subclinical TB \n followup"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_vst_", this.accession.no, ".png")), 
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

# Create a list to store AUC values and roc objects
GSE193777_res_table <- data.frame()
GSE193777_roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_vst[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  GSE193777_res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  GSE193777_res_table <- rbind(GSE193777_res_table, GSE193777_res_current)
  
  GSE193777_roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
}

saveRDS(GSE193777_roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(GSE193777_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))



## 5) ROC Curves -----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(GSE193777_roc_objects), function(comparison) {
  data.frame(
    TPR = rev(GSE193777_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - GSE193777_roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(GSE193777_roc_objects[[comparison]]$auc)
  )
}))


roc_data$Comparison <- factor(roc_data$Comparison, levels = c(
  "HealthyvsActive_TB",
  "HealthyvsIcp_TB_bl",
  "HealthyvsIcp_TB_fu",
  "HealthyvsSub_TB_bl",
  "HealthyvsSub_TB_fu",
  "Active_TBvsIcp_TB_bl",
  "Active_TBvsIcp_TB_fu",
  "Active_TBvsSub_TB_bl",
  "Active_TBvsSub_TB_fu"
))

roc_data$Comparison_plotlabel <- roc_data$Comparison

levels(roc_data$Comparison_plotlabel) <-  c(
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

roc_data$ci <- GSE193777_res_table[match(roc_data$Comparison, GSE193777_res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data[which(roc_data$Comparison_plotlabel == "Healthy vs Active TB" | 
                                     roc_data$Comparison_plotlabel == "Healthy vs Incipient TB baseline" |
                                     roc_data$Comparison_plotlabel == "Healthy vs Incipient TB followup" |
                                     roc_data$Comparison_plotlabel == "Healthy vs Subclinical TB baseline"|
                                     roc_data$Comparison_plotlabel == "Healthy vs Subclinical TB followup"),]

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - Control vs TB",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")


# Timepoint ROC
timepoint_roc_data <- roc_data[which(roc_data$Comparison_plotlabel == "Active TB vs Incipient TB baseline" | 
                                       roc_data$Comparison_plotlabel == "Active TB vs Incipient TB followup" |
                                       roc_data$Comparison_plotlabel == "Active TB vs Subclinical TB baseline" |
                                       roc_data$Comparison_plotlabel == "Active TB vs Subclinical TB followup"), ]

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
    title = "ROC - TB treatment timepoints",
    x = "FPR (1 - Specificity)",
    y = "TPR(Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP")


ggsave(timepoint_roc, filename = file.path(this.figure.dir, "timepoint_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")




# GSE79362 -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE79362"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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



## 3) GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_vst), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("Healthy", "Active TB")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group")

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))



boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 8)+
  
  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")), 
       width = 3500, 
       height = 3200, 
       units = "px" )


## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Healthy", "Active TB")
)

# Create a list to store AUC values and roc objects
GSE79362_res_table <- data.frame()
GSE79362_roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_vst[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_voom$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  GSE79362_res_current <-cbind(
    comparison = paste0(group1,"vs",group2),
    samples_group1 = paste(group1, "=", sum(glm_data$Group == group1)),
    samples_group2 = paste(group2, "=", sum(glm_data$Group == group2)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  GSE79362_res_table <- rbind(GSE79362_res_table, GSE79362_res_current)
  
  GSE79362_roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
}

saveRDS(GSE79362_roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(GSE79362_res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))




## 5) ROC Curves -----------------------------------------------------------

# Convert ROC data to a format suitable for ggplot
roc_data <- do.call(rbind, lapply(names(GSE79362_roc_objects), function(comparison) {
  data.frame(
    TPR = rev(GSE79362_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - GSE79362_roc_objects[[comparison]]$specificities),  # False Positive Rate
    Comparison = comparison,
    auc = rev(GSE79362_roc_objects[[comparison]]$auc)
  )
}))

roc_data$Comparison <- factor(roc_data$Comparison)

roc_data$Comparison_plotlabel <- roc_data$Comparison

levels(roc_data$Comparison_plotlabel) <-  c("Healthy vs Active TB")  

roc_data$ci <- GSE79362_res_table[match(roc_data$Comparison,GSE79362_res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - Control vs TB",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc.png"), 
       width = 3000, 
       height = 3200, 
       units = "px")




# ================================================================================== #
# 2.  OTHER DISEASES DATASETS ======================================================
# ================================================================================== #

# GSE42834  -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE42834"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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

#annotations
# library("illuminaHumanv4.db") 
# library("AnnotationDbi")

gpl <- getGEO("GPL10558", destdir = ".")
gene_annot <- Table(gpl)  # Convert to a data frame
gene_annot <- gene_annot[,c("ID", "Symbol")]
# 
# my_probe_ids <- row.names(counts_norm)
# gene_annot <- data.frame(ProbeID = my_probe_ids)
# 
# gene_annot <- gene_annot %>%
#   mutate(GeneSymbol = mapIds(illuminaHumanv4.db, keys = ProbeID, 
#                              column = "SYMBOL", keytype = "PROBEID", 
#                              multiVals = "first"))


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

# counts_vst <- vst(as.matrix(raw_counts))
counts_norm <- as.matrix(raw_counts)

counts_norm <- counts_norm[,row.names(clinical)]
all(row.names(raw_clinical) == colnames(raw_counts))




## 3) GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1C"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)

x_order <- c("Control", "TB", "Sarcoidosis", "Pneumonia", "Lung cancer")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 7)+

  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")), 
       width = 3600, 
       height = 3200, 
       units = "px" )








## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest
pairwise_comparisons <- list(
  c("Control", "TB"),
  c("TB", "Sarcoidosis"),
  c("TB", "Pneumonia"),
  c("TB", "Lung cancer")
)

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  roc_objects[[paste0(group1,"vs",group2)]] <- roc_obj
  
}

saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))



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


roc_data$Comparison <- factor(roc_data$Comparison)

roc_data$Comparison_plotlabel <- roc_data$Comparison

levels(roc_data$Comparison_plotlabel) <-  c(
  "Control vs TB",
  "TB vs Lung cancer",
  "TB vs Pneumonia",
  "TB vs Sarcoidosis")

roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison_plotlabel,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data

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
    title = "ROC - TB vs Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")



# GSE19491  -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE19491"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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
library("illuminaHumanv4.db") 
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




## 3) ALL GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"]) 

#https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=503
#according to platform GPL6947, FCGR1CP is both of those probs?

"ILMN_2261600" %in% row.names(counts_norm)
"ILMN_2176063" %in%  row.names(counts_norm)

signature_geneid[7] <- "ILMN_2261600"
# signature_geneid[7] <- "ILMN_2176063"

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)

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

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)
# 
# stat.table.gsva <- boxplot_gsva  %>%
#   wilcox_test(V1 ~ group,
#               paired = FALSE) %>%
#   add_xy_position(x = "group") 
# 
# stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
# lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
# stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
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
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_all_studies_", this.accession.no, ".png")), 
       width = 5000, 
       height = 3200, 
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
boxplot_sep_subset <-  boxplot_gsva[boxplot_gsva$group %in% listof_sep_plots[[i]],]

boxplot_sep <- ggplot(
  boxplot_sep_subset, 
  aes(x = factor(group, level = listof_sep_plots[[i]]),
  y = as.numeric(boxplot_sep_subset[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  

  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 22
                                   #,
                                   # angle = 90, 
                                   # vjust = 0.5,
                                   # hjust=1
                                   ))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")


ggsave(boxplot_sep, filename = file.path(this.figure.dir, paste0("gsva_plot_", i , this.accession.no, ".png")), 
       width = 2500, 
       height = 3200, 
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

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"]) 
  signature_geneid[7] <- "ILMN_2261600"
  # signature_geneid[7] <- "ILMN_2176063"
  
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
}

saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))



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


roc_data$Comparison <- factor(roc_data$Comparison)


roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")



### 1) South Africa TB VS Healthy  ---------------------------------
levels(roc_data$Comparison)

SA_roc_tb_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[18])

SA_roc_tb <- ggplot(SA_roc_tb_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - South Africa Latent TB vs TB",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(SA_roc_tb, filename = file.path(this.figure.dir, "SA_roc_tb.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")

### 2) South Africa TB VS Other lung disease   ---------------------------------
levels(roc_data$Comparison)

SA_roc_lungdx_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[19:23])

SA_roc_lungdx <- ggplot(SA_roc_lungdx_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - South Africa TB vs Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(SA_roc_lungdx, filename = file.path(this.figure.dir, "SA_roc_TBvsLungdx.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")



### 3) London TB VS Healthy   ---------------------------------
levels(roc_data$Comparison)

LON_roc_tb_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[c(6:9,12)])

LON_roc_tb <- ggplot(LON_roc_tb_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - London TB vs Healthy",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(LON_roc_tb, filename = file.path(this.figure.dir, "LON_roc_tb.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")


### 4) London across timepoints   ---------------------------------
levels(roc_data$Comparison)

LON_roc_tb_time_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[c(10:11)])

LON_roc_tb_time <- ggplot(LON_roc_tb_time_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - London TB across treatment timepoints",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(LON_roc_tb_time, filename = file.path(this.figure.dir, "LON_roc_tb_time.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")


### 5) London TB VS Other lung disease   ---------------------------------
levels(roc_data$Comparison)

LON_roc_lungdx_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[13:17])

LON_roc_lungdx <- ggplot(LON_roc_lungdx_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - London TB vs Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(LON_roc_lungdx, filename = file.path(this.figure.dir, "LON_roc_TBvsLungdx.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")


### 6) Other lung disease   ---------------------------------
levels(roc_data$Comparison)

roc_lungdx_data <- roc_data %>%
  filter(Comparison %in% levels(roc_data$Comparison)[1:5])

roc_lungdx <- ggplot(roc_lungdx_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(roc_lungdx, filename = file.path(this.figure.dir, "roc_lungdx.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")






# # GSE40738 MICRORNA -----------------------------------------------------------------------------------------------------------------------------------
# this.accession.no <- "GSE40738"
# setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))
# 
# ## 1) Load public data -------------------------------------------------------------
# #Ran in HPC
# # urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# # path <- paste(urld, "acc=GSE193777", "file=GSE193777_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
# # tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
# # gse=getGEO(filename="GSE19491_series_matrix.txt.gz", getGPL = FALSE)
# 
# 
# gset <- getGEO("GSE40738", GSEMatrix =TRUE, getGPL=FALSE)
# gse=gset$GSE40738_series_matrix.txt.gz
# 
# #  The 'normexp' method with offset=10 was used for background-correction, and followed with within-array global 'loess' normalization, with span=1/3 and values from probe-spots with flag-values >1 ignored, and then, between-array 'quantile' normalization. For probe-set summarization, the median of the values from multiple spots was used when the ratio of the maximum and minimum values was >1.5; else, the mean was used."
# raw_metadata <- gse@phenoData@data
# raw_counts <- exprs(gse)
# gene_annot <- read.table(file.path(my_directory,"TB", "data", "public", "GSE89403", "gene_annot.txt"))
# 
# 
# # boxplot(raw_counts_prelog2, main = "Boxplot of Average Chip Normalized Values",
# #         las = 2, col = "lightblue", outline = FALSE)
# # 
# # 
# # boxplot(raw_counts, main = "Boxplot of Average Chip Normalized Values",
# #         las = 2, col = "lightblue", outline = FALSE)
# # 
# # hist(raw_counts)
# counts_norm <- raw_counts
# 
# 
# 
# colnames(raw_metadata) <- make.names(colnames(raw_metadata))
# 
# 
# raw_clinical <- raw_metadata[,c("geo_accession",
#                                 "title",
#                                 "disease.state.ch1")] 
# 
# colnames(raw_clinical) <- c("sample_id",
#                             "title",
#                             "disease")
# 
# clinical <-raw_clinical
# all(row.names(raw_clinical) == colnames(raw_counts))
# 
# 
# ## 2) Get groups to be compared --------------------------------------------
# table(clinical$disease)
# table(clinical$group)
# 
# #NOTE: all patients are household contacts, none are community controls
# 
# clinical$group <- clinical$disease
# 
# clinical[which(clinical$disease == "lung without lung cancer or non-cancerous nodule"), "group"] <- "Control"
# 
# clinical[which(clinical$disease == "granuloma of lung" |
#                                    clinical$disease == "hamartoma of lung"), "group"] <- "Benign"
# 
# clinical[which(clinical$disease == "adenocarcinoma of lung"), "group"] <- "Adenocarcinoma"
# clinical[which(clinical$disease == "squamous cell carcinoma of lung"), "group"] <- "Squamous Cell Carcinoma"
# 
# clinical$group[!clinical$group %in% c("Control", "Benign", "Adenocarcinoma", "Squamous Cell Carcinoma")] <- "Other Cancer"
# #58 controls had a high risk for developing lung cancer because of age >60 years (n = 40) and/or a history of cigarette smoking
# # (n = 57)
# 
# if(any(is.na(clinical$disease))){
#   clinical <- clinical[-which(is.na(clinical$disease)),]
# }
# table(clinical$group)
# 
# # counts_vst <- vst(as.matrix(raw_counts))
# counts_norm <- as.matrix(raw_counts)
# 
# counts_norm <- counts_norm[,row.names(clinical)]
# all(row.names(raw_clinical) == colnames(raw_counts))
# 
# 
# 
# 
# ## 3) GSVA and boxplot to see comparisons ---------------------------
# gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
# 
# signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "GeneID"])
# gene_set_list <- list(c(signature_geneid))
# 
# gsvapar <- gsvaParam(as.matrix(counts_norm), 
#                      gene_set_list, 
#                      maxDiff = TRUE, 
#                      minSize = 1)
# 
# gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway
# 
# 
# all(row.names(gsva_res) == row.names(clinical))
# 
# boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
#                                     group = clinical$group))





# GSE42826 -----------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE42826"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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

## 3) GSVA and boxplot to see comparisons ---------------------------
gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1C"))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
if(length(signature_geneid) < 6){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))


gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)

x_order <- c("Control", "TB", "Non-active Sarcoidosis", "Active Sarcoidosis", "Pneumonia", "Lung Cancer")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 7)+
  
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  scale_x_discrete(labels = function(x) ifelse(x == "Non-active Sarcoidosis", "Non-active \n Sarcoidosis", 
                                               ifelse(x == "Active Sarcoidosis", "Active \n Sarcoidosis",
                                                      x)))+

  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")), 
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

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
}

saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))



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


roc_data$Comparison <- factor(roc_data$Comparison)


roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data

disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 4)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - TB vs Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc2.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")





# GSE83456  ----------------------------------------------------------------------------------------------------------------------------------------

this.accession.no <- "GSE83456"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

## 1) Load public data -------------------------------------------------------------

gset <- getGEO("GSE83456", GSEMatrix =TRUE, getGPL=FALSE)
gse=gset$GSE83456_series_matrix.txt.gz

raw_metadata <- gse@phenoData@data
raw_counts <- exprs(gse)

write.csv(raw_metadata, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_metadata.csv")))
write.csv(raw_counts, file.path(my_directory,"TB", "data", "public", this.accession.no,  paste0(this.accession.no, "_raw_counts.csv")))

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

## 3) GSVA and boxplot to see comparisons (genesig_D_7)---------------------------

gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1C"))
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
if(length(signature_geneid) < 6){
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))


gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("HC", "PTB", "EPTB", "Sarcoidosis")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  {if(nrow(stat.table.gsva) >0 )
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 7)
  } + 
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP"
       # caption = "Signature:S100A8"
       # caption = paste0("Signature: ", paste0(gene_set_list[[1]], collapse = " "))

       
       ) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
       units = "px" )



## 3) GSVA and boxplot to see comparisons (One-by-one) ---------------------------

tT_TBT0vsHCT0 <- read.csv(file.path(my_directory, "TB", "results", "tT", "contrast1.csv"), row.names = 1)
all_genes <- row.names(tT_TBT0vsHCT0)

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)
if(!exists(file.path(this.figure.dir, "boxplots_all_46_genes"))) dir.create(file.path(this.figure.dir, "boxplots_all_46_genes"))


genesig_D_7 <- c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1C")

all_genes[!all_genes %in% gene_annot$Symbol]

gene_annot$Symbol_alias <- gene_annot$Symbol
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "FCGR1C")] <- "FCGR1CP"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "Septin 4")] <- "SEPTIN4"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "FAM26F")] <- "CALHM6"
gene_annot$Symbol_alias[which(gene_annot$Symbol_alias == "TRMT2A")] <- "TRAMT2A"

for (i in all_genes){ # for all 46 genes
# for (i in genesig_D_7){  # for 7 genes 
    print(i)

  gene <- gene_annot[which(gene_annot$Symbol_alias == i), "Symbol"]

gene_set_list <- list(c(gene))

signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])

if(length(signature_geneid) < 1){
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))


gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("HC", "PTB", "EPTB", "Sarcoidosis")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  {if(nrow(stat.table.gsva) >0 )
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 7)
  } + 
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       # caption = paste0("Signature: ", paste0(i))
       caption = paste0("Signature: ", paste(i, tT_TBT0vsHCT0[i, "legend"]))
       ) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")


# genesig_D_7
# ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, "_", i,".png")), 
#        width = 3600, 
#        height = 3200, 
#        units = "px" )

# all 46 genes
ggsave(boxplotfinal2, filename = file.path(this.figure.dir, "boxplots_all_46_genes", paste0("gsva_plot_", this.accession.no, "_", i,".png")), 
       width = 3600, 
       height = 3200, 
       units = "px" )
}





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

# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$Symbol), "ID"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
}

saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))



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


roc_data$Comparison <- factor(roc_data$Comparison)


roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")




# Disease ROC 
disease_roc_data <- roc_data

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
    title = "ROC - TB vs Other lung diseases",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(this.figure.dir, "disease_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")

## Differential expression TB vs Sarcoid Microarray ------------------------------------------------------------------------------------------

diffexp.dir <- file.path(this.accession.res.dir, "diffexp")
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

tT$gene <- gene_annot[match(row.names(tT), gene_annot$ID), "Symbol_alias"]

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




## Heatmap --------------------------------------------------------------------------------------------------------------------------










# GSE147690 (German Identification cohort) ----------------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE147690"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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

# ## 1) Load public data -------------------------------------------------------------
this.accession.no <- "GSE147690"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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


## 2) Get groups to be compared --------------------------------------------
#MDR = Multi-drug resistant TB
#DS-TB = drug-susceptible TB
table(clinical$disease)
table(clinical$group)

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

## 3) GSVA and boxplot to see comparisons (genesig_D_7---------------------------

gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1B")) #NR_027484 is FCGR1C on GPL13497, but has main name FCGR1B
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
#CANNOT FIND FCGR1C OR ITS ALIASES


if(length(signature_geneid) < 6){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway

# gsva_res_normdata <- gsva_res


all(row.names(gsva_res) == row.names(clinical))


boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    disease = clinical$disease,
                                    day = clinical$day,
                                    months = clinical$months))
boxplot_gsva[which(boxplot_gsva$disease == "Healthy Controls"), "months"] <- 0

gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.text = element_text(size = 18)) 

# my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)
# 
# x_order <- c("Control", "TB", "Non-active Sarcoidosis", "Active Sarcoidosis", "Pneumonia", "Lung Cancer")
# 
# boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
# boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)
# 
# stat.table.gsva <- boxplot_gsva  %>%
#   wilcox_test(V1 ~ group,
#               paired = FALSE) %>%
#   add_xy_position(x = "group") 
# 
# stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
# lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
# stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))
# 

#geom_point, split by disease
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.numeric(months),
  y = as.numeric(boxplot_gsva[,1]),
  group = as.factor(disease))) +
  
  theme_bw()+
    gsva_theme +
    geom_point(aes(colour=disease)) +
  geom_smooth(method = "lm", se = FALSE, aes(colour=disease))+
  
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
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       color = "Disease", #legend title
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1B") +
  scale_color_manual(values = c("Healthy Controls" = "#7CAE00",
                     "non MDR" = "#00BFC4",
                     "MDR" = "#F8766D"))+
  ylab (label = "Enrichment Score") +
  xlab (label = "Months of Treatment")


library("ggfortify")
### pca plot -----
# pca_res <- prcomp(t(as.data.frame(counts_norm)), scale. = TRUE, center = TRUE) #center = TRUE
# autoplot(pca_res,data=clinical, colour = "months")

#### heatmap -----



this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, "lm.png")), 
       width = 3600, 
       height = 3200, 
       units = "px" )








# GSE147689 (German validation cohort) ----------------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE147689"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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
this.accession.no <- "GSE147689"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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

library(reshape2)
# ggplot(melt(raw_counts), aes(x = value)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Density Plot of Microarray Data",
#        x = "Expression Value",
#        y = "Density") +
#   theme_minimal()

hist(raw_counts)

## 3) GSVA and boxplot to see comparisons (genesig_D_7---------------------------

gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1B")) #NR_027484 is FCGR1C on GPL13497, but has main name FCGR1B
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
#CANNOT FIND FCGR1C OR ITS ALIASES



if(length(signature_geneid) < 6){ 
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway

# gsva_res_normdata <- gsva_res


all(row.names(gsva_res) == row.names(clinical))


boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    disease = clinical$disease,
                                    day = clinical$day,
                                    months = clinical$months))

gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.text = element_text(size = 18)) 

# my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)
# 
# x_order <- c("Control", "TB", "Non-active Sarcoidosis", "Active Sarcoidosis", "Pneumonia", "Lung Cancer")
# 
# boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
# boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)
# 
# stat.table.gsva <- boxplot_gsva  %>%
#   wilcox_test(V1 ~ group,
#               paired = FALSE) %>%
#   add_xy_position(x = "group") 
# 
# stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
# lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
# stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))
# 

#geom_point, split by disease
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.numeric(months),
  y = as.numeric(boxplot_gsva[,1]),
  group = as.factor(disease))) +
  
  theme_bw()+
  gsva_theme +
  geom_point(aes(colour=disease)) +
  geom_smooth(method = "loess", se = FALSE, aes(colour=disease))+
  
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
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       color = "Disease", #legend title
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1B") +
  ylab (label = "Enrichment Score") +
  xlab (label = "Months of Treatment")

 


this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")), 
       width = 3600, 
       height = 3200, 
       units = "px" )






# Study visits were performed at 
## (ideally) before treatment initiation, 
## at 14 days of therapy, 
## at the times of smear conversion and following culture conversion (not available in the MDR-RVC), 
## at 6 months and/or therapy end in patients with DS-TB, 
## and additionally at 10, 15 and 20 months of therapy in patients with MDR-TB. 
# After completion of 4 weeks of therapy, an additional study visit was performed in patients from the MDR-RVC. 
# All patients completed 12 months of evaluation following the end of therapy to capture disease recurrence. 


# GSE147691 (Romanian Validation cohort) ----------------------------------------------------------------------------------------------------------------------------------------
this.accession.no <- "GSE147691"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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
this.accession.no <- "GSE147690"
setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

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

gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1B")) #NR_027484 is FCGR1C on GPL13497, but has main name FCGR1B
signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
signature_geneid %in% row.names(counts_norm)
#CANNOT FIND FCGR1C OR ITS ALIASES, used FCGR1B because the glp shows them as same


if(length(signature_geneid) < 6){
  print("Missing gene in signature after genone_annot conversion")
  stop() }

gene_set_list <- list(c(signature_geneid))

gsvapar <- gsvaParam(as.matrix(counts_norm),
                     gene_set_list,
                     maxDiff = TRUE,
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical$group))


gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 

my_comparisons <- combn(unique(clinical$group), 2, simplify = FALSE)


x_order <- c("0", "14 Days", "28 Days", "6 Mths",  "10 Mths",  "15 Mths",  "20 Mths")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)
boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(V1 ~ group,
              paired = FALSE) %>%
  add_xy_position(x = "group") 

stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  y = as.numeric(boxplot_gsva[,1]),
  group = group)) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5, 
              width = 0.3) +
  
  {if(nrow(stat.table.gsva) >0 )
    stat_pvalue_manual(stat.table.gsva,
                       label = "p",
                       tip.length = 0.01,
                       size = 5)
  } + 
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  stat_summary(fun = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  theme(axis.text.x = element_text(size = 18))+
  labs(title = paste0("Signature Analysis: ", this.accession.no),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1B"
       # caption = "Signature:S100A8"
       # caption = paste0("Signature: ", paste0(gene_set_list[[1]], collapse = " "))
       
       
  ) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

this.accession.res.dir <- file.path(output.dir, this.accession.no)
if(!exists(this.accession.res.dir)) dir.create(this.accession.res.dir)

this.figure.dir <- file.path(this.accession.res.dir, "figures")
if(!exists(this.figure.dir)) dir.create(this.figure.dir)

ggsave(boxplotfinal2, filename = file.path(this.figure.dir, paste0("gsva_plot_", this.accession.no, ".png")),
       width = 3600,
       height = 3200,
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


# Create a list to store AUC values and roc objects
res_table <- data.frame()
roc_objects <- list()

# Loop through each pairwise comparison
for (pair in pairwise_comparisons) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Subset data to omly include the 2 rgroups of interest
  subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
  subset_counts <- counts_norm[, row.names(subset_clinical)]
  
  subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
  
  # GSVA
  gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1B")) #NR_027484 is FCGR1C on GPL13497, but has main name FCGR1B
  signature_geneid <- as.character(gene_annot[match(gene_set_list[[1]], gene_annot$GENE_SYMBOL), "REFSEQ"])
  gene_set_list <- list(c(signature_geneid))
  
  gsvapar <- gsvaParam(as.matrix(subset_counts), #counts_vst$E is the same
                       gene_set_list, 
                       maxDiff = TRUE, 
                       minSize = 1)
  
  gsva_res <- gsva(gsvapar)
  
  glm_data <- data.frame(Score = gsva_res[1,], Group = subset_clinical$group)
  
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
  
  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
}

saveRDS(roc_objects, file.path(this.accession.res.dir, paste0(this.accession.no, "_roc_objects.rds")))
write.csv(res_table, file.path(this.accession.res.dir, paste0(this.accession.no,"_res_table.csv")))


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


roc_data$Comparison <- factor(roc_data$Comparison)


roc_data$ci <- res_table[match(roc_data$Comparison, res_table$comparison), "ci"]

roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
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
    title = "ROC - MDR-TB (Romanian Validation Cohort",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1B") 

ggsave(timepoint_roc, filename = file.path(this.figure.dir, "timepoint_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")
