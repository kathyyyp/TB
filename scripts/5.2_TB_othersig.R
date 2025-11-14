#Gene sets from publicly available data

# Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
# my_directory <- "/Volumes/One Touch/RBMB"

setwd(file.path(my_directory,"TB"))

.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")
# .libPaths("/Volumes/One Touch/rlibrary")

# Load packages
library(GSVA)
library(pROC)
library(caret)
library(stats)
library(ranger)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(rstatix)
library(ggpubr)

output.dir<- file.path(my_directory,"TB", "output", "other_signatures")



# Define results directories HERE
# current_roc_dir <- TBT0vsHCT0.dir
# load("workspaces/03_03_25_bestandworstauc.Rdata")
# write.table(clinical, file.path("data", "processed", "clinical.txt"))
# write.table(expression, file.path("data", "processed", "expression.txt"))

clinical <- read.delim(file.path("data", "processed", "clinical.txt"), sep = " ")
expression <- as.matrix(read.delim(file.path("data", "processed", "expression.txt"), sep = " ", check.names = FALSE)) #normalised by cyclicloess

Sweeney <- c("GBP5", "DUSP3", "KLF2")
Penn_RISK6 <- c("GBP2", "TUBGCP6", "TRMT2A", "SDR39U1", "SERPING1", "FCGR1B")
Jen_Ho <- c("FCGR1CP", "IFI35", "IFI6", "GZMA", "DHRS9", "TGIF1", "APOL6")
Zak <- c("GBP1", "GBP4", "GBP5", "FCGR1A", "FCGR1B", "TRAFD1", "SERPING1", "SEPTIN4", "ETV7", "BATF2", "TAP1", "ANKRD22", "APOL1", "ID3")
Thakur <- c("GBP1", "FCGR1A", "SERPING1", "BATF2", "BCL6", "AIM2", "SMARCD3", "ANXA3", "SOCS3") #replaced BATF1 with BATF2 bc our study only has BATF2
Darboe <- c("GBP1", "GBP2", "GBP5", "FCGR1B", "TRAFD1", "SERPING1", "ETV7", "BATF2", "TAP1", "ANKRD22", "APOL1")
Maertzdorf <- c("GBP1", "IFITM3", "P2RY14", "ID3")
Chendi <- c("GBP1", "MPPE1", "BATF2") # CD1C replaced with MPPE1

gene_sets <- list(Sweeney, Penn_RISK6, Jen_Ho, Zak, Thakur, Darboe, Maertzdorf, Chendi)
names(gene_sets) <- c("Sweeney", "Penn_RISK6", "Jen_Ho", "Zak", "Thakur", "Darboe", "Maertzdorf", "Chendi")

row.names(expression)[which(row.names(expression) == "TRAMT2A")] <- "TRMT2A" #this is a type in our expression file

#Counts in the signatures that are in our data
overlap <- lapply(gene_sets, function(genes) {
  intersect(genes, row.names(expression)) # dataset_genes = list of genes in your dataset
})

#Genes that are missing
missing_genes <- sapply(gene_sets, function(genes) {
  genes[!(genes %in% row.names(expression))]
})

missing_counts <- sapply(gene_sets, function(genes) {
  sum(!(genes %in% row.names(expression)))
})

data.frame(gene_sets = names(gene_sets), Missing_Genes = missing_counts)



# DISEASE -----------------------------------------
## 1) Subset
#Make Disease a factor
Disease <- clinical$Disease
expression_auc <- expression #rows are samples instead of genes here

#Subset
all(colnames(expression_auc) == row.names(clinical))
T0_samples <- row.names(clinical)[which(clinical$condition == "HC_T0" | clinical$condition == "TB_T0")]


expression_auc <- expression_auc[,T0_samples]
clinical_auc <- clinical[T0_samples,]
all(colnames(expression_auc) == row.names(clinical_auc))


## 2) GSVA  ---------------- 
#genes must be rows
# GSVA generally works well with log2-transformed intensities from microarrays or log-CPM/RPKM/TPM values from RNA-seq.
#Your cyclic loess-normalized data appears to be already in log scale, so you do NOT need to log-transform it again.
gsvapar <- gsvaParam(expression_auc,
                     gene_sets, 
                     maxDiff = TRUE, 
                     minSize = 1)


gsva_res <- gsva(gsvapar)


## 3) GSVA BOXPLOT - DISEASE ---------------

figures.dir <- file.path(output.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)

gsva.res.dir <- file.path(figures.dir, "gsva")
if(!exists(gsva.res.dir)) dir.create(gsva.res.dir)

this.res.dir <- file.path(gsva.res.dir, "disease")
if(!exists(this.res.dir)) dir.create(this.res.dir)


all(colnames(gsva_res) == row.names(clinical_auc))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical_auc$Disease))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- combn(unique(boxplot_gsva$group), 2, simplify = FALSE)


x_order <- c("HC", "TB")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)

boxplot_gsva[,1:(ncol(boxplot_gsva)-1)] <- sapply(boxplot_gsva[,1:(ncol(boxplot_gsva)-1)], as.numeric)

for (i in colnames(boxplot_gsva)[1:(ncol(boxplot_gsva) - 1)]){
  

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(as.formula(paste(i ,"~ group")),
              paired = FALSE) %>%
  add_xy_position(x = "group")

# stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva[, i]) + 0.05*(max(boxplot_gsva[, i]))
stat.table.gsva$y.position <- lowest_bracket


boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,i]),
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
  
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+

  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  theme(axis.text.x = element_text(size = 15))+
  labs(title = paste0("TB Signature Analysis: ", i),
       caption = str_wrap(paste("Signature:", paste(overlap[[i]], collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")


ggsave(boxplotfinal2, filename = file.path(this.res.dir, paste0("gsva_plot_", i, ".png")), 
       width = 2300, 
       height = 2500, 
       units = "px" )

}


## 4) ROC  ------

disease_full_tableofres <- data.frame()
disease_tableofres <- data.frame()
disease_roc_objects <- list()
for (i in 1:nrow(gsva_res)){
glm_data <- data.frame(Score = gsva_res[i,], Group = clinical_auc$Disease)

table(glm_data$Group)


glm_data$Group <- factor(glm_data$Group, levels = c("HC", "TB"))
glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 

test_probs <- predict(glm_model, type = "response")

roc_obj <- roc( glm_data$Group, test_probs)
auc(roc_obj)

auc_ci <- ci.auc(roc_obj)

optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))

res <-cbind(
    signature = row.names(gsva_res)[i],
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
 
disease_tableofres <- rbind(disease_tableofres,res)
disease_roc_objects[[row.names(gsva_res)[i]]] <- roc_obj

}

disease_full_tableofres <- rbind(disease_tableofres, "TBT0vsHCT0")

write.csv(disease_full_tableofres, file.path(output.dir,"disease_full_tableofres.csv"))

## 5) ROC Curves -----------------------------------------------------------
roc.res.dir <- file.path(figures.dir, "roc")
if(!exists(roc.res.dir)) dir.create(roc.res.dir)

this.res.dir <- file.path(roc.res.dir, "disease")
if(!exists(this.res.dir)) dir.create(this.res.dir)


 # Convert ROC data to a format suitable for ggplot
 
roc_data <- do.call(rbind, lapply(names(disease_roc_objects), function(comparison) {
   data.frame(
     TPR = rev(disease_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
     FPR = rev(1 - disease_roc_objects[[comparison]]$specificities),  # False Positive Rate
     Signature = comparison,
     auc = rev(disease_roc_objects[[comparison]]$auc)
   )
 }))


 roc_data$ci <- disease_tableofres[match(roc_data$Signature, disease_full_tableofres$signature), "ci"]

 roc_data$legend <- paste0(roc_data$Signature,": \n AUC = ", 
                                  round(roc_data$auc, 2), " (", roc_data$ci, ")")

 
 disease_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
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
     title = paste("ROC (Publicly available signatures) - Control vs TB"),
     x = "FPR (1 - Specificity)",
     y = "TPR (Sensitivity)",
     color = "Signature") 
 
 ggsave(disease_roc, filename = file.path(this.res.dir, "disease_roc.png"), 
        width = 3000, 
        height = 3200, 
        units = "px")
 
 
 
 
# TIMEPOINT------------------------------------------------------------------------------------------------------------------
# A) FOR BOXPLOT -----------
## 1) Subset (all timepoints) ------
Timepoint <- clinical$Timepoint
expression_auc <- expression #rows are samples instead of genes here

#Subset
all(colnames(expression_auc) == row.names(clinical))
TB_samples <- row.names(clinical)[which(clinical$Disease == "TB")]


expression_auc <- expression_auc[,TB_samples]
clinical_auc <- clinical[TB_samples,]

## 2) GSVA (all timepoints) --------
gsvapar <- gsvaParam(expression_auc,
                     gene_sets, 
                     maxDiff = TRUE, 
                     minSize = 1)


gsva_res <- gsva(gsvapar)

## 3) GSVA BOXPLOT (all timepoints) --------

this.res.dir <- file.path(gsva.res.dir, "Timepoint")
if(!exists(this.res.dir)) dir.create(this.res.dir)


all(colnames(gsva_res) == row.names(clinical_auc))

boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = clinical_auc$Timepoint))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- combn(unique(boxplot_gsva$group), 2, simplify = FALSE)


x_order <- c("T0", "T2", "T4", "T6")

boxplot_gsva$group <- factor(boxplot_gsva$group, levels = x_order)

boxplot_gsva[,1:(ncol(boxplot_gsva)-1)] <- sapply(boxplot_gsva[,1:(ncol(boxplot_gsva)-1)], as.numeric)

for (i in colnames(boxplot_gsva)[1:(ncol(boxplot_gsva) - 1)]){
  

stat.table.gsva <- boxplot_gsva  %>%
  wilcox_test(as.formula(paste(i ,"~ group")),
              paired = FALSE) %>%
  add_xy_position(x = "group")

if(any(stat.table.gsva$p < 0.05)){
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva[, i]) + 0.05*(max(boxplot_gsva[, i]))
stat.table.gsva$y.position <- lowest_bracket
stat.table.gsva$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva))

} else{
  stat.table.gsva <- stat.table.gsva[0,]
}



boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group, level = x_order),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,i]),
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

  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  theme(axis.text.x = element_text(size = 15))+
  labs(title = paste0("TB Signature Analysis: ", i),
       caption = str_wrap(paste("Signature:", paste(overlap[[i]], collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Disease")

if (nrow(stat.table.gsva) > 0) {
  boxplotfinal2 <- boxplotfinal2 + 
    stat_pvalue_manual(stat.table.gsva, label = "p", tip.length = 0.01, size = 6)
}

ggsave(boxplotfinal2, filename = file.path(this.res.dir, paste0("gsva_plot_", i, ".png")), 
       width = 2300, 
       height = 2500, 
       units = "px" )

}




# B) FOR ROC ----------------------

roc.res.dir <- file.path(figures.dir, "roc")
if(!exists(roc.res.dir)) dir.create(roc.res.dir)

this.res.dir <- file.path(roc.res.dir, "timepoint")
if(!exists(this.res.dir)) dir.create(this.res.dir)

## 1) Subset (per comparison) ------------

timepoint_full_tableofres <- data.frame()

for (study in names(gene_sets)){

  this_gene_set <- list(gene_sets[[study]])
case_and_control <- data.frame(case = c("TB_T2", "TB_T4", "TB_T6", "HC_T6"), control = c("TB_T0", "TB_T0", "TB_T0", "HC_T0"))

  
timepoint_tableofres <- data.frame()
timepoint_roc_objects <- list()

for (i in 1:4){
print(paste("Comparison", i))
  
case = case_and_control$case[i]
control = case_and_control$control[i]
  
TB_timepoint_subset <- row.names(clinical)[which(clinical$condition == case | clinical$condition == control )]

#Subset
expression_auc <- expression[,TB_timepoint_subset]
clinical_auc <- clinical[TB_timepoint_subset,]

## 2) GSVA (per comparison) and ROC --------
gsvapar <- gsvaParam(expression_auc,
                     this_gene_set, 
                     maxDiff = TRUE, 
                     minSize = 1)


gsva_res <- gsva(gsvapar)


  glm_data <- data.frame(Score = gsva_res[1,], Group = clinical_auc$condition)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c(control,case))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(glm_data$Group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  


  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  res <-cbind(
    comparison = paste0(control,"vs",case),
    samples_case = paste(case, "=", sum(glm_data$Group == case)),
    samples_control = paste(control, "=", sum(glm_data$Group == control)),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
  )
  
  
 timepoint_tableofres <- rbind(timepoint_tableofres, res) #contains results for all 4 comparisons for this study
 
timepoint_roc_objects[[paste0(control,"vs",case)]] <- roc_obj

  
} #End of this comparison 
timepoint_tableofres <- cbind(timepoint_tableofres, study)




 ## 3) ROC curve --------
roc_data <- do.call(rbind, lapply(names(timepoint_roc_objects), function(comparison) {
   data.frame(
     TPR = rev(timepoint_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
     FPR = rev(1 - timepoint_roc_objects[[comparison]]$specificities),  # False Positive Rate
     Comparison = comparison,
     auc = rev(timepoint_roc_objects[[comparison]]$auc)
   )
 }))
 
 roc_data$Comparison <- factor(roc_data$Comparison)
 
 roc_data$ci <- timepoint_tableofres[match(roc_data$Comparison, timepoint_tableofres$comparison), "ci"]

 roc_data$legend <- paste0(roc_data$Comparison,": \n AUC = ", 
                                  round(roc_data$auc, 2), " (", roc_data$ci, ")")

 

timepoint_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
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
     title = paste0("ROC - TB treatment timepoints (",study, ")"),
     x = "FPR (1 - Specificity)",
     y = "TPR(Sensitivity)",
     color = "Comparison",
     caption = str_wrap(paste("Signature:", paste(overlap[[study]], collapse = " "))))

 ggsave(timepoint_roc, filename = file.path(this.res.dir, paste(study,"timepoint_roc.png")), 
        width = 3200, 
        height = 3400, 
        units = "px")
 

 timepoint_full_tableofres <- rbind(timepoint_full_tableofres, timepoint_tableofres)
 
}

write.csv(timepoint_full_tableofres, file.path(output.dir,"timepoint_full_tableofres.csv"))

