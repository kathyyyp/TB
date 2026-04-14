# SCRIPT 3: Testing the performance of gene sets from publicly available data (provided by Tess) in our discovery cohort data
# Using both GSVA method and also Mean_zscores method
# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

# Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")


# Mac
my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("/Volumes/One Touch/RBMB/rlibrary")

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


# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
output.dir <- file.path(main.dir, "output_v2")
output.othersig.dir<- file.path(my_directory,"TB", "output_v2", "2.2_other_published_signatures")
if(!exists(output.othersig.dir)) dir.create(output.othersig.dir)

# ================================================================================== #
# 1. LOAD IN DATA & GET UPREGULATED GENES ==========================================================
# ================================================================================== #
clinical <- read.csv(file.path(main.dir, "data/processed/post-QC/clinical.csv"), row.names = 1)
expression <- read.csv(file.path(main.dir, "data/processed/post-QC/expression.csv"), row.names = 1, check.names = FALSE)

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

row.names(expression)[which(row.names(expression) == "TRAMT2A")] <- "TRMT2A" #this is a typo in our expression file

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




#==============================================================================================================#
#======================================= DISEASE (TB_T0 vs HC_T0) =============================================
#============================================================================================================== #
### 1) Subset
#Make Disease a factor
Disease <- clinical$Disease
expression_auc <- expression #rows are samples instead of genes here

#Subset
all(colnames(expression_auc) == row.names(clinical))
T0_samples <- row.names(clinical)[which(clinical$condition == "HC_T0" | clinical$condition == "TB_T0")]

expression_auc <- expression_auc[,T0_samples]
expression_auc <- as.matrix(expression_auc)
clinical_auc <- clinical[T0_samples,]
all(colnames(expression_auc) == row.names(clinical_auc))






#Get scores for test set-----------------------------------------------------------------------------------------
#Subset for our x number of genes 
gene_set_exp <- t(expression_auc)

# Below is how to claculate mean of z scores for ONE signature
# I have multiple signatures and want a mean score per sample, for every signature
# Intended output is table with rows as samples and columns are signature name , the mean scores are the values
# gene_set_exp <- gene_set_exp[,gene_set_list] #genes as columns
# gene_set_zscore<-scale(gene_set_exp, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
# #take the mean of the scaled+centred data for each gene and use this as each sample's score
# mean_sig_zscore<- rowMeans(gene_set_zscore)


# All scores collected
score_list <- lapply(gene_sets, function(set) {
  this_gene_set_exp <- gene_set_exp[, set]
  zscore <- scale(this_gene_set_exp, center = TRUE, scale = TRUE)
  rowMeans(zscore)
})

mean_sig_zscore <- as.data.frame(score_list)

## 3) Mean Z-Score BOXPLOT - DISEASE ---------------
meanzscore.res.dir <- file.path(output.othersig.dir, "mean_zscore")
if(!exists(meanzscore.res.dir)) dir.create(meanzscore.res.dir)

disease.res.dir <- file.path(meanzscore.res.dir, "disease")
if(!exists(disease.res.dir)) dir.create(disease.res.dir)

this.boxplot.dir <- file.path(disease.res.dir, "boxplot")
if(!exists(this.boxplot.dir)) dir.create(this.boxplot.dir)

all(row.names(mean_sig_zscore) == row.names(clinical_auc))

boxplot<- as.data.frame(cbind(mean_sig_zscore,
                                    group = clinical_auc$Disease))



boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None",
                    plot.caption = element_text( size = 12, hjust=0),
                    plot.caption.position = "plot")  #align caption to entire plot including margins


my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)


x_order <- c("HC", "TB")

boxplot$group <- factor(boxplot$group, levels = x_order)

boxplot$group <- factor(boxplot$group, levels = x_order)

boxplot[,1:(ncol(boxplot)-1)] <- sapply(boxplot[,1:(ncol(boxplot)-1)], as.numeric)

for (i in colnames(boxplot)[1:(ncol(boxplot) - 1)]){
  
  stat.table <- boxplot  %>%
    wilcox_test(as.formula(paste(i ,"~ group")),
                paired = FALSE) %>%
    add_xy_position(x = "group")
  
  # stat.table <- stat.table[which(stat.table$p < 0.05),]
  lowest_bracket <- max(boxplot[, i]) + 0.05*(max(boxplot[, i]))
  stat.table$y.position <- lowest_bracket
  
  
  boxplotfinal2 <- ggplot(boxplot, aes(
    x = factor(group, level = x_order),
    # x = factor(group),
    y = as.numeric(boxplot[,i]),
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
    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 6)+
    
    # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
    # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
    
    theme(axis.text.x = element_text(size = 15))+
    labs(title = paste0("TB Signature Analysis: ", i),
         caption = paste0(str_wrap(paste("Signature:", paste(overlap[[i]], collapse = ", "))), "\n",
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "P values from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Disease")
  
  
  ggsave(boxplotfinal2, filename = file.path(this.boxplot.dir, paste0("meanzscore_plot_", i, ".png")), 
         width = 2300, 
         height = 2500, 
         units = "px" )
  
}

#Combined boxplot
boxplot_long <- pivot_longer(boxplot, cols = colnames(boxplot)[1:(ncol(boxplot) - 1)],
             names_to = "signature",
             values_to = "mean_sig_zscore")

x_order <- names(gene_sets)
boxplot_long$signature <- factor(boxplot_long$signature, levels = x_order)

#Get P-values for WITHIN signature comparisons (HC vs TB within each signature)
stat.table <- boxplot_long  %>%
  group_by(signature) %>% #this is how you want the plot to be grouped ie. the x-axis label
  wilcox_test(mean_sig_zscore ~ group, #this is the variable that will be in the legend
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "signature") #this is the x-axis/the grouping variable

y_pos <- boxplot_long %>%
  group_by(signature) %>%
  summarise(y.position = max(mean_sig_zscore))

  # # stat.table <- stat.table[which(stat.table$p < 0.05),]
  stat.table$y.position <- y_pos$y.position + 0.12

    boxplot_long <- as.data.frame(boxplot_long)

  stat.table$signature <- factor(stat.table$signature, levels = x_order)
  
boxplotcombined <- ggplot(boxplot_long, aes(
    x = signature,
    y = mean_sig_zscore,
    color = group)) +
    
    theme_bw()+
    
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 24),
          axis.text.x = element_text(angle = 45, hjust = 1),
          title = element_text(size = 20),
          legend.position = "bottom",
          legend.text = element_text(size = 18),
          plot.caption = element_text( size = 18, hjust=0),
          plot.caption.position = "plot") + 
    
    geom_boxplot(position = position_dodge(0.75),
                 width = 0.4) +
    
    # geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
    #            alpha=0.5,
    #           size = 2) +

    stat_summary(fun.y = mean, fill = "red",
                 geom = "point", shape = 21, size =4,
                 position = position_dodge(width = 0.75),
                 show.legend = FALSE) +
    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 6)+
    
    labs(title = paste0("Performance of Published Signatures"),
         color = "Disease",
         caption = paste0(
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "P values from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Signature")
  
  
  ggsave(boxplotcombined, filename = file.path(this.boxplot.dir, paste0("mean_zscore_plot_combined.png")), 
         width = 5000, 
         height = 3000, 
         units = "px" )
  
  
  
## 4) ROC  ------
this.roc.dir <- file.path(disease.res.dir, "roc")
if(!exists(this.roc.dir)) dir.create(this.roc.dir)

this.forest.dir <- file.path(disease.res.dir, "forest")
if(!exists(this.forest.dir)) dir.create(this.forest.dir)

disease_full_tableofres <- data.frame()
disease_tableofres <- data.frame()
disease_roc_objects <- list()

forestplot_full_tableofres <- data.frame()
forestplot_tableofres <- data.frame()

mean_sig_zscore_t <- t(mean_sig_zscore)
for (i in 1:nrow(mean_sig_zscore_t)){
  glm_data <- data.frame(Score = mean_sig_zscore_t[i,], Group = clinical_auc$Disease)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c("HC", "TB"))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc( glm_data$Group, test_probs)
  auc(roc_obj)
  
  auc_ci <- ci.auc(roc_obj)
  
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  res <-cbind(
    signature = row.names(mean_sig_zscore_t)[i],
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  disease_tableofres <- rbind(disease_tableofres,res)
  disease_roc_objects[[row.names(mean_sig_zscore_t)[i]]] <- roc_obj
  
  
    # Sensitivity confidence interval (get the sensitivity at the optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"]))
  
  # Specificity confidence interval (get the specificity at the optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  
  forestplot_res<- cbind(
    signature = row.names(mean_sig_zscore_t)[i],
    
    auc = auc(roc_obj),
    auc_ci_low = as.numeric(auc_ci[1]),
    auc_ci_high = as.numeric(auc_ci[3]),
    
    sensitivity = optimal_threshold_coords$sensitivity, 
    sensitivity_ci_low = ci_sens[, "2.5%"],
    sensitivity_ci_high = ci_sens[, "97.5%"],
    
    specificity = optimal_threshold_coords$specificity,
    specificity_ci_low = ci_spec[, "2.5%"],
    specificity_ci_high = ci_spec[, "97.5%"])
  
  
  forestplot_tableofres <- rbind(forestplot_tableofres,forestplot_res)
  
}

disease_full_tableofres <- rbind(disease_tableofres, "TBT0vsHCT0")
forestplot_full_tableofres <- rbind(forestplot_tableofres, "TBT0vsHCT0")

write.csv(disease_full_tableofres, file.path(this.roc.dir,"mean_zscore_disease_other_published_signatures_tableofauc.csv"))
write.csv(forestplot_full_tableofres, file.path(this.forest.dir,"mean_zscore_disease_other_published_signatures_forestplotres.csv"))

## 5) ROC Curves -----------------------------------------------------------
this.roc.dir <- file.path(disease.res.dir, "roc")
if(!exists(this.roc.dir)) dir.create(this.roc.dir)

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


#reorder signature so colours are how i want
roc_data$Signature <- factor(roc_data$Signature, levels = names(gene_sets))

#reorder legend to match
roc_data <- roc_data %>%
  mutate(
    #get first word of signature and remove the :
    legend_first_word = str_remove(word(legend, 1),":"),
    # rorder legend factor based on signature factor order
    legend = factor(legend, levels = legend[match(levels(Signature), legend_first_word)])
  )
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
        title = element_text(size = 20),
        plot.caption = element_text( size = 12, hjust=0),
        plot.caption.position = "plot" 
) +
  labs(
    title = paste("ROC (Publicly available signatures) - Control vs TB"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Signature",
    caption = paste0("ROC showing performance of a logistic regression model predicting TB_T0 vs HC_TO using the signature score \n",
    "Signature scores calculated as mean of Z-scored expression of signature genes"))

ggsave(disease_roc, filename = file.path(this.roc.dir, "disease_roc.png"), 
       width = 3000, 
       height = 3200, 
       units = "px")



## 6) Forest plot -----------------------------------------------------------

forestplot_full_tableofres

  forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 
  

  # other gene signatures don't need to be plotted, can rename axis to HC_T0 vs TB_T0
  res_table <- forestplot_full_tableofres[-nrow(forestplot_full_tableofres),]
  res_table$signature <- factor(res_table$signature, levels = names(gene_sets))
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)

  auc_plot <-  res_table %>% 
    ggplot(aes(y = signature)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = signature), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = signature),
                   size = 1) +
    forestplot_theme +
    scale_x_continuous(limits = c(0.1, 1) )+
    ylab(NULL)+
    xlab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(y = signature)) + 
      theme_bw() +
      geom_point(aes(x=sensitivity, color = signature), shape=15, size=3) +
      geom_linerange(aes(xmin=sensitivity_ci_low, 
                         xmax=sensitivity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
    scale_x_continuous(limits = c(0.1, 1) )+      
       ylab(NULL)+
      xlab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(y = signature)) + 
      theme_bw() +
      geom_point(aes(x=specificity, color = signature), shape=15, size=3) +
      geom_linerange(aes(xmin=specificity_ci_low, 
                         xmax=specificity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
    scale_x_continuous(limits = c(0.1, 1) )+ 
      ylab(NULL)+
      xlab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3)
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 24),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating TB_T0 vs HC_T0 \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_z_score_forestplot.png"),
         width = 15, height = 20, units = "cm",   bg = "white"  )
  

#FLIPPED PLOT
  auc_plot <-  res_table %>% 
    ggplot(aes(x = signature)) + 
    theme_bw() +
    geom_point(aes(y=auc, color = signature), shape=15, size=3) +
    geom_linerange(aes(ymin=auc_ci_low, 
                       ymax=auc_ci_high, 
                       color = signature),
                   size = 1) +
    forestplot_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0.1, 1) )+
    xlab(NULL)+
    ylab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(x = signature)) + 
      theme_bw() +
      geom_point(aes(y=sensitivity, color = signature), shape=15, size=3) +
      geom_linerange(aes(ymin=sensitivity_ci_low, 
                         ymax=sensitivity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
           theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0.1, 1) )+      
       xlab(NULL)+
      ylab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(x = signature)) + 
      theme_bw() +
      geom_point(aes(y=specificity, color = signature), shape=15, size=3) +
      geom_linerange(aes(ymin=specificity_ci_low, 
                         ymax=specificity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0.1, 1) )+ 
      xlab(NULL)+
      ylab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3,
                              heights = c(1,1,1.4))
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 18),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating TB_T0 vs HC_T0 \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_z_score_forestplot_flip.png"),
         width = 12, height = 25, units = "cm",   bg = "white"  )
  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#==============================================================================================================#
#============================================= HC_T6 vs HC_T0) ================================================
#============================================================================================================== #
### 1) Subsetem
#Make condition a factor
Condition <- clinical$condition
expression_auc <- expression #rows are samples instead of genes here

#Subset
all(colnames(expression_auc) == row.names(clinical))
HC_samples <- row.names(clinical)[which(clinical$condition == "HC_T0" | clinical$condition == "HC_T6")]

expression_auc <- expression_auc[,HC_samples]
expression_auc <- as.matrix(expression_auc)
clinical_auc <- clinical[HC_samples,]
all(colnames(expression_auc) == row.names(clinical_auc))






#Get scores for test set-----------------------------------------------------------------------------------------
#Subset for our x number of genes 
gene_set_exp <- t(expression_auc)

# Below is how to claculate mean of z scores for ONE signature
# I have multiple signatures and want a mean score per sample, for every signature
# Intended output is table with rows as samples and columns are signature name , the mean scores are the values
# gene_set_exp <- gene_set_exp[,gene_set_list] #genes as columns
# gene_set_zscore<-scale(gene_set_exp, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
# #take the mean of the scaled+centred data for each gene and use this as each sample's score
# mean_sig_zscore<- rowMeans(gene_set_zscore)


# All scores collected
score_list <- lapply(gene_sets, function(set) {
  this_gene_set_exp <- gene_set_exp[, set]
  zscore <- scale(this_gene_set_exp, center = TRUE, scale = TRUE)
  rowMeans(zscore)
})

mean_sig_zscore <- as.data.frame(score_list)

## 3) Mean Z-Score BOXPLOT - HCT0 vs HCT6 ---------------
meanzscore.res.dir <- file.path(output.othersig.dir, "mean_zscore")
if(!exists(meanzscore.res.dir)) dir.create(meanzscore.res.dir)

hc.res.dir <- file.path(meanzscore.res.dir, "HCT6_vs_HCT0")
if(!exists(hc.res.dir)) dir.create(hc.res.dir)

this.boxplot.dir <- file.path(hc.res.dir, "boxplot")
if(!exists(this.boxplot.dir)) dir.create(this.boxplot.dir)

all(row.names(mean_sig_zscore) == row.names(clinical_auc))

boxplot<- as.data.frame(cbind(mean_sig_zscore,
                                    group = clinical_auc$condition))



boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None",
                    plot.caption = element_text( size = 12, hjust=0),
                    plot.caption.position = "plot")  #align caption to entire plot including margins


my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)


x_order <- c("HC_T6", "HC_T0")

boxplot$group <- factor(boxplot$group, levels = x_order)

boxplot[,1:(ncol(boxplot)-1)] <- sapply(boxplot[,1:(ncol(boxplot)-1)], as.numeric)

for (i in colnames(boxplot)[1:(ncol(boxplot) - 1)]){
  
  stat.table <- boxplot  %>%
    wilcox_test(as.formula(paste(i ,"~ group")),
                paired = FALSE) %>%
    add_xy_position(x = "group")
  
  # stat.table <- stat.table[which(stat.table$p < 0.05),]
  lowest_bracket <- max(boxplot[, i]) + 0.05*(max(boxplot[, i]))
  stat.table$y.position <- lowest_bracket
  
  
  boxplotfinal2 <- ggplot(boxplot, aes(
    x = factor(group, level = x_order),
    # x = factor(group),
    y = as.numeric(boxplot[,i]),
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
    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 6)+
    
    # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
    # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
    
    theme(axis.text.x = element_text(size = 15))+
    labs(title = paste0("TB Signature Analysis: ", i),
         caption = paste0(str_wrap(paste("Signature:", paste(overlap[[i]], collapse = ", "))), "\n",
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "P values from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Disease")
  
  
  ggsave(boxplotfinal2, filename = file.path(this.boxplot.dir, paste0("meanzscore_plot_", i, ".png")), 
         width = 2300, 
         height = 2500, 
         units = "px" )
  
}

#Combined boxplot
boxplot_long <- pivot_longer(boxplot, cols = colnames(boxplot)[1:(ncol(boxplot) - 1)],
             names_to = "signature",
             values_to = "mean_sig_zscore")

x_order <- names(gene_sets)
boxplot_long$signature <- factor(boxplot_long$signature, levels = x_order)

#Get P-values for WITHIN signature comparisons (HC vs TB within each signature)
stat.table <- boxplot_long  %>%
  group_by(signature) %>% #this is how you want the plot to be grouped ie. the x-axis label
  wilcox_test(mean_sig_zscore ~ group, #this is the variable that will be in the legend
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "signature") #this is the x-axis/the grouping variable

y_pos <- boxplot_long %>%
  group_by(signature) %>%
  summarise(y.position = max(mean_sig_zscore))

  # # stat.table <- stat.table[which(stat.table$p < 0.05),]
  stat.table$y.position <- y_pos$y.position + 0.12

    boxplot_long <- as.data.frame(boxplot_long)

  stat.table$signature <- factor(stat.table$signature, levels = x_order)
  
boxplotcombined <- ggplot(boxplot_long, aes(
    x = signature,
    y = mean_sig_zscore,
    color = group)) +
    
    theme_bw()+
    
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 24),
          axis.text.x = element_text(angle = 45, hjust = 1),
          title = element_text(size = 20),
          legend.position = "bottom",
          legend.text = element_text(size = 18),
          plot.caption = element_text( size = 18, hjust=0),
          plot.caption.position = "plot") + 
    
    geom_boxplot(position = position_dodge(0.75),
                 width = 0.4) +
    # 
    # geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
    #            alpha=0.5,
    #           size = 2) +

    stat_summary(fun.y = mean, fill = "red",
                 geom = "point", shape = 21, size =4,
                 position = position_dodge(width = 0.75),
                 show.legend = FALSE) +
    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 6)+
    
    labs(title = paste0("Performance of Published Signatures"),
         color = "Disease",
         caption = paste0(
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "P values from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Signature")
  
  
  ggsave(boxplotcombined, filename = file.path(this.boxplot.dir, paste0("mean_zscore_plot_combined.png")), 
         width = 5000, 
         height = 3000, 
         units = "px" )
  
  
  
## 4) ROC  ------
this.roc.dir <- file.path(hc.res.dir, "roc")
if(!exists(this.roc.dir)) dir.create(this.roc.dir)

this.forest.dir <- file.path(hc.res.dir, "forest")
if(!exists(this.forest.dir)) dir.create(this.forest.dir)

hc_full_tableofres <- data.frame()
hc_tableofres <- data.frame()
hc_roc_objects <- list()

forestplot_full_tableofres <- data.frame()
forestplot_tableofres <- data.frame()

mean_sig_zscore_t <- t(mean_sig_zscore)
for (i in 1:nrow(mean_sig_zscore_t)){
  print(i)
  glm_data <- data.frame(Score = mean_sig_zscore_t[i,], Group = clinical_auc$condition)
  
  table(glm_data$Group)
  
  
  glm_data$Group <- factor(glm_data$Group, levels = c("HC_T0", "HC_T6"))
  glm_model <- glm(Group ~ Score, data = glm_data, family = binomial) 
  
  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc( glm_data$Group, test_probs)
  auc(roc_obj)
  
  auc_ci <- ci.auc(roc_obj)
  
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  res <-cbind(
    signature = row.names(mean_sig_zscore_t)[i],
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  hc_tableofres <- rbind(hc_tableofres,res)
  hc_roc_objects[[row.names(mean_sig_zscore_t)[i]]] <- roc_obj
  
  
    # Sensitivity confidence interval (get the sensitivity at the optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"]))
  
  # Specificity confidence interval (get the specificity at the optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  
  forestplot_res<- cbind(
    signature = row.names(mean_sig_zscore_t)[i],
    
    auc = auc(roc_obj),
    auc_ci_low = as.numeric(auc_ci[1]),
    auc_ci_high = as.numeric(auc_ci[3]),
    
    sensitivity = optimal_threshold_coords$sensitivity, 
    sensitivity_ci_low = ci_sens[, "2.5%"],
    sensitivity_ci_high = ci_sens[, "97.5%"],
    
    specificity = optimal_threshold_coords$specificity,
    specificity_ci_low = ci_spec[, "2.5%"],
    specificity_ci_high = ci_spec[, "97.5%"])
  
  
  forestplot_tableofres <- rbind(forestplot_tableofres,forestplot_res)
  
}

hc_full_tableofres <- cbind(hc_tableofres, "HCT6vsHCT0")
forestplot_full_tableofres <- cbind(forestplot_tableofres, "HCT6vsHCT0")

write.csv(hc_full_tableofres, file.path(this.roc.dir,"mean_zscore_hc_other_published_signatures_tableofauc.csv"))
write.csv(forestplot_full_tableofres, file.path(this.forest.dir,"mean_zscore_hc_other_published_signatures_forestplotres.csv"))

## 5) ROC Curves -----------------------------------------------------------
this.roc.dir <- file.path(hc.res.dir, "roc")
if(!exists(this.roc.dir)) dir.create(this.roc.dir)

# Convert ROC data to a format suitable for ggplot

roc_data <- do.call(rbind, lapply(names(hc_roc_objects), function(comparison) {
  data.frame(
    TPR = rev(hc_roc_objects[[comparison]]$sensitivities),  # True Positive Rate
    FPR = rev(1 - hc_roc_objects[[comparison]]$specificities),  # False Positive Rate
    Signature = comparison,
    auc = rev(hc_roc_objects[[comparison]]$auc)
  )
}))


roc_data$ci <- hc_tableofres[match(roc_data$Signature, hc_full_tableofres$signature), "ci"]

roc_data$legend <- paste0(roc_data$Signature,": \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")


#reorder signature so colours are how i want
roc_data$Signature <- factor(roc_data$Signature, levels = names(gene_sets))

#reorder legend to match
roc_data <- roc_data %>%
  mutate(
    #get first word of signature and remove the :
    legend_first_word = str_remove(word(legend, 1),":"),
    # rorder legend factor based on signature factor order
    legend = factor(legend, levels = legend[match(levels(Signature), legend_first_word)])
  )
hc_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20),
        plot.caption = element_text( size = 12, hjust=0),
        plot.caption.position = "plot" 
) +
  labs(
    title = paste("ROC (Publicly available signatures) - HC_T6 vs HC_T0"),
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Signature",
    caption = paste0("ROC showing performance of a logistic regression model predicting HC_T6 vs HC_TO using the signature score \n",
    "Signature scores calculated as mean of Z-scored expression of signature genes"))

ggsave(hc_roc, filename = file.path(this.roc.dir, "hc_roc.png"), 
       width = 3000, 
       height = 3200, 
       units = "px")



## 6) Forest plot -----------------------------------------------------------

forestplot_full_tableofres

  forestplot_theme <- theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "None") 
  

  # other gene signatures don't need to be plotted, can rename axis to HC_T0 vs TB_T0
  res_table <- forestplot_full_tableofres[-nrow(forestplot_full_tableofres),]
  res_table$signature <- factor(res_table$signature, levels = names(gene_sets))
  
res_table[,-1] <- lapply(res_table[,-1], as.numeric)

  auc_plot <-  res_table %>% 
    ggplot(aes(y = signature)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = signature), shape=15, size=3) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = signature),
                   size = 1) +
    forestplot_theme +
    scale_x_continuous(limits = c(0, 1) )+
    ylab(NULL)+
    xlab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(y = signature)) + 
      theme_bw() +
      geom_point(aes(x=sensitivity, color = signature), shape=15, size=3) +
      geom_linerange(aes(xmin=sensitivity_ci_low, 
                         xmax=sensitivity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
    scale_x_continuous(limits = c(0, 1) )+      
       ylab(NULL)+
      xlab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(y = signature)) + 
      theme_bw() +
      geom_point(aes(x=specificity, color = signature), shape=15, size=3) +
      geom_linerange(aes(xmin=specificity_ci_low, 
                         xmax=specificity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
    scale_x_continuous(limits = c(0, 1) )+ 
      ylab(NULL)+
      xlab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3)
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 24),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating HC_T6 vs HC_T0 \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_z_score_forestplot.png"),
         width = 15, height = 20, units = "cm",   bg = "white"  )
  

#FLIPPED PLOT
  auc_plot <-  res_table %>% 
    ggplot(aes(x = signature)) + 
    theme_bw() +
    geom_point(aes(y=auc, color = signature), shape=15, size=3) +
    geom_linerange(aes(ymin=auc_ci_low, 
                       ymax=auc_ci_high, 
                       color = signature),
                   size = 1) +
    forestplot_theme +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0, 1) )+
    xlab(NULL)+
    ylab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(x = signature)) + 
      theme_bw() +
      geom_point(aes(y=sensitivity, color = signature), shape=15, size=3) +
      geom_linerange(aes(ymin=sensitivity_ci_low, 
                         ymax=sensitivity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
           theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0, 1) )+      
       xlab(NULL)+
      ylab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(x = signature)) + 
      theme_bw() +
      geom_point(aes(y=specificity, color = signature), shape=15, size=3) +
      geom_linerange(aes(ymin=specificity_ci_low, 
                         ymax=specificity_ci_high, 
                         color = signature),
                     size = 1) +
      forestplot_theme +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0, 1) )+ 
      xlab(NULL)+
      ylab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3,
                              heights = c(1,1,1.4))
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 18),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating HC_T6 vs HC_T0 \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_z_score_forestplot_flip.png"),
         width = 12, height = 25, units = "cm",   bg = "white"  )
  

  
  
  
  
  
  
  
#==============================================================================================================#
#============================================== TIMEPOINTS ====================================================
#============================================================================================================== #
# A) FOR BOXPLOT -----------
## 1) Subset (all timepoints) ------
Timepoint <- clinical$Timepoint
expression_auc <- expression #rows are samples instead of genes here

#Subset
all(colnames(expression_auc) == row.names(clinical))
TB_samples <- row.names(clinical)[which(clinical$Disease == "TB")]


expression_auc <- expression_auc[,TB_samples]
expression_auc <- as.matrix(expression_auc)
clinical_auc <- clinical[TB_samples,]


#Get scores for test set-----------------------------------------------------------------------------------------
#Subset for our x number of genes 
gene_set_exp <- t(expression_auc)

# All scores collected
score_list <- lapply(gene_sets, function(set) {
  this_gene_set_exp <- gene_set_exp[, set]
  zscore <- scale(this_gene_set_exp, center = TRUE, scale = TRUE)
  rowMeans(zscore)
})

mean_sig_zscore <- as.data.frame(score_list)

all(row.names(mean_sig_zscore) == row.names(clinical_auc))

## 3) Mean Z-Score BOXPLOT - DISEASE ---------------
meanzscore.res.dir <- file.path(output.othersig.dir, "mean_zscore")
if(!exists(meanzscore.res.dir)) dir.create(meanzscore.res.dir)

timepoint.res.dir <- file.path(meanzscore.res.dir, "timepoint")
if(!exists(timepoint.res.dir)) dir.create(timepoint.res.dir)

this.boxplot.dir <- file.path(timepoint.res.dir, "boxplot")
if(!exists(this.boxplot.dir)) dir.create(this.boxplot.dir)


boxplot_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None",
                    plot.caption = element_text( size = 14, hjust=0),
                    plot.caption.position = "plot")  #align caption to entire plot including margins



boxplot<- as.data.frame(cbind(mean_sig_zscore,
                                    group = clinical_auc$Timepoint))


my_comparisons <- combn(unique(boxplot$group), 2, simplify = FALSE)

boxplot$group <- paste0("TB_", boxplot$group)
x_order <- c("TB_T0", "TB_T2", "TB_T4", "TB_T6")

boxplot$group <- factor(boxplot$group, levels = x_order)

boxplot[,1:(ncol(boxplot)-1)] <- sapply(boxplot[,1:(ncol(boxplot)-1)], as.numeric)

for (i in colnames(boxplot)[1:(ncol(boxplot) - 1)]){
  
  
  stat.table <- boxplot  %>%
    wilcox_test(as.formula(paste(i ,"~ group")),
                paired = FALSE) %>%
    add_xy_position(x = "group")
  
  if(any(stat.table$p < 0.05)){
    stat.table <- stat.table[which(stat.table$p < 0.05),]
    lowest_bracket <- max(boxplot[, i]) + 0.05*(max(boxplot[, i]))
    stat.table$y.position <- lowest_bracket
    stat.table$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table))
    
  } else{
    stat.table <- stat.table[0,]
  }
  
  
  
  boxplotfinal2 <- ggplot(boxplot, aes(
    x = factor(group, level = x_order),
    # x = factor(group),
    y = as.numeric(boxplot[,i]),
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
    
    # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
    # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
    
    theme(axis.text.x = element_text(size = 15))+
    labs(title = paste0("TB Signature Analysis: ", i),
         caption = paste0(str_wrap(paste("Signature:", paste(overlap[[i]], collapse = ", "))), "\n",
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "P < 0.05 from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Disease")
  
  if (nrow(stat.table) > 0) {
    boxplotfinal2 <- boxplotfinal2 + 
      stat_pvalue_manual(stat.table, label = "p", tip.length = 0.01, size = 6)
  }
  
  
  ggsave(boxplotfinal2, filename = file.path(this.boxplot.dir, paste0("meanzscore_plot_", i, ".png")), 
         width = 2300, 
         height = 2500, 
         units = "px" )
  
}



#Combined boxplot
boxplot_long <- pivot_longer(boxplot, cols = colnames(boxplot)[1:(ncol(boxplot) - 1)],
             names_to = "signature",
             values_to = "mean_sig_zscore")


x_order <- names(gene_sets)
boxplot_long$signature <- factor(boxplot_long$signature, levels = x_order)

#Get P-values for WITHIN signature comparisons (HC vs TB within each signature)
stat.table <- boxplot_long  %>%
  group_by(signature) %>% #this is how you want the plot to be grouped ie. the x-axis label
  wilcox_test(mean_sig_zscore ~ group, #this is the variable that will be in the legend
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "signature") #this is the x-axis/the grouping variable


#max y pos for each group
y_pos <- boxplot_long %>%
  group_by(signature) %>%
  summarise(ymax = max(mean_sig_zscore))

  stat.table <- stat.table[which(stat.table$p < 0.05),]
  
stat.table <- stat.table %>%
  left_join(y_pos, by = "signature") %>%
  group_by(signature) %>%
  arrange(signature) %>%  
  mutate(
    y.position = (ymax+ymax*0.1) + 0.5 * (row_number() - 1)
  ) %>%
  ungroup()

    boxplot_long <- as.data.frame(boxplot_long)

  stat.table$signature <- factor(stat.table$signature, levels = x_order)
  
boxplotcombined <- ggplot(boxplot_long, aes(
    x = signature,
    y = mean_sig_zscore,
    color = group)) +
    
    theme_bw()+
    
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 24),
          axis.text.x = element_text(angle = 45, hjust = 1),
          title = element_text(size = 20),
          legend.position = "bottom",
          legend.text = element_text(size = 18),
          plot.caption = element_text( size = 18, hjust=0),
          plot.caption.position = "plot") + 
    
    geom_boxplot(position = position_dodge(0.75),
                 width = 0.4) +
    # 
    # geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
    #            alpha=0.5,
    #           size = 2) +

    stat_summary(fun.y = mean, fill = "red",
                 geom = "point", shape = 21, size =2,
                 position = position_dodge(width = 0.75),
                 show.legend = FALSE) +
    
    stat_pvalue_manual(stat.table,
                       label = "p",
                       tip.length = 0.01,
                       size = 6)+
    
    labs(title = paste0("Performance of Published Signatures"),
         color = "Disease",
         caption = paste0(
         "Signature scores calculated as mean of z-scored expression of signature genes \n",
         "Mean shown as red points on boxplot",
         "P < 0.05 values from Mann-Whitney U test shown")
         )+
    ylab (label = "Signature Score") +
    xlab (label = "Signature")
  
  
  ggsave(boxplotcombined, filename = file.path(this.boxplot.dir, paste0("mean_zscore_plot_combined.png")), 
         width = 7000, 
         height = 3000, 
         units = "px" )
  
  

# B) FOR ROC ----------------------
this.roc.dir <- file.path(timepoint.res.dir, "roc")
if(!exists(this.roc.dir)) dir.create(this.roc.dir)


## 1) Subset (per comparison) ------------
this.forest.dir <- file.path(timepoint.res.dir, "forest")
if(!exists(this.forest.dir)) dir.create(this.forest.dir)

timepoint_full_tableofres <- data.frame()
forestplot_full_tableofres <- data.frame()


for (study in names(gene_sets)){
  print(study)
  this_gene_set <- c(gene_sets[[study]])
  case_and_control <- data.frame(case = c("TB_T2", "TB_T4", "TB_T6"), control = c("TB_T0", "TB_T0", "TB_T0"))
  
  
  timepoint_tableofres <- data.frame()
  timepoint_roc_objects <- list()
  forestplot_tableofres <- data.frame()
  
  for (i in 1:3){
    print(paste("Comparison", i))
    
    case = case_and_control$case[i]
    control = case_and_control$control[i]
    
    TB_timepoint_subset <- row.names(clinical)[which(clinical$condition == case | clinical$condition == control )]
    
    #Subset
    expression_auc <- expression[,TB_timepoint_subset]
    expression_auc <- as.matrix(expression_auc)
    clinical_auc <- clinical[TB_timepoint_subset,]
    
    ## 2) GSVA (per comparison) and ROC --------
    
gene_set_exp <- t(expression_auc)
gene_set_exp <- gene_set_exp[,this_gene_set] #genes as columns
gene_set_zscore<-scale(gene_set_exp, center=T, scale=T)   
mean_sig_zscore<- rowMeans(gene_set_zscore)
mean_sig_zscore <- t(mean_sig_zscore)

    glm_data <- data.frame(Score = mean_sig_zscore[1,], Group = clinical_auc$condition)
    
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
      comparison = paste0(case," vs ",control),
      samples_case = paste(case, "=", sum(glm_data$Group == case)),
      samples_control = paste(control, "=", sum(glm_data$Group == control)),
      auc = auc(roc_obj),
      ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
      sensitivity = optimal_threshold_coords$sensitivity, 
      specificity = optimal_threshold_coords$specificity
    )
    
    
    timepoint_tableofres <- rbind(timepoint_tableofres, res) #contains results for all 4 comparisons for this study
    
    timepoint_roc_objects[[paste0(case," vs ",control)]] <- roc_obj
    
    

    # Sensitivity confidence interval (get the sensitivity at the optimal specificity)
  ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"]))
  
  # Specificity confidence interval (get the specificity at the optimal sensitivity)
  ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))
  
  
  forestplot_res<- cbind(
    comparison = paste0(case," vs ",control),
    auc = auc(roc_obj),
    auc_ci_low = as.numeric(auc_ci[1]),
    auc_ci_high = as.numeric(auc_ci[3]),
    
    sensitivity = optimal_threshold_coords$sensitivity, 
    sensitivity_ci_low = ci_sens[, "2.5%"],
    sensitivity_ci_high = ci_sens[, "97.5%"],
    
    specificity = optimal_threshold_coords$specificity,
    specificity_ci_low = ci_spec[, "2.5%"],
    specificity_ci_high = ci_spec[, "97.5%"])
  
  
  forestplot_tableofres <- rbind(forestplot_tableofres,forestplot_res)
  
    
  } #End of this comparison 
  timepoint_tableofres <- cbind(timepoint_tableofres, study)
  forestplot_tableofres <- cbind(forestplot_tableofres,study)

  
  
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
  
  #reorder signature so colours are how i want
roc_data$Comparison <- factor(roc_data$Comparison, levels = names(timepoint_roc_objects))

#reorder legend to match
roc_data <- roc_data %>%
  mutate(
    #get first word of signature and remove the :
    legend_first_word = str_remove(word(legend, 1, 3),":"),
    # rorder legend factor based on signature factor order
    legend = factor(legend, levels = legend[match(levels(Comparison), legend_first_word)])
  )

  
  timepoint_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
    geom_line(size = 1.2) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
    guides(colour = guide_legend(nrow = 2)) +
    scale_color_manual(values = c("#F564E3","#619CFF", "#B79F00")) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 24),
          legend.text = element_text(size = 16),
          title = element_text(size = 20),
          plot.caption = element_text( size = 14, hjust=0),
                    plot.caption.position = "plot") +
    labs(
      title = paste0("ROC - TB treatment timepoints (",study, ")"),
      x = "FPR (1 - Specificity)",
      y = "TPR(Sensitivity)",
      color = "Comparison",
      caption = paste0(str_wrap(paste("Signature:", paste(overlap[[study]], collapse = " "))),"\n",
                       "ROC showing performance of a logistic regression model predicting timepoint using the signature score \n",
    "Signature scores calculated as mean of Z-scored expression of signature genes")
    )
  
  ggsave(timepoint_roc, filename = file.path(this.roc.dir, paste0(study,"_timepoint_roc.png")), 
         width = 3200, 
         height = 3400, 
         units = "px")
  
  
  timepoint_full_tableofres <- rbind(timepoint_full_tableofres, timepoint_tableofres)
  
  timepoint_full_tableofres <- rbind(timepoint_full_tableofres, timepoint_tableofres)
  forestplot_full_tableofres <- rbind(forestplot_full_tableofres, forestplot_tableofres)

} #end of this signature/study

write.csv(timepoint_full_tableofres, file.path(this.roc.dir,"mean_zscore_timepoint_other_published_signatures_tableofauc.csv"))
write.csv(forestplot_full_tableofres, file.path(this.forest.dir,"mean_zscore_timepoint_other_published_signatures_forestplotres.csv"))

## 6) Forest plot -----------------------------------------------------------

forestplot_full_tableofres



  # other gene signatures don't need to be plotted, can rename axis to HC_T0 vs TB_T0
  res_table <- forestplot_full_tableofres
  res_table$study <- factor(res_table$study, levels = names(gene_sets))
  
res_table[,-c(1, ncol(res_table))] <- lapply(res_table[,-c(1, ncol(res_table))], as.numeric)


    colours<- c(
      "HC_T0 vs TB_T0" = "#F8766D",
      "HC_T0 vs HC_T6" = "#00BA38",
      "TB_T6 vs TB_T0" = "#B79F00",
      "TB_T4 vs TB_T0" = "#619CFF",
      "TB_T2 vs TB_T0" = "#F564E3"
    )
        colours <- colours[3:5]

  auc_plot <-  res_table %>% 
    ggplot(aes(x = study)) + 
    theme_bw() +
    geom_point(aes(y=auc, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +
    geom_linerange(aes(ymin=auc_ci_low, 
                       ymax=auc_ci_high, 
                       color = comparison,
                       ),
                   size = 1,
                   position = position_dodge(0.5)) +
          scale_colour_manual(values = colours)+

      forestplot_theme +
           theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0, 1) )+
    xlab(NULL)+
    ylab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(x = study)) + 
      theme_bw() +
       
       geom_point(aes(y=sensitivity, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +
       geom_linerange(aes(ymin=sensitivity_ci_low, 
                         ymax=sensitivity_ci_high, 
                         color = comparison),
                     size = 1,
                     position = position_dodge(0.5)) +
                 scale_colour_manual(values = colours)+

      forestplot_theme +
       theme(axis.text.x = element_blank(),
             axis.ticks.x = element_blank()) +
    scale_y_continuous(limits = c(0, 1) )+      
       xlab(NULL)+
      ylab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(x = study)) + 
      theme_bw() +
      geom_point(aes(y=specificity, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +      
      geom_linerange(aes(ymin=specificity_ci_low, 
                         ymax=specificity_ci_high, 
                         color = comparison),
                     size = 1,
                     position = position_dodge(0.5)) +
          scale_colour_manual(values = colours)+
  guides(colour = guide_legend(nrow = 2)) +

      theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
                            legend.position = "bottom",
                        legend.text = element_text(size = 14),
                        legend.title  = element_text(size = 0),
            axis.text.x = element_text(angle = 45, hjust = 1
            ))    +  

    scale_y_continuous(limits = c(0, 1) )+ 
      xlab(NULL)+
      ylab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3,
                              heights = c(1,1,1.8))
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 18),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating treatment timepoints \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_zscore_forestplot.png"),
         width = 12, height = 25, units = "cm",   bg = "white"  )
  
  
  
  #Flipped plot
  
  auc_plot <-  res_table %>% 
    ggplot(aes(y = study)) + 
    theme_bw() +
    geom_point(aes(x=auc, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +
    geom_linerange(aes(xmin=auc_ci_low, 
                       xmax=auc_ci_high, 
                       color = comparison,
                       ),
                   size = 1,
                   position = position_dodge(0.5)) +
          scale_colour_manual(values = colours)+

      forestplot_theme +
    scale_x_continuous(limits = c(0, 1) )+
    ylab(NULL)+
    xlab("AUC")
  
     sens_plot <- res_table %>% 
      ggplot(aes(y = study)) + 
      theme_bw() +
       
       geom_point(aes(x=sensitivity, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +
       geom_linerange(aes(xmin=sensitivity_ci_low, 
                         xmax=sensitivity_ci_high, 
                         color = comparison),
                     size = 1,
                     position = position_dodge(0.5)) +
                 scale_colour_manual(values = colours)+

      forestplot_theme +
    scale_x_continuous(limits = c(0, 1) )+      
       ylab(NULL)+
      xlab("Sensitivity")

     
    spec_plot <- res_table %>% 
      ggplot(aes(y = study)) + 
      theme_bw() +
      geom_point(aes(x=specificity, color = comparison), shape=15, size=3,
               position = position_dodge(0.5)) +      
      geom_linerange(aes(xmin=specificity_ci_low, 
                         xmax=specificity_ci_high, 
                         color = comparison),
                     size = 1,
                     position = position_dodge(0.5)) +
          scale_colour_manual(values = colours)+

      theme(axis.title = element_text(size =  16),
                            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
                        legend.title  = element_text(size = 0),

                            legend.position = "bottom")    +  
            # theme(axis.text.y = element_text(angle = 45, hjust = 1))+
    scale_x_continuous(limits = c(0, 1) )+ 
      ylab(NULL)+
      xlab("Specificity")

    
    panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
                              ncol = 1,
                              nrow = 3,
                              heights = c(1,1,1.4))
    
    panel_forest <- annotate_figure(
      panel_forest,
      top = text_grob("Performance of published signatures", size = 18),
      bottom = text_grob(paste0("Forest plot showing AUCs for differentiating treatment timepoints \nin our dataset",
                                "\nSenstivity and specificity calculated at Youden threshold \n",
                                "95% confidence intervals shown"),
                                size = 10, hjust = 0, x = 0)
    )

  ggsave(panel_forest, filename= file.path(this.forest.dir, "mean_zscore_forestplot_flip.png"),
         width = 18, height = 40, units = "cm",   bg = "white"  )
  

