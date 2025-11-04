
# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

# #Mac
my_directory <- "/Users/kathyphung/Library/CloudStorage/OneDrive-UTS/Documents/RBMB"
setwd(paste0(my_directory,"/TB"))
.libPaths("/Volumes/One Touch/rlibrary")

# Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB"
setwd(file.path(my_directory))
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")


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
library(readxl)
library(limma)

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory
data.dir <- file.path(my_directory, "data")
processed.dir <- file.path(my_directory, "data", "processed")
output.dir <- file.path(my_directory,"output")
validation.dir <- file.path(output.dir,"validation")
if(!exists(validation.dir)) dir.create(validation.dir)

results.dir <- file.path(output.dir, "results")
signature.dir <- file.path(results.dir, "signature")


# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(data.dir,"raw"))

#clinical
clinical <- read.csv(file.path(processed.dir, "pre-QC/clinical.csv"), row.names = 1)
HC_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Data for Analyse v1 17102025.xlsx"), sheet = "HCT0", skip = 1, col_names = TRUE)
TB_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Data for Analyse v1 17102025.xlsx"), sheet = "TBT0", skip = 1, col_names = TRUE)

setwd(file.path(main.dir))


# ================================================================================== #
# 2. CLEAN UP EXPRESSION DATA ======================================================
# ================================================================================== #
listoffiles <- list(HC_T0 = HC_T0_ct, TB_T0 = TB_T0_ct)



odd_indexes <- seq(from = 1, to= 17,by = 2)


listoffiles2 <- lapply(X = listoffiles, function(list_item){
  
  new_list <- data.frame() #make new list
  for (i in odd_indexes){
    print(i)
    
    # take mean of the duplicates f, disregard column with the gene name
    # disregard NAs (if there are NA, the remaining duplicate is used eg. alldata_genes_ct["CD274","3205046"] is 38.01081)
    colmeans <- colMeans(list_item[c(i, i +1), -1, drop = FALSE], na.rm = TRUE) 
    gene_data <- c(list_item[i,1], colmeans) #recombine gene name and averaged ct values
    gene_df <- data.frame(gene_data, check.names = FALSE) #turn into a data frame
    new_list <- rbind(new_list, gene_df) #bind this gene to new_list
  } #close for loop
  
  colnames(new_list) <- colnames(list_item)
  rownames(new_list) <- new_list[,1] #make gene names the rownames
  new_list <- new_list[,-1]
  return(new_list)
  
} #close lapply function
) #close lapply

  
  
  
  # ================================================================================== #
  # 4. COMBINE ALL DATA ==============================================================
  # ================================================================================== #
  listofdata_beforenorm <- c(listoffiles2)
  
  ##Reorder to ensure all matrices are in same order before binding them  ---------------------------------------------
  #loop incase there are more files in future
  geneorder <- row.names(listoffiles2[[1]])
  for (i in 1:length(listofdata_beforenorm)){
    df <- listofdata_beforenorm[[i]]
    df <- df[geneorder,]
    
    listofdata_beforenorm[[i]] <- df
  }
  
  #Create R objects from the list
  for(i in 1:length(listofdata_beforenorm)){
    assign(names(listofdata_beforenorm)[i], listofdata_beforenorm[[i]])
  }
  

  all_data <- cbind(HC_T0, TB_T0)
  

# ================================================================================== #
# 5. NORMALISE (delta Ct) ==========================================================
# ================================================================================== #
  
  alldata_genes_ct<- all_data[!(row.names(all_data) == "B2M" |row.names(all_data) == "GAPDH" ),]
  alldata_genes_ct <- alldata_genes_ct[,-which(colnames(alldata_genes_ct) == "3201400")]
  
  ### Using B2M as the housekeeping gene --------------------------------------------------------------------
  b2m_means <- all_data[row.names(all_data) == "B2M",]
  delta_ct_b2m <- alldata_genes_ct
  # Change to delta ct, which is roughly log2-scaled (since PCR amplification is exponential). every unit difference = 2x change in expression
  for (s in colnames(delta_ct_b2m)) {
    delta_ct_b2m[[s]] <- delta_ct_b2m[[s]] - b2m_means[1, s] 
  }
  #### Transform to linear scale (higher value = higher expression)
  normdata_b2m_hk <- 2^(-delta_ct_b2m)
   
  ### Using GAPDH as the housekeeping gene --------------------------------------------------------------------
  gapdh_means <- all_data[row.names(all_data) == "GAPDH",]
  delta_ct_gapdh <- alldata_genes_ct
  
  for (s in colnames(delta_ct_gapdh)) {
    delta_ct_gapdh[[s]] <- delta_ct_gapdh[[s]] - gapdh_means[1, s]
  }
  #### Transform to linear scale (higher value = higher expression)
  normdata_gapdh_hk <- 2^(-delta_ct_gapdh)
  

write.csv(normdata_b2m_hk, file.path(output.dir, "normdata_b2m_hk.csv"))
write.csv(normdata_gapdh_hk, file.path(output.dir, "normdata_gapdh_hk.csv"))
write.csv(all_data, file.path(output.dir, "all_data_beforenorm.csv"))

which(is.na(as.matrix(delta_ct_b2m)), arr.ind = TRUE)
#3205033 has NA for CD274 as there were no values for either duplicate


qc.dir <- file.path(validation.dir, "qc")
if(!exists(qc.dir)) dir.create(qc.dir)

listof_delta_ct <- list(B2M = delta_ct_b2m, GAPDH = delta_ct_gapdh)

for (i in 1:length(listof_delta_ct)){
#use delta_ct values, pca expects log-like scaling (ct values are exponential)
  
delta_ct_data <- listof_delta_ct[[i]]
pca_res <- prcomp(t(as.data.frame(delta_ct_data[,-which(colnames(delta_ct_data) == "3205033")])), scale. = TRUE, center = TRUE) #center = TRUE

colnames(delta_ct_data) %in% colnames(HC_T0)
colnames(delta_ct_data) %in% colnames(TB_T0)

clinical <- data.frame(sample = colnames(delta_ct_data))
clinical$disease <- ifelse(clinical$sample %in% colnames(HC_T0),
                           "HC_T0",
                           "TB_T0")

row.names(clinical) <- clinical$sample
clinical_pca <- cbind(pca_res$x, clinical[-which(row.names(clinical) == "3205033"),])

png(filename = file.path(qc.dir, paste0(names(listof_delta_ct)[[i]],"_pca5_plot.png")), width = 26, height = 23, units = "cm", res = 1200)

variable <- as.factor(clinical_pca$disease)

pcaplot <- pairs(clinical_pca[,c(1:5)],
                 main = paste("PCA (with delta Ct values and", names(listof_delta_ct)[[i]],"as HK gene)"),
                 pch = 19,
                 cex = 0.75,
                 oma=c(3,3,6,20),
                 col = variable)

par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c(levels(variable)), title = "Disease")
dev.off()
}

# ================================================================================== #
# 6. Z-Score and Mean Signature performance  =========================================================================
# ================================================================================== #
# use 2^-delta ct values (linear)  (ie. very right skewed distribution, most values are around zero and some are higher)

figures.dir <- file.path(validation.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)

## 2) Get groups to be compared --------------------------------------------
normdata <- normdata_b2m_hk
if(any(is.na(clinical$disease))){
  clinical <- clinical[-which(is.na(clinical$disease)),]
}

clinical$group <-clinical$disease
table(clinical$group)


expression <- as.matrix(normdata)

colnames(expression) == row.names(clinical)

# write.table(expression, file.path("expression.txt"))
# write.table(clinical, file.path("clinical.txt"))

setwd(file.path(my_directory,"TB", "data", "public", this.accession.no))

# shiny.data.dir <- file.path(my_directory,"TB", "shiny", "data", "public", this.accession.no)
# if(!exists(shiny.data.dir)) dir.create(shiny.data.dir, recursive = TRUE)
# write.table(expression, file.path(shiny.data.dir, "expression.txt"))
# write.table(clinical, file.path(shiny.data.dir, "clinical.txt"))


clinical <- clinical[-which(row.names(clinical) == "3205033"),]
expression <- expression[,-which(colnames(expression) == "3205033")]


# ================================================================================== #
# 6.1. BOXPLOTS FOR ALL 2^-delta Ct  =================================================
# ================================================================================== #
expr_long <- as.data.frame(expression) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample_id", values_to="expression") 


expr_long <- cbind(expr_long, clinical[match(expr_long$sample_id, clinical$sample),])

png(filename = file.path(figures.dir, "all_genes_boxplot.png"),
    width = 15, height = 15, units = "cm", res = 1200)

ggplot(expr_long, aes(x = disease, y = expression, color=disease)) +
  geom_boxplot(outlier.shape=NA, alpha=0.2) +
  geom_jitter(width=0.2, size=1) +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =2,
               show.legend = TRUE) +
  stat_compare_means(method = "wilcox.test", 
                     size = 3,
                     label.y.npc = 0.9,
                     aes(label = ..p.format..))+
  facet_wrap(~gene, scales="free_y") +
  ylab("Expression (2^-delta Ct)") +
  xlab("Disease") +
  theme_minimal() +
  theme(legend.position = "none")


dev.off()


# ================================================================================== #
# 6A) BOXPLOTS FOR Z-SCORE and MEAN  =================================================
# ================================================================================== #
# option A: mean linear expression (2^-delta Ct)
sig_score_linear <- colMeans(expression)

# Option B: z-score per gene before averaging (to give equal weight)
expr_z <- t(scale(t(expr_linear)))  # z-score each gene across samples
sig_score_z <- colMeans(expr_z, na.rm = TRUE)

# Combine and test ---
results <- data.frame(
  sample = clinical$sample,
  disease = clinical$disease,
  sigscore_linear = sig_score_linear[clinical$sample],
  sigscore_z = sig_score_z[clinical$sample]
)

# Boxplot of mean 2^-ΔCt scores
ggplot(results, aes(x = disease, y = sigscore_linear, color = disease)) +
  geom_jitter(width = 0.15, size = 2) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  theme_minimal(base_size = 13) +
  labs(
    title = "7-gene TB signature validation",
    y = "Mean 2^-delta Ct (signature score)",
    x = NULL
  )

# Wilcoxon test
wilcox.test(sigscore_linear ~ disease, data = results)


png(filename = file.path(figures.dir, "auc_plot.png"),
    width = 15, height = 15, units = "cm", res = 1200)

## ROC
roc_obj <- roc(results$disease, results$sigscore_linear, levels = c("HC_T0", "TB_T0"), direction = ">") #direction > means that we expect higher sigscore_linear to indiciate the poisitve class (TB_T0)
auc(roc_obj)
plot(roc_obj, print.auc = TRUE, print.auc.y = 0.8, col = "#2E86AB", main = "ROC curve - TB signature")

# Plot with z-scored version
roc_obj_z <- roc(results$disease, results$sigscore_z, levels = c("HC_T0", "TB_T0"), direction = ">")
auc(roc_obj_z)
plot(roc_obj_z, print.auc = TRUE,  add = TRUE, col = "darkorange")
legend("bottomright", legend = c("Linear mean", "Z-score mean"),
       col = c("#2E86AB", "darkorange"), lwd = 2, bty = "n")

dev.off()




# ================================================================================== #
# 7) SCALED MEAN  =================================================
# ================================================================================== #

###Insert each of the signatures where "il17sig1" is and run this code for each of the signatures. 
# il17ind1<-match(il17sig1,genes)

row.names(expression)
expr_set<-expression[row.names(expression),]

#transpose for scaling
expr_set<-t(expr_set)

summary(expr_set)
boxplot(expr_set)
hist(expr_set)

#center and scale for one score, center only for the other
#centers and/or scales the data based on its mean and standard deviation (so that it has a mean of 0 and a standard deviation of 1)
#If center = TRUE, each value is adjusted by subtracting the mean of the dataset
#If scale = TRUE, each value is divided by the standard deviation after centering
# If center is TRUE then centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns

expr_scale<-scale(expr_set, center=T, scale=T)   # This results in a standardized dataset with mean = 0 and standard deviation = 1 (z-score transformation).
expr_center<-scale(expr_set, center=T, scale=F) #This shifts the data so that its new mean is 0, but the values retain their original scale.


#Transpose back to original position for heatmap
expr_scale_t <-t(expr_scale)
expr_center_t <-t(expr_center)


#Create scores
#take the mean of the scaled+centred data for each gene
sig_scale<-c()
for(i in 1:nrow(expr_scale)) {
  sig_scale[i]<-mean(expr_scale[i,])
}

#take the mean of the centred data for each gene
sig_center<-c()
for(i in 1:nrow(expr_center)) {
  sig_center[i]<-mean(expr_center[i,])
}




mean_standardised<-rbind(sig_scale,sig_center)





# ================================================================================== #
# 7. PLOT BOXPLOTS FOR VALIDATION =================================================
# ================================================================================== #

## SCALED & CENTERED =================================================
boxplot_data <- as.data.frame(cbind(score = mean_standardised["sig_scale",],
                                    group = as.character(clinical$disease),
                                    PID = as.character(clinical$sample)))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- list(
  c("HC_T0", "TB_T0") #,
  # c("TB_T0", "TB_T2"),
  # c("TB_T0", "TB_T4"),
  # c("TB_T0", "TB_T6")
)


boxplot_data$score <- as.numeric(boxplot_data$score)
boxplot_data$group <- factor(boxplot_data$group)

stat.table <- boxplot_data  %>%
  wilcox_test(score ~ group,
              paired = FALSE,
              comparisons = my_comparisons) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]

#FOR SCALED+CENTERED
lowest_bracket <- max(boxplot_data$score) + 0.05*(max(boxplot_data$score))
stat.table$y.position <- seq(lowest_bracket, by= 0.3, length.out = nrow(stat.table))




boxplotfinal2 <- ggplot(boxplot_data, aes(
  x = factor(group),
  # x = factor(group),
  y = as.numeric(boxplot_data[,1]))) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(aes(color = group),position = position_dodge(1)) +
  
  # Unpaired boxplot
  # geom_jitter(aes(color = group),
  #             alpha = 0.5,
  #             size = 2.5,
  #             width = 0.3) +
  
  #Paired boxplot
  geom_point(aes(color = group))+
  geom_line(aes(group = PID), color = "black", alpha = 0.2) +
  
  stat_pvalue_manual(stat.table,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+
  
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  labs(title = paste0("Mean of scaled & centered expression"),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Mean of scaled & centered expression") +
  xlab (label = "Condition")

ggsave(boxplotfinal2, filename = file.path(figures.dir, "scaledcentered_validation_all_paired.png"), 
       width = 3500, 
       height = 3200, 
       units = "px" )




## CENTERED =================================================
boxplot_data <- as.data.frame(cbind(score = mean_standardised["sig_center",],
                                    group = as.character(clinical$disease),
                                    PID = as.character(clinical$sample)))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- list(
  # c("HC_T0", "HC_T6"),
  c("HC_T0", "TB_T0")#,
  # c("TB_T0", "TB_T2"),
  # c("TB_T0", "TB_T4"),
  # c("TB_T0", "TB_T6")
)


boxplot_data$score <- as.numeric(boxplot_data$score)
boxplot_data$group <- factor(boxplot_data$group)

stat.table <- boxplot_data  %>%
  wilcox_test(score ~ group,
              paired = FALSE,
              comparisons = my_comparisons) %>%
  add_xy_position(x = "group")

stat.table <- stat.table[which(stat.table$p < 0.05),]

# FOR CENTERED
lowest_bracket <- max(boxplot_data$score)
stat.table$y.position <- seq(lowest_bracket, by= 0.8, length.out = nrow(stat.table)) #centered


boxplotfinal2 <- ggplot(boxplot_data, aes(
  x = factor(group),
  # x = factor(group),
  y = as.numeric(boxplot_data[,1]))) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(aes(color = group),position = position_dodge(1)) +
  
  # Unpaired boxplot
  # geom_jitter(aes(color = group),
  #             alpha = 0.5,
  #             size = 2.5,
  #             width = 0.3) +
  
  #Paired boxplot
  geom_point(aes(color = group))+
  geom_line(aes(group = PID), color = "black", alpha = 0.2) +
  
  stat_pvalue_manual(stat.table,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+
  
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
  # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
  
  labs(title = paste0("Mean of centered expression"),
       caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") +
  ylab (label = "Mean of centered expression") +
  xlab (label = "Condition")

ggsave(boxplotfinal2, filename = file.path(figures.dir, "centered_validation_all_paired.png"), 
       width = 3500, 
       height = 3200, 
       units = "px" )



# ================================================================================== #
# 8. PLOT BOXPLOTS FOR ROC/AUC curves =================================================
# ================================================================================== #

# GSVA gives a score per sample - enrichment score of how enriched the gene set is in that sample (how up or down the gene set is in the sample relative to other genes)
# Since we can't do GSVA here (our new data has the same 7 genes), we can take our mean scaled and centered data (we already calculated the mean expression of the new 7 genes per sample, scaled and centered)
# Mean of expression (scaled & centered) = the values are comparable bc they reflect how high or low expression of the genes are in that sample compared to other samples

genes <- row.names(expression)
mean_standardised<-rbind(sig_scale,sig_center)
colnames(mean_standardised) <- colnames(expression)

mean_standardised_scores <- as.data.frame(t(mean_standardised))
mean_standardised_scores$group <- clinical[colnames(mean_standardised), "group"]

mean_standardised_scores$group <- factor(mean_standardised_scores$group, levels = c("HC_T0", "TB_T0"))

#scaled+centered
glm_model <- glm(group ~ sig_scale, data = mean_standardised_scores, family = binomial)

# Predict probabilities
test_probs <- predict(glm_model, type = "response")

# Compute ROC and AUC
roc_obj <- roc(mean_standardised_scores$group, test_probs)
plot(roc_obj, print.auc = TRUE)
auc(roc_obj)
auc_ci <- ci.auc(roc_obj)


roc_data <- data.frame(
  TPR = rev(roc_obj$sensitivities),  # True Positive Rate
  FPR = rev(1 - roc_obj$specificities),  # False Positive Rate
  # Comparison = comparison,
  auc = rev(roc_obj$auc),
  ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)))

roc_data$legend <- paste0("HC_T0 vs TB_T0: \n AUC = ", 
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")

disease_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - HC_T0 vs TB_T0",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 

ggsave(disease_roc, filename = file.path(figures.dir, "scaledcentered_roc.png"), 
       width = 3200, 
       height = 3500, 
       units = "px")
