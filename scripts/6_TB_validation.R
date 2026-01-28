#This script ONLY contains comparison of HC_T0 vs TB_T0 unpaired data
#Tess sent this batch of data first in October 2025, then in December sent the next batch of data (longitudinal data)
#As of 12/12/2025, I have not combined the unpaired HC_T0 and TB_T0 (Batch1) with the paired HC_T0 and TB_T0 samples (Batch2) due to disparity in housekeeping genes
# Batch1 = B2M is more stable
# Batch2 = GAPDH is more stable
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


library(pROC)
library(stats)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(rstatix)
library(ggpubr)
library(readxl)

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory
data.dir <- file.path(my_directory, "data")
processed.dir <- file.path(my_directory, "data", "processed")
output.dir <- file.path(my_directory,"output")
validation.dir <- file.path(output.dir,"validation", "unpaired_HCT0_TBT0")
if(!exists(validation.dir)) dir.create(validation.dir)



# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(data.dir,"raw"))

#clinical
HC_T0_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "HCT0", col_names = TRUE) 
TB_T0_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "TBT0", col_names = TRUE) 

HC_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse with new GAPDH 22012026.xlsx"), sheet = "HCT0", skip = 1, col_names = TRUE)
TB_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse with new GAPDH 22012026.xlsx"), sheet = "TBT0", skip = 1, col_names = TRUE) 
setwd(file.path(main.dir))

# ================================================================================== #
# 2. CLEAN UP EXPRESSION DATA ======================================================
# ================================================================================== #

clinical <- as.data.frame(rbind(HC_T0_meta,
                                TB_T0_meta))


#Check duplicates
clinical[duplicated(clinical$sample),"sample"]

#3201009 is in unpaired TBT0 twice - remove one of them
which(clinical$sample == "3201009") #55 and 62
clinical <- clinical[-which(clinical$sample == "3201009")[2],]

#Check duplicates now again
clinical[duplicated(clinical$sample),"sample"]

clinical$disease_simple <- clinical$disease
clinical$disease <- paste0(clinical$disease_simple, "_", clinical$timepoint)


#Average the duplicated genes

odd_indexes <- seq(from = 1, to= 19,by = 2)

listoffiles <- list(HC_T0 = HC_T0_ct, TB_T0 = TB_T0_ct)


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

###  NORMALISE 
all_data <- cbind(
  HC_T0,
  TB_T0)

colnames(all_data)[duplicated(colnames(all_data))]

#Remove NAs
which(is.na(as.matrix(all_data)), arr.ind = TRUE)

# "3201400" "3205033" Has missing data for at least one gene. Remove
colnames(all_data)[c(13, 72)]
all_data[,c(13, 72)]
all_data <- all_data[,-c(13, 72)]  #NA values for expression


which(is.na(as.matrix(all_data)), arr.ind = TRUE)


#Match up with clinical file
row.names(clinical) <- clinical$sample
dim(clinical) #103
clinical <- clinical[intersect(colnames(all_data), clinical$sample),]
dim(clinical) #101 

all_data <- all_data[,clinical$sample]

write.csv(all_data, file.path(validation.dir, "all_data_beforenorm.csv"))


# ================================================================================== #
# 5. NORMALISE (delta Ct) ==========================================================
# ================================================================================== #


#Remove old GAPDH data
all_data <- all_data[-which(row.names(all_data) == "GAPDH old"),]
row.names(all_data)[which(row.names(all_data) == "GAPDH new")] <- "GAPDH"

  #Calculate relative expression normalized to housekeeping genes
  alldata_genes_ct<- all_data[!(row.names(all_data) == "B2M" |row.names(all_data) == "GAPDH" ),]

  ### 1) Using B2M as the housekeeping gene --------------------------------------------------------------------
  b2m_means <- all_data[row.names(all_data) == "B2M",]
  delta_ct_b2m <- alldata_genes_ct
  # Change to delta ct, which is roughly log2-scaled (since PCR amplification is exponential). every unit difference = 2x change in expression
  for (s in colnames(delta_ct_b2m)) {
    delta_ct_b2m[[s]] <- delta_ct_b2m[[s]] - b2m_means[1, s] 
  }
  #### Transform to linear scale (higher value = higher expression)
  normdata_b2m_hk <- 2^(-delta_ct_b2m)
   
  ### 2) Using GAPDH as the housekeeping gene --------------------------------------------------------------------
  gapdh_means <- all_data[row.names(all_data) == "GAPDH",]
  delta_ct_gapdh <- alldata_genes_ct
  
  for (s in colnames(delta_ct_gapdh)) {
    delta_ct_gapdh[[s]] <- delta_ct_gapdh[[s]] - gapdh_means[1, s]
  }
  #### Transform to linear scale (higher value = higher expression)
  normdata_gapdh_hk <- 2^(-delta_ct_gapdh)
  
  ### 3) Using average of B2M & GAPDH as the housekeeping gene --------------------------------------------------------------------
  twohk_means <- as.data.frame(t(colMeans(all_data[row.names(all_data) == "GAPDH" | row.names(all_data) == "B2M",])))
  delta_ct_twohk <- alldata_genes_ct
  
  for (s in colnames(delta_ct_twohk)) {
    delta_ct_twohk[[s]] <- delta_ct_twohk[[s]] - twohk_means[1, s]
  }
  #### Transform to linear scale (higher value = higher expression)
  normdata_twohk <- 2^(-delta_ct_twohk)
  

write.csv(normdata_b2m_hk, file.path(validation.dir, "normdata_b2m_hk.csv"))
write.csv(normdata_gapdh_hk, file.path(validation.dir, "normdata_gapdh_hk.csv"))
write.csv(normdata_twohk, file.path(validation.dir, "normdata_avg_b2m_gapdh_hk.csv"))


# ================================================================================== #
# PCA  =================================================
# ================================================================================== #
qc.dir <- file.path(validation.dir, "qc")
if(!exists(qc.dir)) dir.create(qc.dir)

listof_delta_ct <- list(B2M = delta_ct_b2m, GAPDH = delta_ct_gapdh, avg_B2M_GAPDH = delta_ct_twohk )

for (i in 1:length(listof_delta_ct)){
  #use delta_ct values, pca expects log-like scaling (ct values are exponential)
  
  delta_ct_data <- listof_delta_ct[[i]]
  pca_res <- prcomp(t(as.data.frame(delta_ct_data)), scale. = TRUE, center = TRUE) #center = TRUE
  
  
  clinical_pca <- cbind(pca_res$x, clinical) # has NA
  
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
# 6.1. BOXPLOTS FOR ALL 2^-delta Ct  =================================================
# ================================================================================== #

listof_normdata <- list(B2M = normdata_b2m_hk, GAPDH = normdata_gapdh_hk, avg_B2M_GAPDH = normdata_twohk )
# 
# #Remove outliers !!!!!!!!
# #When plotting boxplot, this was an extreme GAPDH outlier, so went back and removed it here
outliers <- c("3200865")

listof_normdata <- lapply(
  listof_normdata,
  function(df) df[,-which(colnames(df) %in% outliers)]
)


## 2) Get groups to be compared --------------------------------------------

listofresults <- list()
for (hk in names(listof_normdata)){
  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(figures.dir)) dir.create(figures.dir)
  
  normdata <- listof_normdata[[hk]]
  expression <- as.matrix(normdata)
  clinical <- clinical[colnames(expression),]
  clinical$group <-clinical$disease
  table(clinical$group)

row.names(clinical) <- clinical$sample

clinical$group <-clinical$disease
table(clinical$group)

colnames(expression) == row.names(clinical)

write.csv(expression, file.path(this.output.dir, paste0(hk,"_relative_expression.csv")))
write.csv(clinical, file.path(this.output.dir, "clinical.csv"))



expr_long <- as.data.frame(expression) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to="sample_id", values_to="expression") 


expr_long <- cbind(expr_long, clinical[match(expr_long$sample_id, clinical$sample),])


boxplot_theme <- theme(axis.title = element_text(size = 20),
                       axis.text = element_text(size = 15),
                       title = element_text(size = 20),
                       strip.text.x = element_text(size = 20),
                       legend.text = element_text(size = 15),
                       legend.position = "bottom",
                       axis.text.x = element_blank(),      
                       axis.ticks.x = element_blank()) 

plot <- ggplot(expr_long, aes(x = disease, y = expression, color=disease)) +
  geom_boxplot(outlier.shape=NA, alpha=0.2) +
  geom_jitter(width=0.2, size=1) +
  theme_bw() +
  boxplot_theme +
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =2,
               show.legend = TRUE) +
  stat_compare_means(method = "wilcox.test", 
                     size = 8,
                     comparisons = list(c("HC_T0", "TB_T0")),
                     label.y.npc = 0.9,
                     aes(label = ..p.format..))+
  facet_wrap(~gene, scales="free_y", ncol = 7) +
  
  guides(color = guide_legend(nrow = 1))+
  ylab("Expression (2^-delta Ct)") +
  xlab("Disease") +
  labs(caption = paste("Relative expression normalized to", hk)) +

  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


ggsave(plot, filename= file.path(figures.dir, paste0(hk,"_all_genes_boxplot.png")),
      width = 60, height = 20, units = "cm")


plot_labelled <- plot +
  
  geom_text_repel(
    aes(label = sample),
    size = 4,
    show.legend = FALSE)   

ggsave(plot_labelled, filename= file.path(figures.dir,paste0(hk,"_all_genes_boxplot_labelled.png")),
       width = 60, height = 20, units = "cm")

# ================================================================================== #
# 7) Z-Score Transformation/ SCALED+CENTEREDMEAN  ====================================
# ================================================================================== #


scaledcentered_mean_func <- function(number_of_genes = "7"){
  
expr_set<-expression[row.names(expression),]
# 
# if(number_of_genes == "6"){
#   expr_set <- expr_set[!row.names(expr_set) == "S100A8",]
# }
# 
# if(number_of_genes == "5"){
#   expr_set <- expr_set[-c(which(row.names(expr_set) == "S100A8" | row.names(expr_set) == "CD274")),]
# }
# 
# if(number_of_genes == "4"){
#   expr_set <- expr_set[-c(which(row.names(expr_set) == "S100A8" | row.names(expr_set) == "CD274" | row.names(expr_set) == "IFITM1")),]
# }
#transpose for scaling
expr_set<-t(expr_set)

#Two methods
# 1) center and scale 
# 2) center only 
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
colnames(mean_standardised) <- colnames(expression)

return(list(scores = mean_standardised,
            gene_list = colnames(expr_set)))
}

listofstandardised_scores <- list()
listofstandardised_scores[["7_genes"]] <- scaledcentered_mean_func()
# listofstandardised_scores[["6_genes"]] <- scaledcentered_mean_func(number_of_genes = "6")
# listofstandardised_scores[["5_genes"]] <- scaledcentered_mean_func(number_of_genes = "5")
# listofstandardised_scores[["4_genes"]] <- scaledcentered_mean_func(number_of_genes = "4")

listofresults[[hk]] <- listofstandardised_scores

mean_sigcenter_scores <- do.call(rbind,
                              lapply(names(listofstandardised_scores), 
                                     function(x){
                                       data.frame(sig_scale = listofstandardised_scores[[x]][["scores"]]["sig_scale",],
                                                  numberofgenes = x,
                                                  sample = colnames(listofstandardised_scores[[x]][["scores"]])
                                       )
                                     }
                              )
) %>%  pivot_wider(names_from = "sample",
              values_from = "sig_scale")
write.csv(mean_sigcenter_scores, file.path(this.output.dir, paste0(hk,"_mean_ztransformed_scores.csv" )))
}

saveRDS(listofresults, file.path(validation.dir, "listofresults.RDS"))


# ================================================================================== #
# 7. PLOT BOXPLOTS FOR VALIDATION =================================================
# ================================================================================== #

for (hk in names(listofresults)){
  
listofstandardised_scores <- listofresults[[hk]]
listofboxplots_scaledcentered <- list()
listofboxplots_centered <- list()

this.output.dir <- file.path(validation.dir, hk)
if(!exists(this.output.dir)) dir.create(this.output.dir)

this.figures.dir <- file.path(this.output.dir, "figures")
if(!exists(this.figures.dir)) dir.create(this.figures.dir)

## SCALED & CENTERED =================================================
#loop over 4, 5, 6 and 7genes
# for(g in names(listofstandardised_scores)){
  
g = "7_genes"
  score_data <- listofstandardised_scores[[g]][["scores"]]
  
  gene_list <- listofstandardised_scores[[g]][["gene_list"]]
  
  boxplot_data <- as.data.frame(cbind(score = score_data["sig_scale",],
                                    group = as.character(clinical$disease),
                                    PID = as.character(clinical$sample)))



  
  gsva_theme <- theme(axis.title = element_text(size = 20),
                      axis.text = element_text(size = 20),
                      title = element_text(size = 20),
                      legend.position = "None") 
  


my_comparisons <- list(
  c("HC_T0", "TB_T0") 
)


boxplot_data$score <- as.numeric(boxplot_data$score)
boxplot_data$group <- factor(boxplot_data$group)

stat.table <- boxplot_data  %>%
  wilcox_test(score ~ group,
              paired = FALSE,
              comparisons = my_comparisons) %>%
  add_xy_position(x = "group")


#FOR SCALED+CENTERED
lowest_bracket <- max(boxplot_data$score) + 0.05*(max(boxplot_data$score))
stat.table$y.position <- seq(lowest_bracket, by= 0.3, length.out = nrow(stat.table))

# stat.table <- stat.table[which(stat.table$p <= 0.05),]

boxplot_scaledcentered <- ggplot(boxplot_data, aes(
  x = factor(group),
  # x = factor(group),
  y = as.numeric(boxplot_data[,1]))) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(aes(color = group),position = position_dodge(1)) +
  
  # Unpaired boxplot
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5,
              width = 0.3) +

  #Paired boxplot
  # geom_point(aes(color = group))+
  # geom_line(aes(group = PID), color = "black", alpha = 0.2) +
  
  stat_pvalue_manual(stat.table,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+
  
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  

  labs(title = paste0("HB_T0 vs TB_T0 expression"),
       caption = paste("Signature:", paste0(gene_list, collapse = ","))) +
  ylab (label = "Mean of scaled & centered expression") +
  xlab (label = "Condition") 


ggsave(boxplot_scaledcentered, filename = file.path(this.figures.dir, paste0(hk,"_boxplot_mean_ztransformed.png")),
       width = 2000,
       height = 2200,
       units = "px" )


# listofboxplots_scaledcentered[[g]] <- ggplotGrob(boxplot_scaledcentered) #ggGrob freezes the image in place, otherwise the pvalue brackets move when put into the list


# } #close loop for number of genes
# 
# boxplot_scaledcentered_panel <- annotate_figure(
#   ggarrange(
#     plotlist = list(listofboxplots_scaledcentered[["7_genes"]],
#                     listofboxplots_scaledcentered[["6_genes"]],
#                     listofboxplots_scaledcentered[["5_genes"]],
#                     listofboxplots_scaledcentered[["4_genes"]]),
#     ncol = 4),
#   bottom = text_grob(paste("Relative expression normalized to", hk),
#                      hjust = 1, x = 1, size = 16))
# 
# ggsave(boxplot_scaledcentered_panel, filename = file.path(this.figures.dir, paste0(hk,"_boxplot_mean_ztransformed_panel.png")),
#        width = 8000,
#        height = 2000,
#        units = "px" )


} #close loop for hk gene



# ================================================================================== #
# 8. PLOT ROC/AUC curves =================================================
# ================================================================================== #

for (hk in names(listofresults)){ 

  listofstandardised_scores <- listofresults[[hk]]
  listofboxplots_scaledcentered <- list()
  listofboxplots_centered <- list()

  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  this.figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(this.figures.dir)) dir.create(this.figures.dir)
  
  res_table <- data.frame() 
  forestplot_res_table <- data.frame()
  roc_objects <- list()
  
  
  # for(g in names(listofstandardised_scores)){ 

g = "7_genes"

  score_data <- listofstandardised_scores[[g]][["scores"]]
  gene_list <- paste(
    "7_genes:", paste0(listofstandardised_scores[["7_genes"]][["gene_list"]], collapse = ","), "\n" #,
    # "6_genes:", paste0(listofstandardised_scores[["6_genes"]][["gene_list"]], collapse = ","), "\n",
    # "5_genes:", paste0(listofstandardised_scores[["5_genes"]][["gene_list"]], collapse = ","), "\n",
    # "4_genes:", paste0(listofstandardised_scores[["4_genes"]][["gene_list"]], collapse = ","), "\n"
    
  )

  
  
mean_standardised_scores <- as.data.frame(t(as.matrix(score_data)))
mean_standardised_scores$group <- clinical[colnames(score_data), "group"]

mean_standardised_scores$group <- factor(mean_standardised_scores$group, levels = c("HC_T0", "TB_T0"))


#scaled+centered
glm_model <- glm(group ~ sig_scale, data = mean_standardised_scores, family = binomial)

# Predict probabilities
test_probs <- predict(glm_model, type = "response")

# Compute ROC and AUC
roc_obj <- roc(mean_standardised_scores$group, test_probs)

auc_ci <- ci.auc(roc_obj) # # 2.5% and 97.5% 

optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))

if(nrow(optimal_threshold_coords) > 1) {
  optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
}


# Sensitivity confidence interval
ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"]))

# Specificity confidence interval
ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))


res_current <-cbind(
  numberofgenes = g,
  auc = auc(roc_obj),
  ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
  sensitivity = optimal_threshold_coords$sensitivity, 
  specificity = optimal_threshold_coords$specificity)


res_table <- rbind(res_table, res_current)
write.csv(res_table, file.path(this.output.dir, paste0(hk,"_mean_ztransformed_scores_res_table.csv")))


forestplot_res_table <- rbind(forestplot_res_table, 
                              cbind(numberofgenes = g,
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

write.csv(forestplot_res_table, file.path(this.output.dir, paste0(hk,"_mean_ztransformed_scores_forestplot_res_table.csv")))

roc_objects[[g]] <- roc_obj


#do.call to execute a function, using a list of arguments (here we are saying, apply the function to all roc_objects names, then rbind the results together )
roc_data <- do.call(rbind, 
                    lapply(names(roc_objects), function(numberofgenes) {
                      data.frame(
                        TPR = rev(roc_objects[[numberofgenes]]$sensitivities),  # True Positive Rate
                        FPR = rev(1 - roc_objects[[numberofgenes]]$specificities),  # False Positive Rat
                        numberofgenes = numberofgenes,
                        auc = rev(roc_objects[[numberofgenes]]$auc))
                      }) #close lapply
                    ) #close do.call


roc_data$numberofgenes <- factor(roc_data$numberofgenes)

roc_data$ci <- res_table[match(roc_data$numberofgenes, res_table$numberofgenes), "ci"]

#When we want to show 4,5,6 and 7 gene signature
# roc_data$legend <- paste0(roc_data$numberofgenes,": \n AUC = ", 
#                           round(roc_data$auc, 2), " (", roc_data$ci, ")")

# When we want to show only 7 gene signature
roc_data$legend <- paste0("HC_T0 vs TB_T0",": \n AUC = ",
                          round(roc_data$auc, 2), " (", roc_data$ci, ")")



disease_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - HC_T0 vs TB_T0",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = paste(g, "Signature:", paste0(listofstandardised_scores[[g]][["gene_list"]], collapse = ",")))

ggsave(disease_roc, filename = file.path(this.figures.dir, paste0(hk, "_ROC_mean_z_transformed_scores.png")),
       width = 2500,
       height = 3000,
       units = "px" )


  # } #close number of genes loop
  
} 
#scaled + centered is better



### Forest plots ####
library(tidyverse)


for (hk in names(listofresults)){ # housekeeping gene loop
  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  this.figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(this.figures.dir)) dir.create(this.figures.dir)
  
  listofstandardised_scores <- listofresults[[hk]]
  
res_table <- read.csv(file.path(this.output.dir, paste0(hk,"_mean_ztransformed_scores_forestplot_res_table.csv")), row.names = 1)

# other gene signatures don't need to be plotted, can rename axis to HC_T0 vs TB_T0
res_table$numberofgenes <- "HC_T0 vs TB_T0"
auc_plot <- res_table %>% 
  ggplot(aes(y = numberofgenes)) + 
  theme_bw() +
  geom_point(aes(x=auc), shape=15, size=3) +
  geom_linerange(aes(xmin=auc_ci_low, xmax=auc_ci_high)) +
  ylab(NULL)+
  xlab("AUC")+
coord_cartesian( xlim=c(0.1, 1))

sens_plot <- res_table %>% 
  ggplot(aes(y = numberofgenes)) + 
  theme_bw() +
  geom_point(aes(x=sensitivity), shape=15, size=3) +
  geom_linerange(aes(xmin=sensitivity_ci_low, xmax=sensitivity_ci_high)) +
  ylab(NULL)+
  xlab("Specificity")+
coord_cartesian( xlim=c(0.1, 1))

spec_plot <- res_table %>% 
  ggplot(aes(y = numberofgenes)) + 
  theme_bw() +
  geom_point(aes(x=specificity), shape=15, size=3) +
  geom_linerange(aes(xmin=specificity_ci_low, xmax=specificity_ci_high)) +
  ylab(NULL)+
  xlab("Sensitivity")+
  coord_cartesian( xlim=c(0.1, 1))
  

panel_forest <- ggarrange(plotlist = list(auc_plot, sens_plot, spec_plot ),
          ncol = 1,
          nrow = 3)

panel_forest <- annotate_figure(
  panel_forest,
  bottom = text_grob(paste0("Normalised to ",hk,"\nSenstivity and specificity calculated at Youden threshold \n",
                            gene_list), 
                      size = 8, hjust = 0, x = 0)
)


ggsave(panel_forest, filename= file.path(this.figures.dir, paste0(hk, "_forestplot_panel.png")),
       width = 10, height = 15, units = "cm",   bg = "white"  )

}

