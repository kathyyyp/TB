
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
    # disregard NAs
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
# 3. REMOVE 18S AND SUBTRACT FROM 40 ===============================================
# ================================================================================== #


for (i in 1:length(listoffiles2)){
    df <- listoffiles2[[i]]
    df_numeric <- as.data.frame(apply(df, 2, as.numeric))
    row.names(df_numeric) <- row.names(df)
    df_numeric <- as.matrix(df_numeric)
    df_numeric[df_numeric == 999] <- 40
    df_numeric <- 40 - df_numeric
    listoffiles2[[i]] <- df_numeric
}



  #Create R objects from the list
  for(i in 1:length(listoffiles2)){
    assign(names(listoffiles2)[i], listoffiles2[[i]])
  }

  #Check if there are any NAs
# YES NAs
  any(is.na(listoffiles2[[1]]))
  any(is.na(listoffiles2[[2]]))
  
  
  
  
  # ================================================================================== #
  # 4. COMBINE ALL DATA ==============================================================
  # ================================================================================== #
  listofdata_beforenorm <- c(listoffiles2)
  
  ##Reorder to ensure all matrices are in same order before binding them  ---------------------------------------------
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
  
  # 
  # # ================================================================================== #
  # # 5. NORMALISE (CYCLICLOESS) =======================================================
  # # ================================================================================== #
  
  # # Not great method for when we only have 7 genes 
  
  # It adjusts each sample’s intensity profile so that systematic biases (like non-linear intensity effects between samples) are removed.
  # It does this iteratively (“cyclic”) — normalizing pairs of samples repeatedly until the system stabilizes
  # The method assumes that most genes should behave similarly across samples (i.e., not all genes are differentially expressed)
  # Upweight B2M and GAPDH - essentially saying that these genes are reliable and to trust them more heavily when fitting the loess curves that define how samples should be normalized.
    # Ensures their values remain almost identical across all samples (they become the “baseline” for correction)
  #set weight of housekeeping genes as 50 and other genes as 1
  weights = rep(1, times = nrow(all_data))
  weights[row.names(all_data) == "GAPDH"] <- 50
  weights[row.names(all_data) == "B2M"] <- 50


  norm_data <- normalizeCyclicLoess(all_data, weights = weights) #other arguments, = span, iterations, method
  
  all_data_onlyb2m <- all_data[-which(row.names(all_data)=="GAPDH"),]
  weights_onlyb2m = rep(1, times = nrow(all_data_onlyb2m))
  weights_onlyb2m[row.names(all_data_onlyb2m) == "B2M"] <- 50
  
  all_data_onlygapdh <- all_data[-which(row.names(all_data)=="B2M"),]
  weights_onlygapdh = rep(1, times = nrow(all_data_onlygapdh))
  weights_onlygapdh[row.names(all_data_onlygapdh) == "GAPDH"] <- 50
  
  norm_data_onlyb2m <- normalizeCyclicLoess(all_data_onlyb2m, weights = weights_onlyb2m) #other arguments, = span, iterations, method
  norm_data_onlygapdh <- normalizeCyclicLoess(all_data_onlygapdh, weights = weights_onlygapdh) #other arguments, = span, iterations, method
  
write.csv(norm_data_onlyb2m, file.path(output.dir, "norm_data_onlyb2m.csv"))
write.csv(norm_data_onlygapdh, file.path(output.dir, "norm_data_onlygapdh.csv"))

  # pca_res <- prcomp(t(as.data.frame(norm_data[-which(row.names(norm_data) %in% c("B2M", "GAPDH")),])), scale. = TRUE, center = TRUE) #center = TRUE

  # variance before normalisation
  apply(all_data, 1, var, na.rm = TRUE) #1 means rows
  apply(norm_data, 1, function(x) var(x, na.rm = TRUE)) #1 means rows
  
  # variance after normalisation
  # notice that b2m and gapdh have 0 variance - they were used as anchor for normalising
  # S100A8 also has 0 variance, already not much variance in original data and was therefore flattened by loess curve
  apply(norm_data_onlyb2m, 1, function(x) var(x, na.rm = TRUE)) #1 means rows
  apply(norm_data_onlygapdh, 1, function(x) var(x, na.rm = TRUE)) #1 means rows
  

