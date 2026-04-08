# SCRIPT 1: Differential Expression
# Raw data files were cleaned up and manipulated + QC checks were conducted prior to this script
# QC checks: "3204316" and "3201323" were identified as outliers in PCA plots and removed
# Samples with missing gene expression data excluded

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #
#Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
main.dir <- file.path(my_directory, "TB_v2")
setwd(file.path(main.dir))
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")

# Mac
my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("/Volumes/One Touch/RBMB/rlibrary")

library(readxl)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(GSVA)


# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
output.dir <- file.path(main.dir, "output_v2")
results.dir <- file.path(output.dir, "1_differential_expression")
if(!exists(results.dir)) dir.create(results.dir)
figures.dir <- file.path(results.dir, "figures")
if(!exists(figures.dir)) dir.create(figures.dir)

#tT and tT2 results dir
if(!exists(file.path(results.dir, "tT"))) dir.create(file.path(results.dir, "tT"))
if(!exists(file.path(results.dir, "tT2"))) dir.create(file.path(results.dir, "tT2"))

#volcano dir
volcano.dir <- file.path(figures.dir, "volcano")
if(!exists(volcano.dir)) dir.create(volcano.dir)

#heatmap dir
heatmap.dir <- file.path(figures.dir, "heatmap")
if(!exists(heatmap.dir)) dir.create(heatmap.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #

clinical <- read.csv(file.path(main.dir, "data/processed/post-QC/clinical.csv"), row.names = 1)
expression <-  read.csv(file.path(main.dir, "data/processed/post-QC/expression.csv"), row.names = 1, check.names = F)

clinical_paired <- read.csv(file.path(main.dir, "data/processed/post-QC/clinical_paired.csv"), row.names = 1)
expression_paired <- read.csv(file.path(main.dir, "data/processed/post-QC/expression_paired.csv"), row.names = 1, check.names = F)

# ================================================================================== #
# 1. FORMAT DATA ==================================================================
# ================================================================================== #
clinical$condition <- as.factor(clinical$condition)
clinical_paired$condition <- as.factor(clinical_paired$condition)


# ================================================================================== #
# 2. LIMMA DIFFERENTIAL EXPRESSION =================================================
# ================================================================================== #

# LIMMA -----------------------------------------------------------------

label <- c(paste0("contrast", seq(1,10, 1)))

contrast <- c(
  "TB_T0 - HC_T0", #Plate 1 and 2: Use all data (100 HC T0 samples vs 100 TB T0 samples)
  "HC_T6 - HC_T0", #Plate 3&4 
  "TB_T6 - TB_T0", #Plate 5-7
  "TB_T6 - HC_T6",
  "TB_T6 - TB_T2",
  "TB_T4 - TB_T0",
  "TB_T2 - TB_T0",
  "TB_T6 - TB_T4",
  "TB_T4 - TB_T2",
  "(TB_T6 - TB_T0) - (HC_T6 - HC_T0)"
)

nameconvert <- as.data.frame(cbind(label = label, contrast = contrast))

nameconvert$c1 <- sub(" -.*", "", nameconvert$contrast)
nameconvert$c2 <- sub(".*- ", "", nameconvert$contrast)

####  LIMMA -----------------
design <- model.matrix(~0 + condition + age + sex,
                       data = clinical)

design_paired <- model.matrix(~0 + condition + age + sex, #lost 43 samples as there is no sex data
                              data = clinical_paired)


colnames(design)[1:length(levels(clinical$condition))] <- levels(clinical$condition)
colnames(design_paired)[1:length(levels(clinical_paired$condition))] <- levels(clinical_paired$condition)

clinical2 <- clinical[row.names(design),]
expression2 <- expression[,row.names(clinical2)]

clinical_paired2 <- clinical_paired[row.names(design_paired),]
expression_paired2 <- expression_paired[,row.names(clinical_paired2)]

#Block by plate for treatment plates
corfit_paired <- duplicateCorrelation(
  expression_paired2,
  design_paired,
  block = clinical_paired2$Plate) 


fit <- lmFit(expression2, 
             design = design)  

fit_paired <- lmFit(expression_paired2, 
                    design = design_paired,
                    block = clinical_paired2$Plate,
                    correlation = corfit_paired$consensus)  

#taking all TB_T0 and HC_T0 and comparing them
cont.matrix <- makeContrasts(
  contrast1 = TB_T0 - HC_T0, 
  levels = design) 

cont.matrix_paired <- makeContrasts(
  contrast2 = HC_T6 - HC_T0,
  contrast3 = TB_T6 - TB_T0,
  contrast4 = TB_T6 - HC_T6,
  contrast5 = TB_T6 - TB_T2,
  contrast6 = TB_T4 - TB_T0,
  contrast7 = TB_T2 - TB_T0,
  contrast8 = TB_T6 - TB_T4,
  contrast9 = TB_T4 - TB_T2,
  #which genes respond differently over time in TB compared to household
  contrast10 = (TB_T6 - TB_T0) - (HC_T6 - HC_T0),
  levels = design_paired)

listoffit <- list(fit = fit, fit_paired = fit_paired)
listofcontmatrix <- list(cont.matrix = cont.matrix, cont.matrix_paired = cont.matrix_paired)

listoftT <- list()
listoftT2 <- list()

for (x in 1:2){
  
  fit <- listoffit[[x]]
  cont.matrix <- listofcontmatrix[[x]]
  
  for (i in colnames(cont.matrix)){
    
    fit2 <- contrasts.fit(fit, contrast=cont.matrix[,i])
    
    fit2 <- eBayes(fit2, robust = FALSE, trend = FALSE)
    
    tT<- topTable(fit2, adjust="BH", sort.by="P", number=nrow(fit2))
    
    tT$genename <- row.names(tT)
    selection <- which(tT$adj.P.Val<0.05)

    tT$legend <-  ifelse(
      tT$adj.P.Val< 0.05 & tT$logFC > 1, "Upregulated",
      ifelse(
        tT$adj.P.Val < 0.05 & tT$logFC < -1, "Downregulated",
        "Not Sig"))
    
    tT$legend[is.na(tT$legend)]="Not Sig"
    
    tT2 <- tT[selection,] #Only includes sites that have significant FDR (Adj p val < 0.05)
    
    listoftT[[i]] <- tT
    listoftT2[[i]] <- tT2
    

    write.csv(tT, file = file.path(results.dir, paste0("tT/", i, ".csv")))
    write.csv(tT2, file = file.path(results.dir, paste0("tT2/", i, ".csv")))
    
  }
}

# ================================================================================== #
## 2.1. VOLCANO PLOTS ==============================================================
# ================================================================================== #

listofvolcano<- list()

for (i in names(listoftT)){
  volcano <- ggplot(listoftT[[i]], aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = legend)) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not Sig" = "grey", "Upregulated" = "red"))+
    geom_hline(yintercept =-log10(max(as.data.frame(listoftT2[[i]][,"P.Value"]))),colour="black", linetype="dashed") +
    geom_vline(xintercept =-1,colour="black", linetype="dashed")+
    geom_vline(xintercept =1,colour="black", linetype="dashed")+
    geom_text_repel(data = as.data.frame(listoftT2[[i]]),
                    aes(label= genename),
                    size = 3, 
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines") 
    ) +
    theme_bw(base_size = 12) +
    
    theme(legend.position = "bottom") +
    
    labs(title = nameconvert[which(nameconvert$label == i), "contrast"]) 
  
  listofvolcano[[i]] <- volcano
  
  #save
  ggsave(volcano, filename = file.path(volcano.dir, paste0(i,".png")), width = 6.96, height = 7.18, units = "in")
  
  
}

# Upregulated genes in TB_T0 - HC_T0
listoftT2[["contrast1"]][which(listoftT2[["contrast1"]]$legend == "Upregulated"), "genename"]

