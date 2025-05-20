# ================================================================================== #
# OPTIONAL. LOAD WORKSPACES ========================================================
# ================================================================================== #



# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

#run scripts 1 and 2

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #

output.dir <- file.path(my_directory,"TB", "output")
results.dir <- file.path(output.dir, "results")
figures.dir <- file.path(output.dir, "figures")

if(!exists(output.dir)) dir.create(output.dir)
if(!exists(results.dir)) dir.create(results.dir)
if(!exists(figures.dir)) dir.create(figures.dir)

# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #

#CONTINUE ON FROM SCRIPT 1_TB_DATAWRANGLING.r

# ================================================================================== #
# 2. QUALITY CONTROL ==================================================================
# ================================================================================== #

## PCA Plots-----
exp <- expression
clin <- clinical

#SET FOLDER
filt <- "post-filt"
pca_res <- prcomp(t(as.data.frame(exp)), scale. = TRUE, center = TRUE) #center = TRUE

pca.dir <- file.path(figures.dir, "pca")

if(!exists(pca.dir))dir.create(pca.dir, recursive = TRUE)
if(!exists(file.path(pca.dir, filt))) dir.create(file.path(pca.dir, filt), recursive = TRUE)

# PCs <- pca_res[["x"]]

eigenvalues <- pca_res$sdev^2 #eigenvalues = stdev squared
var_explained = pca_res$sdev^2 / sum(pca_res$sdev^2)

#PCA plot
plot(var_explained[1:20], type = "b",
     xlab = "Principal Component",
     ylab = "Percentage of Variance Explained")


ggsave(autoplot(pca_res,data=clin, colour = "Plate", label = TRUE) + ggtitle("Plate"),
       file = file.path(pca.dir,filt,paste0("outliers",".png")),
       width = 1900,
       height = 1500,
       units = "px")


for (i in c("Plate", "Disease", "Timepoint", "condition", "age", "sex", "smoking_status")){
  ggsave(autoplot(pca_res,data=clin, colour = i) + ggtitle(i),
         file = file.path(pca.dir,filt,paste0(i,".png")),
         width = 1900,
         height = 1500,
         units = "px")
}


ggsave(autoplot(pca_res,data=clin, colour = "PID") + guides(col = "none") + ggtitle("PID"),
       file = file.path(pca.dir,filt,paste0("PID",".png")),
       width = 1900,
       height = 1500,
       units = "px")

if(!exists(file.path(pca.dir, filt, "pca_5"))) dir.create(file.path(pca.dir, filt, "pca_5"), recursive = TRUE)


#PLOT TOP 5 PCAS
for (i in c("Plate", "Disease", "Timepoint", "condition", "age", "sex")){
  png(filename = file.path(pca.dir,filt,"pca_5", paste0(i, ".png")), width = 1200, height = 900, units = "px")
  clinical_pca <- cbind(pca_res$x, clin)
  variable <- as.factor(clinical_pca[,i])
  pairs(clinical_pca[,c(1:5)],
        pch = 19, 
        cex = 1, 
        oma=c(3,3,3,15),
        col = variable)
  par(xpd=TRUE)
  legend("bottomright", fill = unique(variable), legend = c( levels(variable)))
  dev.off()
}


##PCA Plots requested by Tess - Only Plate 1 and 2 ##

expression_plate1.2 <- exp[,c(row.names(clin)[clin$Plate == "1" | clin$Plate == "2"])]
clinical_plate1.2 <- clin[colnames(expression_plate1.2), ]

pca_res <- prcomp(t(
  as.data.frame(expression_plate1.2)), 
  scale. = TRUE, 
  center = TRUE) #center = TRUE

clinical_plate1.2_pca <- cbind(pca_res$x, clinical_plate1.2)

variable <- as.factor(clinical_plate1.2_pca[,"Plate"])

pcaplot <- pairs(clinical_plate1.2_pca[,c(1:5)],
                 pch = 19, 
                 cex = 1, 
                 oma=c(3,3,3,15),
                 col = variable)
par(xpd=TRUE)
legend("bottomright", fill = unique(variable), legend = c( levels(variable)))

# dev.off()
print(pcaplot)

# Correlation Matrix --------

library("PCAtools")


# EIGENCORPLOT -----------------------------------------------
library("PCAtools")

clinical_numeric <- clinical
clinical_numeric$condition <- as.factor(clinical_numeric$condition)
clinical_numeric$sex <- as.factor(clinical_numeric$sex)
clinical_numeric$Plate <- as.factor(clinical_numeric$Plate)
clinical_numeric$smoking_status <- as.factor(clinical_numeric$smoking_status)
clinical_numeric$age <- as.factor(clinical_numeric$age)

clinical_numeric$condition <- as.numeric(clinical_numeric$condition)
clinical_numeric$sex <- as.numeric(clinical_numeric$sex)
clinical_numeric$Plate <- as.numeric(clinical_numeric$Plate)
clinical_numeric$smoking_status <- as.numeric(clinical_numeric$smoking_status)
clinical_numeric$age <- as.numeric(clinical_numeric$age)

clinical_numeric <- clinical_numeric[,c("condition", "sex", "Plate", "smoking_status", "age")]


eigen.pca <- pca(expression,
                 metadata = clinical_numeric, 
                 center = TRUE,
                 scale = TRUE,
)


eigencorplot(eigen.pca,
             metavars = c(colnames(clinical_numeric)))





# ================================================================================== #
# 3. REMOVE 2 OUTLIERS AND EDIT FILES ==============================================
# ================================================================================== #

# #REMOVE OUTLIERS THAT WERE IDENTIFIED IN PCA PLOT FROM CLINICAL AND EXPRESSION FILES
clinical <- clinical[!row.names(clinical) == "3204316",]
clinical <- clinical[!row.names(clinical) == "3201323",]

expression <- expression[,row.names(clinical)]

write.csv(clinical, "data/processed/post-QC/clinical.csv")
write.csv(expression, "data/processed/post-QC/expression.csv")

#Paired samples
unpaired_samples <- c(colnames(HC_Plate2_ct), colnames(TB_Plate1_ct))

clinical_paired <- clinical[c(setdiff(row.names(clinical),unpaired_samples)),] #setdiff is the opposite of intersect
clinical_paired_wide <- pivot_wider(clinical_paired[,-c(5,21,22)], names_from = Timepoint, values_from = PAXGENE)
#-22 because of the 3200266 patient who has data in different plates
duplicated(clinical_paired_wide$PID)

#These 4 patients are meant to have paired T0/T6 samples but are missing expression data for T0 or T6. Removed
PID_missingpair <- clinical_paired_wide[(rowSums(is.na(clinical_paired_wide[,c(18:21)])) >=3),"PID"]
PID_missingpair <- PID_missingpair$PID

clinical_paired <- clinical_paired[-c(match(PID_missingpair, clinical_paired$PID)),]
expression_paired <- expression[,row.names(clinical_paired)]

write.csv(clinical_paired, "data/processed/post-QC/clinical_paired.csv")
write.csv(expression_paired, "data/processed/post-QC/expression_paired.csv")


# ================================================================================== #
# 4. CHANGE QC TO POST_FILT AND RERUN ALL QC ABOVE =================================
# ================================================================================== #




# ================================================================================== #
# 4. EXTRA UMAP PLOTS ==============================================================
# ================================================================================== #
#https://alexslemonade.github.io/refinebio-examples/03-rnaseq/dimension-reduction_rnaseq_02_umap.html
library("umap")

clinical_umap <- clinical
clinical_umap$PAXGENE <- as.character(clinical_umap$PAXGENE)
umap_results <- umap::umap(t(expression))

umap_plot_df <- data.frame(umap_results$layout) %>%
  # Turn sample IDs stored as row names into a column
  tibble::rownames_to_column("PAXGENE") %>%
  # Add the metadata into this data frame; match by sample IDs
  dplyr::inner_join(clinical_umap, by = "PAXGENE")

dir.create(file.path(figures.dir,"umap"), recursive = TRUE)
# for (i in c("Plate", "Disease", "Timepoint", "condition", "age")){

# Plate 
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = Plate)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/Plate.png"),
  width = 1900,
  height = 1500,
  units = "px")

# Disease 
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = Disease)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/Disease.png"),
  width = 1900,
  height = 1500,
  units = "px")

# Timepoiint
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = Timepoint)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/Timepoint.png"),
  width = 1900,
  height = 1500,
  units = "px")

#Sex
# Plate 
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = sex)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/Sex.png"),
  width = 1900,
  height = 1500,
  units = "px")

#smoke
# Plate 
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = smoking_status)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/smoking_status.png"),
  width = 1900,
  height = 1500,
  units = "px")

#Age
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = age)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/age.png"),
  width = 1900,
  height = 1500,
  units = "px")

#Condition
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = condition)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1), # Plot individual points to make a scatterplot
  
  file = paste0("figures/umap/condition.png"),
  width = 1900,
  height = 1500,
  units = "px")

#PID
ggsave(
  ggplot(umap_plot_df,
         aes(x = X1, y = X2, color = PID)) +
    labs(x = "UMAP1", y = "UMAP2") +
    geom_point(size = 1)  # Plot individual points to make a scatterplot
  + theme(legend.position = "none"),
  file = paste0("figures/umap/PID.png"),
  width = 1900,
  height = 1500,
  units = "px")

#Visualise timepoint + disease
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = Timepoint,
    shape = Disease
  )
) +
  geom_point() # Plot individual points to make a scatterplot

#make new file and run from here. dont overwrite old files


