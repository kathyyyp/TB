

# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

# #Mac
my_directory <- "/Users/kathyphung/Library/CloudStorage/OneDrive-UTS/Documents/RBMB/TB"
setwd(file.path(my_directory))
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
library(limma)

# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
main.dir <- my_directory
data.dir <- file.path(my_directory, "data")
processed.dir <- file.path(my_directory, "data", "processed")
output.dir <- file.path(my_directory,"output")
validation.dir <- file.path(output.dir,"validation", "longitudinal")
if(!exists(validation.dir)) dir.create(validation.dir)



# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #
setwd(file.path(data.dir,"raw"))

#clinical
HC_T0_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "HCT0", col_names = TRUE) 
TB_T0_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "TBT0", col_names = TRUE) 
HC_T0_T6_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "HCT0T6", col_names = TRUE) 
T4P_Plate1_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "T4P Plate 1", col_names = TRUE) 
T4P_Plate2_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "T4P Plate 2", col_names = TRUE) 
T4P_Plate3_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "T4P Plate 3", col_names = TRUE) 
T4P_Plate4_meta <- read_excel(file.path(data.dir, "raw", "RNA signature validation study pt infos 16122025 Kathy.xlsx"), sheet = "T4P Plate 4", col_names = TRUE) 



HC_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "HCT0", skip = 1, col_names = TRUE) 
TB_T0_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "TBT0", skip = 1, col_names = TRUE) 
HC_T0_T6_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "HBT0T6", col_names = TRUE) 
T4P_Plate1_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "T4P Plate 1", col_names = FALSE) 
T4P_Plate2_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "T4P Plate 2", col_names = FALSE) 
T4P_Plate3_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "T4P Plate 3", col_names = FALSE) 
T4P_Plate4_ct <- read_excel(file.path(data.dir, "raw", "Biomarkers Raw Data for Analyse v 1.4 - 161225.xlsx"), sheet = "T4P Plate 4", col_names = FALSE) 
# Note that there are 3 wells that were water



setwd(file.path(main.dir))


# ================================================================================== #
# 2. CLEAN UP EXPRESSION DATA ======================================================
# ================================================================================== #
clinical <- as.data.frame(rbind(HC_T0_meta,
                  HC_T0_T6_meta,
                  TB_T0_meta,
                  T4P_Plate1_meta,
                  T4P_Plate2_meta,
                  T4P_Plate3_meta,
                  T4P_Plate4_meta))


#Check duplicates
clinical[duplicated(clinical$sample),"sample"]

#3201009 is in unpaired TBT0 twice - remove one of them
which(clinical$sample == "3201009") #149 and 155
clinical <- clinical[-which(clinical$sample == "3201009")[2],]

#Check duplicates now again
clinical[duplicated(clinical$sample),"sample"]

clinical$disease_simple <- clinical$disease
clinical$disease <- paste0(clinical$disease_simple, "_", clinical$timepoint)

#Check that clinical timepoints match up with data timepoints
cbind(clinical[match(as.character(T4P_Plate1_ct[2,-1]), clinical$sample), "timepoint"], as.character(T4P_Plate1_ct[1,-1]))
cbind(clinical[match(as.character(T4P_Plate2_ct[2,-1]), clinical$sample), "timepoint"], as.character(T4P_Plate2_ct[1,-1]))
cbind(clinical[match(as.character(T4P_Plate3_ct[2,-1]), clinical$sample), "timepoint"], as.character(T4P_Plate3_ct[1,-1]))
cbind(clinical[match(as.character(T4P_Plate4_ct[2,-1]), clinical$sample), "timepoint"], as.character(T4P_Plate4_ct[1,-1]))


wrangling_func <- function(x){
  #quick check that there are 4 samples in a row (no patient missing a timepoint)
  print(all(x[1,-1] == rep(c("T0", "T2", "T4","T6"), ncol(x[-1])/4))) 
  
  #clean up
  colnames(x)[-1]<- x[2,-1]
  x <- x[-c(1,2),]
  row.names(x) <- NULL
  return(x)
}


T4P_Plate1 <- wrangling_func(T4P_Plate1_ct)
T4P_Plate2 <- wrangling_func(T4P_Plate2_ct)
T4P_Plate3 <- wrangling_func(T4P_Plate3_ct)
T4P_Plate4 <- wrangling_func(T4P_Plate4_ct)
 
#check for dups
colnames(TB_T0_ct)[duplicated(colnames(TB_T0_ct))]
colnames(HC_T0_ct)[duplicated(colnames(HC_T0_ct))]
colnames(HC_T0_T6_ct)[duplicated(colnames(HC_T0_T6_ct))]

#T4P_Plate1 has 4 water samples
colnames(T4P_Plate1)[duplicated(colnames(T4P_Plate1))]
T4P_Plate1 <- T4P_Plate1[,-which(colnames(T4P_Plate1) == "water")]

colnames(T4P_Plate2)[duplicated(colnames(T4P_Plate2))]
colnames(T4P_Plate3)[duplicated(colnames(T4P_Plate3))]
colnames(T4P_Plate4)[duplicated(colnames(T4P_Plate4))]


listoffiles <- list(HC_T0 = HC_T0_ct, 
                    TB_T0 = TB_T0_ct,
                    HC_T0_T6 = HC_T0_T6_ct %>% mutate(across(-1, ~ as.numeric(.))),
                    T4P_Plate1 = T4P_Plate1 %>% mutate(across(-1, ~ as.numeric(.))),
                    T4P_Plate2 = T4P_Plate2 %>% mutate(across(-1, ~ as.numeric(.))),
                    T4P_Plate3 = T4P_Plate3 %>% mutate(across(-1, ~ as.numeric(.))),
                    T4P_Plate4 = T4P_Plate4 %>% mutate(across(-1, ~ as.numeric(.))))


#Average the duplicated genes
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

###  NORMALISE 
### 17/12/25 did not include the unpaired HC_T0 and TB_T0 samples yet. Waiting for Tess to resolve housekeeping gene disparity
all_data <- cbind(
                  # HC_T0,
                  # TB_T0,
                  HC_T0_T6,
                  T4P_Plate1,
                  T4P_Plate2 ,
                  T4P_Plate3,
                  T4P_Plate4)

colnames(all_data)[duplicated(colnames(all_data))]

#Remove NAs
which(is.na(as.matrix(all_data)), arr.ind = TRUE)

# "3200648" Has missing data for at least one gene. Remove 
colnames(all_data)[c(187)]
all_data[,c(187)]
all_data <- all_data[,-c(187)]  #NA values for expression

# When HC_T0 and TB_T0 are included
# # "3201400" "3205033" "3200648" Has missing data for at least one gene. Remove 
# colnames(all_data)[c(13, 72, 292)]
# all_data[,c(13, 72, 292)]
# all_data <- all_data[,-c(13, 72, 292)]  #NA values for expression



which(is.na(as.matrix(all_data)), arr.ind = TRUE)


#Match up with clinical file
row.names(clinical) <- clinical$sample
dim(clinical) #445
clinical <- clinical[intersect(colnames(all_data), clinical$sample),]
dim(clinical) # 341 if only including paired (17/12/25). 442 if including unpaired (later)

all_data <- all_data[,clinical$sample]

write.csv(all_data, file.path(validation.dir, "all_data_beforenorm.csv"))

# ================================================================================== #
# 5. NORMALISE (delta Ct) ==========================================================
# ================================================================================== #

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
write.csv(normdata_twohk, file.path(validation.dir, "normdata_twohk.csv"))


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

## 2) Get groups to be compared --------------------------------------------

listofresults <- list()
for (hk in names(listof_normdata)){
  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(figures.dir)) dir.create(figures.dir)
  
  normdata <- listof_normdata[[hk]]
  expression <- as.matrix(normdata)
  
  clinical$group <-clinical$disease
  table(clinical$group)
  
  # View(pivot_wider(clinical[,-2], names_from = "group", values_from = "sample"))
  colnames(expression) == row.names(clinical)
  
  write.csv(expression, file.path(this.output.dir, paste0(hk,"_relative_expression.csv")))
  write.csv(clinical, file.path(this.output.dir, "clinical.csv"))
  
  
  expr_long <- as.data.frame(expression) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to="sample_id", values_to="expression") 
  
  
  expr_long <- cbind(expr_long, clinical[match(expr_long$sample_id, clinical$sample),])
  
  
    #Remove outlier !!!!!!!!
    #When plotting boxplot, this was an extreme outlier for B2M, so went back and removed it here
    expr_long <- expr_long[-which(expr_long$sample_id == "3200200"),]
    
  
  
  ## Per-comparison filtering for Wilcoxin paired statistical test ---------------------------------------------------------------------------------------
  #First, need to create a function for per-comparison filtering
  #Create function to get ONLY paired data between two groups 
  get_paired_data <- function(data, group1, group2) {
    ids1 <- data$PID[data$group == group1]
    ids2 <- data$PID[data$group == group2]
    paired_ids <- intersect(ids1, ids2)
    
    data %>%
      filter(group %in% c(group1, group2), PID %in% paired_ids)
  }
  
  ## Make stat.table.gsva ---------------------------------------------------------------------------------------
  #map_dfr will apply the below function to each element of list my_paired_comparisons, and directly rbind the results
  #In this function, 'groups' is each element of the list. so groups[1] should be "HC_T0"
  #Using map_dfr here is the same as if we wrote stat_table_function <- function(groups){...} and then did stat_table_function(groups = my_paired_comparisons[1]) and stat_table_function(groups = my_paired_comparisons[2]) etc etc for all 4 comparisons, and then rbinded the output together
 
  my_paired_comparisons = list(
    c("HC_T0", "HC_T6"),
    c("TB_T6", "TB_T0"),
    c("TB_T4", "TB_T0"),
    c("TB_T2", "TB_T0"),
    c("TB_T6", "TB_T4"),
    c("TB_T6", "TB_T2"),
    c("TB_T4", "TB_T2"))
  
  
my_unpaired_comparisons <- list(
  c("HC_T0", "TB_T0"), 
  c("HC_T6", "TB_T6")
)

##  MAP TO EACH GENE   ===================================================================

  stat.table.full <- map_dfr(unique(expr_long$gene), function(this_gene) { #close map for each gene
    
    boxplot_data <- expr_long[expr_long$gene == this_gene,]
    
##  MAP TO EACH PAIRES COMPARISON   ===================================================================

   stat.table <- map_dfr(my_paired_comparisons, function(groups) { #close map for each comparison
     
# ================================================================================== #
##  PAIRED COMPARISONS  ===================================================================
# ================================================================================== #

boxplot_data$expression <- as.numeric(boxplot_data$expression)
boxplot_data$group <- factor(boxplot_data$group)

    #1) extract the timepoint name
    g1 <- groups[1]
    g2 <- groups[2]
    
    #2) make a data frame only containing our paired samples, making use of the function we made above
     # wilcox test with paired= TRUE needs two groups. it assumes row order corresponds to pair samples. so patient hc1 T0 must be followed by hc1 T6
    df <- get_paired_data(data = boxplot_data, 
                          group1 = g1, 
                          group2 = g2)
    
    df$group <- as.character(df$group) #don't know why this is needed. HC_T0 vs HC_T6 works without this. but the rest don't. something to do with unused factor levels?
    
    #3) run wilcox_test on the paired samples
    if (n_distinct(df$PID) > 0) {
      wilcox_test(df, expression ~ group, 
                  paired = TRUE) %>%
        add_xy_position(x = "group", step.increase = 0.1) %>% 
        mutate(group1 = g1, group2 = g2)
    } else {
      tibble()  # return empty if no paired data
    }
    
    
    
  }) #close map for each comparison
  
   
   
   
   
   
# ================================================================================== #
##  UNPAIRED COMPARISONS  ===================================================================
# ================================================================================== #
stat.table2 <- wilcox_test(boxplot_data, expression ~ group, 
                                paired = FALSE,
                                comparisons = my_unpaired_comparisons) %>%
        add_xy_position(x = "group")
stat.table2 <- stat.table2[,!colnames(stat.table2) %in%
                                          c("p.adj","p.adj.signif")]

stat.table.all <- rbind(stat.table, stat.table2)

# stat.table <- stat.table[which(stat.table$p < 0.05),]
# lowest_bracket <- max(boxplot_data$expression) + 0.05*(max(boxplot_data$expression))
# stat.table.all$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.all))


   
   # Compute max y for each comparison
   stat.table.all <- stat.table.all %>%
     rowwise() %>% # groupdataframe by row
     mutate(y.max = max(boxplot_data$expression[boxplot_data$group %in% c(group1, group2)], na.rm = TRUE)) %>%
     ungroup()
   
   
   # Assign y positions
   # 1) HC bracket
   hc_id <- which(stat.table.all$group1 == "HC_T0" & stat.table.all$group2 == "HC_T6")
   stat.table.all$y.position <- NA
   if(length(hc_id) == 1){
     step_hc <- 0.08 * stat.table.all$y.max[hc_id]
     stat.table.all$y.position[hc_id] <- stat.table.all$y.max[hc_id] + step_hc
   }
   # 2) TB brackets
   tb_id <- setdiff(seq_len(nrow(stat.table.all)), hc_id)
   if(length(tb_id) > 0){
     # first TB bracket above its own max
     step_tb_first <- 0.08 * stat.table.all$y.max[tb_id[1]]
     stat.table.all$y.position[tb_id[1]] <- stat.table.all$y.max[tb_id[1]] + step_tb_first
     
     # remaining TB brackets stagger incrementally using their own max * 0.1
     if(length(tb_id) > 1){
       for(i in 2:length(tb_id)){
         step_tb <- 0.08 * stat.table.all$y.max[tb_id[i]]
         stat.table.all$y.position[tb_id[i]] <- stat.table.all$y.position[tb_id[i-1]] + step_tb
       }
     }
   } #close TB bracket y pos
   
   
   
   
   # Add gene name to the results
   stat.table.all$gene <- this_gene
   
   stat.table.all
  }) #close map for each gene
  
  #Fix up x axis - all the comparisons were overlapping
  group_levels <- sort(unique(expr_long$group))
  
  stat.table.full <- stat.table.full %>%
    mutate(
      xmin = match(group1, group_levels),
      xmax = match(group2, group_levels)
    )
  
  
  
  
  
  boxplot_theme <- theme(axis.title = element_text(size = 20),
                      axis.text = element_text(size = 15),
                      title = element_text(size = 20),
                      strip.text.x = element_text(size = 20),
                      legend.text = element_text(size = 15),
                      legend.position = "bottom",
                      axis.text.x = element_blank(),      
                      axis.ticks.x = element_blank()) 
  
  plot <- ggplot(expr_long, aes(x = disease, y = expression, color=disease
                                # , label = sample
                                ))+
    geom_boxplot(outlier.shape=NA, alpha=0.2) +
    geom_jitter(width=0.2, size=1) +
    theme_bw() +
    boxplot_theme +
    # geom_text() + #to see which samples are the outliers (add label = sample to aes)
  
    stat_summary(fun.y = mean, fill = "red",
                 geom = "point", shape = 21, size =2,
                 show.legend = TRUE) +
      stat_pvalue_manual(
        stat.table.full,
        label = "p",
        xmin = "xmin",
        xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        size = 4.5
      ) +
    facet_wrap(~gene, scales="free_y", ncol = 7) +
      
    guides(color = guide_legend(nrow = 1))+
    ylab("Expression (2^-delta Ct)") +
    xlab("Disease") +
    labs(caption = paste("Relative expression normalized to", hk)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))
  
  
  ggsave(plot, filename= file.path(figures.dir,paste0(hk,"_all_genes_boxplot.png")),
         width = 60, height = 20, units = "cm")

  
  # ================================================================================== #
  # 7) SCALED+CENTERED AND CENTERED MEAN  =================================================
  # ================================================================================== #
  

  scaledcentered_mean_func <- function(number_of_genes = "7"){
    
    expr_set<-expression[row.names(expression),]
    
    if(number_of_genes == "6"){
      expr_set <- expr_set[!row.names(expr_set) == "S100A8",]
    }
    
    if(number_of_genes == "5"){
      expr_set <- expr_set[-c(which(row.names(expr_set) == "S100A8" | row.names(expr_set) == "CD274")),]
    }
    
    if(number_of_genes == "4"){
      expr_set <- expr_set[-c(which(row.names(expr_set) == "S100A8" | row.names(expr_set) == "CD274" | row.names(expr_set) == "IFITM1")),]
    }
    #transpose for scaling
    expr_set<-t(expr_set)
    

    
    #Two methods
    # 1) center and scale 
    # 2) center only     #centers and/or scales the data based on its mean and standard deviation (so that it has a mean of 0 and a standard deviation of 1)
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
  listofstandardised_scores[["6_genes"]] <- scaledcentered_mean_func(number_of_genes = "6")
  listofstandardised_scores[["5_genes"]] <- scaledcentered_mean_func(number_of_genes = "5")
  listofstandardised_scores[["4_genes"]] <- scaledcentered_mean_func(number_of_genes = "4")
  
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
  for(g in names(listofstandardised_scores)){
    
    
    score_data <- listofstandardised_scores[[g]][["scores"]]
    
    gene_list <- listofstandardised_scores[[g]][["gene_list"]]
    
    boxplot_data <- as.data.frame(cbind(score = score_data["sig_scale",],
                                        group = as.character(clinical$disease),
                                        sample = as.character(clinical$sample),
                                        PID = as.character(clinical$PID)))
    
        boxplot_data <- boxplot_data[-which(boxplot_data$sample == "3200200"),]

    
    
    gsva_theme <- theme(axis.title = element_text(size = 20),
                        axis.text = element_text(size = 20),
                        title = element_text(size = 20),
                        legend.position = "None") 
    
    
    my_paired_comparisons = list(
      c("HC_T0", "HC_T6"),
      c("TB_T6", "TB_T0"),
      c("TB_T4", "TB_T0"),
      c("TB_T2", "TB_T0"),
      c("TB_T6", "TB_T4"),
      c("TB_T6", "TB_T2"),
      c("TB_T4", "TB_T2"))
    
    
  
my_unpaired_comparisons <- list(
  c("HC_T0", "TB_T0"), 
  c("HC_T6", "TB_T6")
)

    boxplot_data$score <- as.numeric(boxplot_data$score)
    boxplot_data$group <- factor(boxplot_data$group)
    

      stat.table.gsva <- map_dfr(my_paired_comparisons, function(groups) {
        
        #1) extract the timepoint name
        g1 <- groups[1]
        g2 <- groups[2]
        
        #2) make a data frame only containing our paired samples, making use of the function we made above
        # wilcox test with paired= TRUE needs two groups. it assumes row order corresponds to pair samples. so patient hc1 T0 must be followed by hc1 T6
        df <- get_paired_data(data = boxplot_data, 
                              group1 = g1, 
                              group2 = g2)
        
        df$group <- as.character(df$group) #don't know why this is needed. HC_T0 vs HC_T6 works without this. but the rest don't. something to do with unused factor levels?
        
        #3) run wilcox_test on the paired samples
        if (n_distinct(df$PID) > 0) {
          wilcox_test(df, score ~ group, 
                      paired = TRUE) %>%
            add_xy_position(x = "group", step.increase = 0.1) %>% 
            mutate(group1 = g1, group2 = g2)
        } else {
          tibble()  # return empty if no paired data
        }
      }) #close map for paired comparison
      
# ================================================================================== #
##  UNPAIRED COMPARISONS  ===================================================================
# ================================================================================== #
stat.table.gsva2 <- wilcox_test(boxplot_data, score ~ group, 
                                paired = FALSE,
                                comparisons = my_unpaired_comparisons) %>%
        add_xy_position(x = "group")
stat.table.gsva2 <- stat.table.gsva2[,!colnames(stat.table.gsva2) %in%
                                          c("p.adj","p.adj.signif")]

stat.table.gsva.all <- rbind(stat.table.gsva, stat.table.gsva2)


   
   # Compute max y for each comparison
   stat.table.gsva.all <- stat.table.gsva.all %>%
     rowwise() %>% # groupdataframe by row
     mutate(y.max = max(boxplot_data$expression[boxplot_data$group %in% c(group1, group2)], na.rm = TRUE)) %>%
     ungroup()
      
      
      
    
    #Fix up x axis - all the comparisons were overlapping
    group_levels <- sort(unique(boxplot_data$group))
    
    stat.table.gsva.all <- stat.table.gsva.all %>%
      mutate(
        xmin = match(group1, group_levels),
        xmax = match(group2, group_levels)
      )
    
    
    #FOR SCALED+CENTERED
    
    # Compute max y for each comparison
    stat.table.gsva.all <- stat.table.gsva.all %>%
      rowwise() %>% # groupdataframe by row
      mutate(y.max = max(boxplot_data$score[boxplot_data$group %in% c(group1, group2)], na.rm = TRUE)) %>%
      ungroup()
    #HC bracket
    # Assign y positions
    # 1) HC bracket
    hc_id <- which(stat.table.gsva.all$group1 == "HC_T0" & stat.table.gsva.all$group2 == "HC_T6")
    stat.table.gsva.all$y.position <- NA
    if(length(hc_id) == 1){
      step_hc <- 0.1 * stat.table.gsva.all$y.max[hc_id]
      stat.table.gsva.all$y.position[hc_id] <- stat.table.gsva.all$y.max[hc_id] + step_hc
    }
    # 2) TB brackets
    tb_id <- setdiff(seq_len(nrow(stat.table.gsva.all)), hc_id)
    if(length(tb_id) > 0){
      # first TB bracket above its own max
      step_tb_first <- 0.1 * stat.table.gsva.all$y.max[tb_id[1]]
      stat.table.gsva.all$y.position[tb_id[1]] <- stat.table.gsva.all$y.max[tb_id[1]] + step_tb_first
      
      # remaining TB brackets stagger incrementally using their own max * 0.1
      if(length(tb_id) > 1){
        for(i in 2:length(tb_id)){
          step_tb <- 0.1 * stat.table.gsva.all$y.max[tb_id[i]]
          stat.table.gsva.all$y.position[tb_id[i]] <- stat.table.gsva.all$y.position[tb_id[i-1]] + step_tb
        }
      }
    }
    
    lowest_bracket <- max(boxplot_data$score) + 0.05*(max(boxplot_data$score))
    stat.table$y.position <- seq(lowest_bracket, by= 0.3, length.out = nrow(stat.table))
    
    
    boxplot_scaledcentered <- ggplot(boxplot_data, aes(
      x = factor(group),
      # x = factor(group),
      y = as.numeric(boxplot_data[,"score"]))) +
      
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
      
      stat_pvalue_manual(stat.table.gsva.all,
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
    
    listofboxplots_scaledcentered[[g]] <- ggplotGrob(boxplot_scaledcentered) #ggGrob freezes the image in place, otherwise the pvalue brackets move when put into the list
    
    
  } #close loop for number of genes
  
  boxplot_scaledcentered_panel <- annotate_figure(
    ggarrange(
      plotlist = list(listofboxplots_scaledcentered[["7_genes"]],
                      listofboxplots_scaledcentered[["6_genes"]],
                      listofboxplots_scaledcentered[["5_genes"]],
                      listofboxplots_scaledcentered[["4_genes"]]),
      ncol = 4),
    bottom = text_grob(paste("Relative expression normalized to", hk),
                       hjust = 1, x = 1, size = 16))
  
  ggsave(boxplot_scaledcentered_panel, filename = file.path(this.figures.dir, paste0(hk,"_boxplot_mean_ztransformed_panel.png")),
         width = 8000,
         height = 2000,
         units = "px" )
  
  
} #close loop for hk gene









## 4) Validation  ------------------------------------------------------

# --- PAIRWISE ROC ANALYSIS --- #
# Define all pairwise comparisons of interest

pairwise_comparisons = list(
  c("HC_T0", "TB_T0"), 
  c("HC_T6", "TB_T6"),
  c("HC_T0", "HC_T6"),
  c("TB_T6", "TB_T0"),
  c("TB_T4", "TB_T0"),
  c("TB_T2", "TB_T0"),
  c("TB_T6", "TB_T4"),
  c("TB_T6", "TB_T2"),
  c("TB_T4", "TB_T2"))

  

for (hk in names(listofresults)){ # housekeeping gene loop
  
  listofstandardised_scores <- listofresults[[hk]]
  listofboxplots_scaledcentered <- list()
  listofboxplots_centered <- list()
  
  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  this.figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(this.figures.dir)) dir.create(this.figures.dir)
  
  # Create a list to store AUC values and roc objects
  
  res_table <- data.frame() 
  roc_objects <- list()

  for(g in names(listofstandardised_scores)){ # number of genes loop
      forestplot_res_table <- data.frame()

    score_data <- listofstandardised_scores[[g]][["scores"]]
    
    # Loop through each pairwise comparison
    for (pair in pairwise_comparisons) { #pairwise comparison loop
      
      
      group1 <- pair[1]
      group2 <- pair[2]
      
      # Subset data to omly include the 2 rgroups of interest
      subset_clinical <- clinical[clinical$group %in% c(group1,group2),]
      subset_counts <- score_data[, row.names(subset_clinical)]
      
      subset_clinical$group <- factor(subset_clinical$group, levels = c(group1, group2))
      
      
      
    
    gene_list <- paste(
      "7_genes:", paste0(listofstandardised_scores[["7_genes"]][["gene_list"]], collapse = ","), "\n",
      "6_genes:", paste0(listofstandardised_scores[["6_genes"]][["gene_list"]], collapse = ","), "\n",
      "5_genes:", paste0(listofstandardised_scores[["5_genes"]][["gene_list"]], collapse = ","), "\n",
      "4_genes:", paste0(listofstandardised_scores[["4_genes"]][["gene_list"]], collapse = ","), "\n"
      
    )
    
    mean_standardised_scores <- as.data.frame(t(as.matrix(subset_counts)))
    mean_standardised_scores$group <- clinical[colnames(subset_counts), "group"]
    
    mean_standardised_scores$group <- factor(mean_standardised_scores$group, levels = c(group1, group2))
    
    #scaled+centered
    normtype = "sig_scale"
    glm_model <- glm(group ~ sig_scale, data = mean_standardised_scores, family = binomial)
    

  test_probs <- predict(glm_model, type = "response")
  
  roc_obj <- roc(mean_standardised_scores$group, test_probs)
  
  plot(roc_obj)
  auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  #  The "optimal threshold" refers to the point on the ROC curve where you achieve the best balance between sensitivity and specificity, or where the classifier is most effective at distinguishing between the positive and negative classes.
  optimal_threshold_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity", best.method = "youden"))
  
  if(nrow(optimal_threshold_coords) > 1) {
    optimal_threshold_coords <- optimal_threshold_coords[1,] # some output have 2 equally optimal thresholds = same AUC. just keep  first one as results are the same
  }
  
  
# Sensitivity confidence interval
ci_sens <- ci.se(roc_obj, specificities =  as.numeric(optimal_threshold_coords["specificity"]))

# Specificity confidence interval
ci_spec <- ci.sp(roc_obj, sensitivities =  as.numeric(optimal_threshold_coords["sensitivity"]))

  
  res_current <-cbind(
    comparison = paste0(group1," vs ",group2),
    auc = auc(roc_obj),
    ci = paste0(round(as.numeric(auc_ci[1]),2), "-", round(as.numeric(auc_ci[3]),2)),
    sensitivity = optimal_threshold_coords$sensitivity, 
    specificity = optimal_threshold_coords$specificity
    
  )
  
  res_table <- rbind(res_table, res_current)
  write.csv(res_table, file.path(this.output.dir, paste0(hk,"_", g, "_mean_ztransformed_scores_res_table.csv")))
  
  

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

write.csv(forestplot_res_table, file.path(this.output.dir, paste0(hk,"_", g, "_mean_ztransformed_scores_forestplot_res_table.csv")))

  roc_objects[[paste0(group1," vs ",group2)]] <- roc_obj
  
}



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
disease_roc_data <- roc_data[which(roc_data$Comparison == "HC_T0 vs TB_T0" | 
                                     roc_data$Comparison == "HC_T6 vs TB_T6"), ]


disease_roc <- ggplot(disease_roc_data, aes(x = FPR, y = TPR, color = legend)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black")  +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 16),
        title = element_text(size = 20)) +
  labs(
    title = "ROC - TB vs HC",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = paste(g, "Signature:", paste0(listofstandardised_scores[[g]][["gene_list"]], collapse = ",")))



ggsave(disease_roc, filename = file.path(this.figures.dir, paste0(hk, "_", g, "_", normtype, "_disease_ROC_mean_z_transformed_scores.png")),
       width = 2500,
       height = 3000,
       units = "px" )



# Timepoints ROC 
timepoints_roc_data <- roc_data[which(roc_data$Comparison != "HC_T0 vs TB_T0" & 
                                      roc_data$Comparison != "HC_T6 vs TB_T6"), ]


timepoints_roc <- ggplot(timepoints_roc_data, aes(x = FPR, y = TPR, color = legend)) +
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
    title = "ROC - TB vs HC",
    x = "FPR (1 - Specificity)",
    y = "TPR (Sensitivity)",
    color = "Comparison",
    caption = paste(g, "Signature:", paste0(listofstandardised_scores[[g]][["gene_list"]], collapse = ",")))

ggsave(timepoints_roc, filename = file.path(this.figures.dir, paste0(hk, "_", g, "_timepoints_ROC_mean_z_transformed_scores.png")),
       width = 2500,
       height = 3000,
       units = "px" )





  } #close number of genes loop
} #close housekeeping gene loop





### Forest plots ####

## Prep table


for (hk in names(listofresults)){ # housekeeping gene loop
  
  this.output.dir <- file.path(validation.dir, hk)
  if(!exists(this.output.dir)) dir.create(this.output.dir)
  
  this.figures.dir <- file.path(this.output.dir, "figures")
  if(!exists(this.figures.dir)) dir.create(this.figures.dir)
  
    listofstandardised_scores <- listofresults[[hk]]

  for(g in names(listofstandardised_scores)){ # number of genes loop
    

res_table <- read.csv(file.path(this.output.dir, paste0(hk,"_",g,"_mean_ztransformed_scores_forestplot_res_table.csv")), row.names = 1)

auc_plot <- res_table %>% 
  ggplot(aes(y = comparison)) + 
  theme_bw() +
  geom_point(aes(x=auc), shape=15, size=3) +
  geom_linerange(aes(xmin=auc_ci_low, xmax=auc_ci_high)) +
  ylab(NULL)+
  xlab("AUC")+
coord_cartesian( xlim=c(0.1, 1))

sens_plot <- res_table %>% 
  ggplot(aes(y = comparison)) + 
  theme_bw() +
  geom_point(aes(x=sensitivity), shape=15, size=3) +
  geom_linerange(aes(xmin=sensitivity_ci_low, xmax=sensitivity_ci_high)) +
  ylab(NULL)+
  xlab("Specificity")+
coord_cartesian( xlim=c(0.1, 1))

spec_plot <- res_table %>% 
  ggplot(aes(y = comparison)) + 
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
                               paste(g, "Signature:", paste0(listofstandardised_scores[[g]][["gene_list"]], collapse = ","))),
                      size = 8, hjust = 0, x = 0)
)


ggsave(panel_forest, filename= file.path(this.figures.dir, paste0(hk, "_", g, "_","_forestplot_panel.png")),
       width = 10, height = 15, units = "cm",   bg = "white"  )

  }
}


#to add to git 
#reran analysis with corrected samples, removed an outlier, commented out centered-only data, ran analysis on paired data
