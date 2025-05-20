library("ggplot2")
library("rstatix")
library("dplyr")
library("ggpubr")

genesig_D_7 <- c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP")

clinical <- read.csv("C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB/data/processed/post-QC/clinical.csv", row.names = 1)
expression <- as.matrix(read.csv("C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB/data/processed/post-QC/expression.csv", row.names = 1, check.names = FALSE))

#Normality tests from 4_TB_signature show not normal distribution

# Explanation of what the code below and scale() function is actually doing, using manual calculationS
# Example, manual calculation
# IFITM1 expression for sample 3200196 = 31.67308

# scale(data, center=T and scale = T) means we do the below:
## sample expression - mean(all IFITM1 expression values)/sd(all IFITM1 expression values)
# (31.67308 - mean(il17set1["IFITM1",]))/sd(il17set1["IFITM1",]) =  0.04839308
## 0.04839308 is the scaled and centered IFITM1 expression for sample 3200196
## Now we take the mean of tthe scaled + centered values of each sampleand plot that 

# scale(data,center=T and scale = F) means we do the below:
## sample expression - mean(all IFITM1 expression values)
# 31.67308 - mean(il17set1["IFITM1",]) = 0.1328462
## 0.04839308 is the centered IFITM1 expression for sample 3200196
## Now we take the mean of the centered values of each sample and plot that 



###Insert each of the signatures where "il17sig1" is and run this code for each of the signatures. 
# il17ind1<-match(il17sig1,genes)

expr_set<-expression[genesig_D_7,]

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
# 7. PLOT BOXPLOTS FOR GENESIG_D_7 =================================================
# ================================================================================== #

## SCALED & CENTERED =================================================
boxplot_data <- as.data.frame(cbind(score = mean_standardised["sig_scale",],
                                    group = as.character(clinical$condition),
                                    PID = as.character(clinical$PID)))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- list(
  c("HC_T0", "HC_T6"),
  c("HC_T0", "TB_T0"),
  c("TB_T0", "TB_T2"),
  c("TB_T0", "TB_T4"),
  c("TB_T0", "TB_T6")
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



genesig_D_7_figures.dir <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB/results/roc_curve/TBT0vsHCT0/figures/genesig_D_7"


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

ggsave(boxplotfinal2, filename = file.path(genesig_D_7_figures.dir, "scaledcentered_signature_all_paired.png"), 
       width = 3500, 
       height = 3200, 
       units = "px" )







## CENTERED =================================================
boxplot_data <- as.data.frame(cbind(score = mean_standardised["sig_center",],
                                    group = as.character(clinical$condition),
                                    PID = as.character(clinical$PID)))



gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_comparisons <- list(
  c("HC_T0", "HC_T6"),
  c("HC_T0", "TB_T0"),
  c("TB_T0", "TB_T2"),
  c("TB_T0", "TB_T4"),
  c("TB_T0", "TB_T6")
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


genesig_D_7_figures.dir <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB/results/roc_curve/TBT0vsHCT0/figures/genesig_D_7"


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

ggsave(boxplotfinal2, filename = file.path(genesig_D_7_figures.dir, "centered_signature_all_paired.png"), 
       width = 3500, 
       height = 3200, 
       units = "px" )
























heatmapsamples <- clinical$PAXGENE[which(clinical$condition == "TB_T0" | clinical$condition == "HC_T0")]
heatmapsamples <- clinical$PAXGENE[which(clinical$Disease == "TB")]

# selection <-which((tT$logFC>1|tT$logFC< -1)& tT$adj.P.Val<0.05)
# 
# tT2=tT[selection,]

clinical_heatmap <- clinical[as.character(heatmapsamples),]
expression_heatmap <- expression[,as.character(heatmapsamples)]

clinical_heatmap_ordered <- clinical_heatmap[order(clinical_heatmap$Timepoint,clinical_heatmap$Disease, clinical_heatmap$sex),]

clinical_heatmap_ordered <- clinical_heatmap_ordered[-which(is.na(clinical_heatmap_ordered$sex)),]

#Make legend
#Make a vector of timepoints in the same order as the clincal_ordered table. Change each timepoint to a different colour
Timepoint = as.character(clinical_heatmap_ordered$Timepoint)
Timepoint[Timepoint == "T0"] <- "lightgrey"
Timepoint[Timepoint == "T2"] <- "lightblue3"
Timepoint[Timepoint == "T4"] <- "skyblue3"
Timepoint[Timepoint == "T6"] <- "steelblue4"

#Make a vector of disease in the same order as the clincal_ordered table. Change each disease to a different colour
Disease = as.character(clinical_heatmap_ordered$Disease)
Disease[Disease == "HC"] <- "seagreen"
Disease[Disease == "TB"] <- "red3"

Sex = as.character(clinical_heatmap_ordered$sex)
Sex[Sex == "male"] <- "purple"
Sex[Sex == "female"] <- "yellow2"

labs = cbind(Timepoint,Disease, Sex)

#reorder expression columns/samples to be the same as clinical_ordered
heatmapdata <- as.matrix(expression_heatmap[,row.names(clinical_heatmap_ordered)])
heatmapdata <- heatmapdata[-which(row.names(heatmapdata) == "B2M" | row.names(heatmapdata) == "GAPDH"),]


library(heatmap3)

png(filename = file.path(heatmap.dir, paste0("heatmap_withsex","TB_T0","and","HC_T0",".png")),
    width = 1200,
    height = 1000,
)

heatmap3(heatmapdata, 
         Colv=NA, 
         scale = "row",
         # lasRow = 
         balanceColor=T,
         labCol=NA,
         showColDendro = F, 
         showRowDendro = F,
         margins=c(1,1),
         ColSideLabs = F, 
         ColSideColors =labs,
         cexRow=1.5,
         legendfun=function()
           showLegend(legend=c("HC","TB","T0", "T2", "T4", "T6", "Male", "Female"),
                      col=c("seagreen","red3","lightgrey", "lightblue","skyblue3", "steelblue4", "purple", "yellow2"),
                      cex=1.5)
)
dev.off()

