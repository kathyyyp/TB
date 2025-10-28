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

#CONTINUE ON FROM SCRIPT 1_TB_DATAWRANGLING.r



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

#Block by plate for paired plates. Can't do for TB vs HC because TB is plate 1 and HC is plate 2 (batch effect and disease effect is the same then)
corfit_paired <- duplicateCorrelation(
  expression_paired2,
  design_paired,
  block = clinical_paired2$Plate) #block by plate? #plate 1 was TB and plate 2 was HC


fit <- lmFit(expression2, 
             design = design
             # block = clinical2$Disease,
             # correlation = corfit$consensus
)  

fit_paired <- lmFit(expression_paired2, 
                    design = design_paired,
                    block = clinical_paired2$Plate,
                    correlation = corfit_paired$consensus
)  

cont.matrix <- makeContrasts(
  contrast1 = TB_T0 - HC_T0, 
  levels = design) #all data (paired + unpaired) - taking ALL TB_T0 and HC_T0 and comparing them

cont.matrix_paired <- makeContrasts(
  contrast2 = HC_T6 - HC_T0,
  contrast3 = TB_T6 - TB_T0,
  contrast4 = TB_T6 - HC_T6,
  contrast5 = TB_T6 - TB_T2,
  contrast6 = TB_T4 - TB_T0,
  contrast7 = TB_T2 - TB_T0,
  contrast8 = TB_T6 - TB_T4,
  contrast9 = TB_T4 - TB_T2,
  #which genes respond differently over time in Tb compared to household
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
    # selection <-which((tT$logFC>1|tT$logFC< -1)&tT$FDR<0.05)
    
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

listoftT2[["contrast1"]][which(listoftT2[["contrast1"]]$legend == "Upregulated"), "genename"]

# save.image("workspaces/07_04_2025_DEresults.RData")




# ================================================================================== #
## 2.2. BOXLOTS ====================================================================
# ================================================================================== #


boxplotdata_unpaired <- as.data.frame(t(expression))

boxplotdata_unpaired <- boxplotdata_unpaired[row.names(clinical),]

boxplotdata_unpaired <- cbind(boxplotdata_unpaired,
                              PID = clinical$PID,
                              Disease = clinical$Disease,
                              Timepoint = clinical$Timepoint,
                              Condition = clinical$condition)


boxplot <- boxplotdata_unpaired[,c("FCGR1CP",
                                   "PID",
                                   "Disease",
                                   "Timepoint",
                                   "Condition")]




colnames(boxplot)[1] <- "gene"

## GET P-VALUES ------------------------------------------

#create a list of length-2 vectors specifying the groups to be compared eg. list(c("A", "B"), c("B", "C"))

comparisonlist <- list()
for(i in 1:(nrow(nameconvert)-1)){
  x <- as.character(nameconvert[i,c("c1","c2")])
  comparisonlist[[i]] <- x
}



boxplot$Condition <- factor(boxplot$Condition, levels = c("HC_T0", "TB_T0", "TB_T2", "TB_T4", "TB_T6", "HC_T6"))

#### STEP 1) THIS IS  A TABLE OF THE COMPARISONS I ACTUALLY WANT ###
stat.table <- boxplot %>%
  t_test(gene ~ Condition,
         comparisons = comparisonlist)
stat.table<- stat.table %>%
  add_xy_position(x = "Timepoint", dodge = 0.8)
stat.table$contrast <- paste0(stat.table$group2, "-" ,stat.table$group1)


# # Can't figure out how to make p-values for grouped plot where only T0 and T6 are grouped
# #### step 2) THIS GIVES US X POSITIONS BUT TOO MANY COMPARISONS ###
# stat.table2 <- boxplot %>%
#   # group_by(Timepoint) %>%
#   t_test(gene ~ Condition)
# stat.table2<- stat.table2 %>%
#   add_xy_position(x = "Time", dodge = 0.8)
# stat.table2$contrast <- paste0(stat.table2$group1, "-" ,stat.table2$group2)
# 

#Non grouped 
#### step 2) THIS GIVES US X POSITIONS BUT TOO MANY COMPARISONS ###
stat.table2 <- boxplot %>%
  # group_by(Timepoint) %>%
  t_test(gene ~ Condition)
stat.table2<- stat.table2 %>%
  add_xy_position(x = "Condition", dodge = 0.8)
stat.table2$contrast <- paste0(stat.table2$group2, "-" ,stat.table2$group1)



#### filter STAT TABLE 2 TO ONLY INCLUDE COMPARISONS OF INTEREST####
stat.table3 <- stat.table2[match(stat.table$contrast, stat.table2$contrast),]

#Check that these are in the same order
gsub(" ", "", nameconvert$contrast[1:9]) == stat.table3$contrast
# one is same comparison just written different order

#### STEP 3) CREATE A TABLE WITH GROUPS NEEDED FOR STAT_P_VALUE MANUAL #####
stat.table4 <- cbind(stat.table3, resultsname = nameconvert$label[1:9])

# ##### manually change y positions - this is for aesthetics after seeing boxplot, all p-values are staggered but some can be next to eachother #####
stat.table4[which(stat.table4$resultsname == "contrast7"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"]
stat.table4[which(stat.table4$resultsname == "contrast9"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"]
stat.table4[which(stat.table4$resultsname == "contrast8"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"]
stat.table4[which(stat.table4$resultsname == "contrast4"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"]

#### STEP 4) REPLACE P VALUES FROM T-TEST WITH THOSE OUR DIFF EXP ANALYSIS  
for (i in stat.table4$resultsname) {
  ifelse("FCGR1CP" %in% listoftT[[i]]$genename,
         stat.table4[which(stat.table4$resultsname == i),"p"] <- listoftT[[i]]["FCGR1CP", "P.Value"],
         stat.table4[which(stat.table4$resultsname == i),"p"] <- "NA")
}
stat.table4$p <- signif(as.numeric(stat.table4$p), digits = 3)
stat.table4[which(is.na(stat.table4$p)), "p"] <- "NA" #the NA characters have become code NAs from the previous code, so reassign them to character NAs - surely there is a better way to do this ??, may as if !na statement in the ggplot line for showing p value




plotdisplay <-  ggplot(boxplot, aes(
  x = Condition,
  y = gene)) +
  
  theme_bw(base_size = 12)+
  
  # geom_boxplot( aes(fill = Disease))+
  
  
  
  geom_dotplot(binaxis = "y",
               stackdir='center',
               position=position_dodge(),
               stackratio = 0.3,
               aes(fill = Disease),
               dotsize = 0.55)+
  
  
  # facet_wrap(~ cse, strip.position = "bottom") +
  # 
  stat_summary(fun.data = "mean_se",
               geom = "errorbar",
               width = 0.2,
               aes(fill = Disease)) +
  
  labs(title = paste0("gene")) +
  ylab (label = "Normalised Expression") +
  xlab (label = "Group") +
  # theme(legend.position = "None")
  # +
  
  stat_pvalue_manual(stat.table4,
                     label = "p",
                     tip.length = 0.02,
                     # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                     #  bracket.shorten = 0.1,
                     size = 3)

print (plotdisplay)




# ================================================================================== #
## 2.3. HEATMAPS ===================================================================
# ================================================================================== #


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






