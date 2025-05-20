#TB vs Household Contact PCR app-7

#Mac
my_directory <- "/Volumes/One Touch/RBMB/TB"
setwd(file.path(my_directory,"shiny","data"))
.libPaths("/Volumes/One Touch/rlibrary")


# INSTRUCTIONS
# STEP 1: Change the line below to your directory (path to where the 'shiny' folder is located on your computer)
# my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB/TB"

# STEP 2: Click "Run App" in the top right of the source pane OR highlight all of the script after this line to the end of the page and press Ctrl + Enter to run

# setwd(file.path(my_directory,"shiny","data"))
# .libPaths(file.path(my_directory,"shiny","library"))


library(shiny)
library(readxl)
library(bslib)
library(DT)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(ggrepel)
library(gplots)
library(tidyverse)
library(rstatix)
library(limma)
library(shinydashboard)
library(GSVA)
library(ggfortify)
library(heatmap3)
# 
# remove.packages("shiny")
# install.packages("shiny", type = "source")
# 
# install.packages("readxl")
# install.packages("bslib")
# install.packages("DT")
# install.packages("ggplot2")
# install.packages("ggprism")
# install.packages("ggpubr")
# install.packages("ggrepel")
# install.packages("gplots")
# install.packages("tidyverse")
# install.packages("rstatix")
# install.packages("shinydashboard")
# BiocManager::install("GSVA")
# install.packages("ggfortify")
# install.packages("heatmap3")



# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #

#File Legend
#Plate 1: TB
#Plate 2: Household Contact
#Plate 3-4: Household Contact T0 and T6
#Plate 5-7: TB patients with treatment 2,4 and 6 months

## Plate 1: T0 TB patients  ---------------------------------------------------------------------------------------------------------------------------------
TB_Plate1_ct <- read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet= "Plate 1", skip=2, col_names = TRUE)
TB_Plate1_info <- read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "TB T0 PT", col_names = TRUE)

## Plate 2: T0 Household Contacts ---------------------------------------------------------------------------------------------------------------------------------
HC_Plate2_ct <- read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet= "Plate 2", skip=2, col_names = TRUE)
HC_Plate2_info <- read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "HC T0 PT", col_names = TRUE)

## Plate 3-4: T0 AND T6 Household Contacts ------------------------------------------------------------------------------------------------------------------------------------------
HC_Plate3_ct <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "Plate 3", skip = 1, col_names = TRUE)
HC_Plate4_ct <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "Plate 4", skip = 1, col_names = TRUE)
HC_Plate3.4_Info <-  read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "P3and P4 HC PT", col_names = TRUE)

## Plate 5-7: T0, T2, T4, T6 TB Patients  -------------------------------------------------------------------------------------------------------------------------------------------
TB_Plate5_ct <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "Plate 5", skip = 1, col_names = TRUE)
TB_Plate6_ct <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "Plate 6", skip = 1, col_names = TRUE)
TB_Plate7_ct <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "Plate 7", skip = 1, col_names = TRUE)
TB_Plate5.6.7_Info <-read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "P5 P7 TB PT", col_names = TRUE)

#Note plates 5, 6 and 7 are not T2, T4 and T6 respectively. They each contain various timepoints. Shown through code below
#xx = intersect(as.character(TB_Plate5.6.7_Info$PAXGENE), colnames(TB_Plate7_ct))
#TB_Plate5.6.7_Info[match(xx, TB_Plate5.6.7_Info$PAXGENE),]

# Additional clinical info
  #More_Info contains extra metadata (hiv status, concomitatn disease, weight, diabetes, smoking etc)
  #HC_samples_place just matches up the PIDs fom More_Info with corresponding PAXGENE ids and timepoints
HC_More_Info <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "HC infos", skip = 1, col_names = TRUE)
HC_samples_place <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "HC samples place", col_names = FALSE)

TB_More_Info <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "index infos", skip = 1, col_names = TRUE)
TB_samples_place <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "index patients samples place", col_names = FALSE)



# ================================================================================== #
# 2. CLEAN UP EXPRESSION DATA ======================================================
# ================================================================================== #

## Plate 1: T0 TB patients  ---------------------------------------------------------------------------------------------------------------------------------
#Remove empty row and column
TB_Plate1_ct <- as.data.frame(TB_Plate1_ct[-1,-ncol(TB_Plate1_ct)])

## Plate 2: T0 Household Contacts ---------------------------------------------------------------------------------------------------------------------------------
#Remove empty row and column
HC_Plate2_ct <- as.data.frame(HC_Plate2_ct[-1,-ncol(HC_Plate2_ct),])

## Plate 3-4: T0 AND T6 Household Contacts ------------------------------------------------------------------------------------------------------------------------------------------
HC_Plate3_ct <- HC_Plate3_ct[-c(1:2),]
HC_Plate4_ct <- HC_Plate4_ct[-1,]

## T0, T2, T4, T6: TB Patients  ------------------------------------------------------------------------------------------------------------------------------------------
TB_Plate5_ct <- TB_Plate5_ct[-c(1:2), -ncol(TB_Plate5_ct)]
TB_Plate6_ct <- TB_Plate6_ct[-c(1:2),]
TB_Plate7_ct <- TB_Plate7_ct[-1,]


## Make genes row names ----------------------------------------------------------------------------------------------------------------------------------------------------------
listoffiles <- list(TB_Plate1_ct = TB_Plate1_ct, HC_Plate2_ct = HC_Plate2_ct, HC_Plate3_ct= HC_Plate3_ct, HC_Plate4_ct = HC_Plate4_ct, TB_Plate5_ct= TB_Plate5_ct, TB_Plate6_ct = TB_Plate6_ct, TB_Plate7_ct = TB_Plate7_ct)
for (i in 1:length(listoffiles)){
  df <- listoffiles[[i]]
  df <- as.data.frame(df)
  row.names(df) <- df[,1]
  df <- df[,-1]
  listoffiles[[i]] <- df
}

#Create R objects from the list (note all of the assign() functions below are purely for visualising the data frames in RStudio)
for(i in 1:length(listoffiles)){
  assign(names(listoffiles)[i], listoffiles[[i]])
}

#3200512 weird values


# ================================================================================== #
# 3. REMOVE 18S AND SUBTRACT FROM 40 ===============================================
# ================================================================================== #


listoffiles2 <- list(TB_Plate1 = TB_Plate1_ct, HC_Plate2 = HC_Plate2_ct, HC_Plate3= HC_Plate3_ct, HC_Plate4 = HC_Plate4_ct, TB_Plate5= TB_Plate5_ct, TB_Plate6 = TB_Plate6_ct, TB_Plate7 = TB_Plate7_ct)
for (i in 1:length(listoffiles2)){
  df <- listoffiles2[[i]]
  df <- df[-which(row.names(df) == "18S VER1"),] #REMOVE as not enough variability in plates 3-7
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
any(is.na(listoffiles2[[7]]))

if(any(is.na(listoffiles2[[7]])) == TRUE){ 
  stop("NA values present in expression data files") }

# ================================================================================== #
# 4. COMBINE ALL DATA ==============================================================
# ================================================================================== #
listofdata_beforenorm <- c(listoffiles2)

##Reorder to ensure all matrices are in same order before binding them  ---------------------------------------------
geneorder <- row.names(listoffiles2[[1]])
for (i in 1:length(listofdata_beforenorm)){
  df <- listofdata_beforenorm[[i]]
  df <- df[geneorder,]
  
  ##For any plates with R in the sample name (it appears its plates 3, 5 and 6), remove the R (denotes that the sample is a repeat) ---------------------------------------------
  colnames(df) <- gsub("R", "", colnames(df)) 
  listofdata_beforenorm[[i]] <- df
}

#Create R objects from the list
for(i in 1:length(listofdata_beforenorm)){
  assign(names(listofdata_beforenorm)[i], listofdata_beforenorm[[i]])
}

any(is.na(listofdata_beforenorm[[7]]))

all_data <- cbind(TB_Plate1, HC_Plate2, HC_Plate3, HC_Plate4, TB_Plate5, TB_Plate6, TB_Plate7)



# ================================================================================== #
# 5. NORMALISE (CYCLICLOESS) =======================================================
# ================================================================================== #

#set weight of housekeeping genes as 50 and other genes as 1
weights = rep(1, times = nrow(all_data))
weights[row.names(all_data) == "GAPDH"] <- 50
weights[row.names(all_data) == "B2M"] <- 50

norm_data <- normalizeCyclicLoess(all_data, weights = weights) #other arguments, = span, iterations, method




# ================================================================================== #
# 6. MAKE CLINICAL FILE ============================================================
# ================================================================================== #

## Plate 1: T0 TB patients  ---------------------------------------------------------------------------------------------------------------------------------
TB_Plate1_info <- as.data.frame(TB_Plate1_info)
TB_Plate1_info <- TB_Plate1_info[,-1]
TB_Plate1_info <- TB_Plate1_info %>% select(-c(EDTA,
                                               Catalogue...9,
                                               Location...10,
                                               Catalogue...12,
                                               Location...13)) #dont need info in these columns
duplicated(TB_Plate1_info$PID)

## Plate 2: T0 Household Contacts ---------------------------------------------------------------------------------------------------------------------------------
HC_Plate2_info <- as.data.frame(HC_Plate2_info)
HC_Plate2_info$PID[duplicated(HC_Plate2_info$PID)] #??? there are two samples in HC_Plate2_info (3201412, 3101204) for this same patient at T0 (07-19-0043-01), but only 3201412 is in the HC_Plate_ct data file and most of the values are 999

#Removed 3101204 
HC_Plate2_info <- HC_Plate2_info[-which(HC_Plate2_info$PAXGENE == "3101204"),]
HC_Plate2_info <- HC_Plate2_info[,-1]
HC_Plate2_info <- HC_Plate2_info %>% select(-c(EDTA,
                                               Catalogue...9,
                                               Location...10,
                                               Catalogue...12,
                                               Location...13)) 

## Plate 3-4: T0 AND T6 Household Contacts ------------------------------------------------------------------------------------------------------------------------------------------
HC_Plate3.4_Info <- as.data.frame(HC_Plate3.4_Info)
HC_Plate3.4_Info <- HC_Plate3.4_Info[,-1]
HC_Plate3.4_Info <- HC_Plate3.4_Info %>% select(-c(EDTA,
                                                   Catalogue...9,
                                                   Location...10,
                                                   Catalogue...12,
                                                   Location...13,
                                                   ...14,
                                                   ...15,
                                                   ...16)) 

## Additional data for Household Contacts ------------------------------------------------------------------------------------------------------------------------------------------
HC_More_Info <- as.data.frame(HC_More_Info)
HC_More_Info <- HC_More_Info[,-1]

HC_samples_place<- as.data.frame(HC_samples_place)
HC_samples_place<- HC_samples_place[,-1]
HC_samples_place<- HC_samples_place%>% select(c(1:6,10)) 

## Plate 5-7: T0, T2, T4, T6 TB Patients  -------------------------------------------------------------------------------------------------------------------------------------------
TB_Plate5.6.7_Info <- as.data.frame(TB_Plate5.6.7_Info)
TB_Plate5.6.7_Info <- TB_Plate5.6.7_Info[,-1]
TB_Plate5.6.7_Info <- TB_Plate5.6.7_Info %>% select(-c(EDTA,
                                                       Catalogue...9,
                                                       Location...10,
                                                       Catalogue...12,
                                                       Location...13,
                                                       ...14,
                                                       ...15,
                                                       ...16)) 
duplicated(TB_Plate5.6.7_Info$PID)

TB_More_Info <- as.data.frame(TB_More_Info)
TB_More_Info <- TB_More_Info[-1,-1]

TB_samples_place<- as.data.frame(TB_samples_place)
TB_samples_place<- TB_samples_place[,-1]
TB_samples_place<- TB_samples_place%>% select(c(1,3,4,5,6,7,11,2)) #dont need info in these columns

#Update colnames FOR ALL CLINICAL
cols <- c("PID",
          "Disease",	
          "Initial",
          "Year.of.birth",
          "Date.collected",
          "Timepoint",
          "PAXGENE")

colnames(TB_Plate1_info) <- cols
colnames(HC_Plate2_info) <- cols
colnames(HC_Plate3.4_Info) <- cols
colnames(TB_Plate5.6.7_Info) <- cols
colnames(HC_samples_place) <- cols
colnames(TB_samples_place) <- c(cols, "status")

listofpatientinfo <- list(TB_Plate1_info = TB_Plate1_info, HC_Plate2_info = HC_Plate2_info, HC_Plate3.4_Info = HC_Plate3.4_Info, TB_Plate5.6.7_Info = TB_Plate5.6.7_Info)

clinical_raw <-as.data.frame(rbind(TB_Plate1_info, HC_Plate2_info, HC_Plate3.4_Info, TB_Plate5.6.7_Info))


#Remove duplicate sample 3204477
#This TB T0 sample has paired T2, T4 and T6 data and only PCR data in plate 7. But it also appears in Plate 1 clinical info. Remove the duplicate. when we match the paxgenes to the expression file later, it will be labelled as plate 7
clinical_raw <- clinical_raw[-which(duplicated(clinical_raw$PAXGENE)),]
row.names(clinical_raw) <- clinical_raw$PAXGENE




# ================================================================================== #
## 6.1. Get AGE.JOINED.STUDY ========================================================
# ================================================================================== #

library("stringr")
clinical_pid_ordered <- clinical_raw[order(clinical_raw$PID, clinical_raw$Timepoint),] #sort by PID and timepoint
length(unique(clinical_pid_ordered$PID)) #199
length(which(clinical_pid_ordered$Timepoint == "T0")) #199

# confirmation each patient has T0
clinical_pid_ordered$PID[which(clinical_pid_ordered$Timepoint == "T0")] == unique(clinical_pid_ordered$PID) 

# Subset only T0 patients, this has all 199 patients, because everyone has T0. Which is the date we want to calculate age from
clinical_pid_ordered_T0 <- clinical_pid_ordered[which(clinical_pid_ordered$Timepoint == "T0"),]

#Extrac year of T0 collection
year <- str_extract(string = as.character(clinical_pid_ordered_T0$Date.collected), pattern = "^(.+?)-")
year <-str_extract(year, pattern = ".*(?=-)")
clinical_pid_ordered_T0$year.collected <- year

#Age is calculated as T0 year of collection - Year that patient joined study
clinical_pid_ordered_T0$age <- as.numeric(clinical_pid_ordered_T0$year.collected) - as.numeric(clinical_pid_ordered_T0$Year.of.birth)

#Add these ages to the full clinical file, duplicated for the same patient
clinical_pid_ordered$age <- clinical_pid_ordered_T0$age[match(clinical_pid_ordered$PID, clinical_pid_ordered_T0$PID)]
#could also do age - date collected but then the same patient will have various ages which gets confusing



# ================================================================================== #
## 6.2. CLEAN UP CLINICAL FILE =====================================================
# ================================================================================== #

#Clean up file which now has age
clinical_raw <- clinical_pid_ordered

colnames(HC_More_Info)[1] <- c("PID")
colnames(HC_More_Info)[3] <- c("sex")
colnames(TB_More_Info)[1:4] <- c("PID", "province", "age", "sex")
colnames(TB_More_Info)[2] <- "province"

#Clean up colnames
colnames(HC_More_Info) <- gsub("_enrol","", colnames(HC_More_Info))
colnames(TB_More_Info) <- gsub("_crf101","", colnames(TB_More_Info))
colnames(TB_More_Info) <- gsub("_crf102","", colnames(TB_More_Info))
colnames(TB_More_Info)[which(colnames(TB_More_Info) == "smokelastmonth")] <- "smokedlastmonth"

#Other metadata that is missing or different in ove file vs the other
intersect(colnames(TB_More_Info), colnames(HC_More_Info))
setdiff(colnames(TB_More_Info), colnames(HC_More_Info)) 
setdiff(colnames(HC_More_Info), colnames(TB_More_Info)) 


clinical_additional_colnames <- intersect(colnames(TB_More_Info), colnames(HC_More_Info))
HC_Additional_Info <- HC_More_Info[,clinical_additional_colnames]
TB_Additional_Info <- TB_More_Info[,clinical_additional_colnames]
colnames(HC_Additional_Info) == colnames(TB_Additional_Info)



# ================================================================================== #
## 6.3. EXTRA CLINICAL INFORMATION  ================================================
# ================================================================================== #

## Create clinical_additional_features file, which contains extra clinical information about each patient
clinical_additional_features <- rbind(HC_Additional_Info, TB_Additional_Info)
row.names(clinical_additional_features) <- NULL

## Create samples_place file, which contains each patient's experimental conditions ie. their PAXGENE and timepoint

samples_place <- as.data.frame(rbind(TB_samples_place[,1:7], HC_samples_place))
samples_place$PAXGENE[duplicated(samples_place$PAXGENE)] #there are 5 duplicates
samples_place_raw <- samples_place[-which(duplicated(samples_place$PAXGENE)),]
samples_place_raw$PAXGENE <- as.character(samples_place_raw$PAXGENE) 

#Visualise patients and their associated paxgenes
samples_place_raw_wide <- pivot_wider(samples_place_raw[,-5], names_from = "Timepoint", values_from = "PAXGENE")
length(unique(samples_place_raw_wide$PID)) #200




# ================================================================================== #
## 6.4. FIXING UP PATIENTS SAMPLES MISMATCHED FILES  ===============================
# ================================================================================== #

# ## THE LAST FEW ARE DECIMALS????
# edit -- the ones with decimals were removed below
# samples_place_raw$PAXGENE <- gsub("\\..*","",samples_place_raw$PAXGENE) ##GET RID OF NUMBERS AFTER DECIMAL
# row.names(samples_place_raw) <- samples_place_raw$PAXGENE # this doesn't work becuase the paxgenes with decimals are still wrong after removing decimals. they are duplicates of other patients
# samples_place_raw$PAXGENE[duplicated(samples_place_raw$PAXGENE)]

# 201 in clinical_additional_features
clinical_additional_features$PID[duplicated(clinical_additional_features$PID)]
all(clinical_additional_features$PID %in% samples_place_raw$PID) #all true, so we have the same patients in both clinical info/experimental info files
# 09-13-0077 THIS SINGULAR PATIENT IS DUPLICATED in clinical_additional_features (province both an giang and lts, but an giang row is missing info)
# They have initials LTS and paxgene ids for T0 up to T6 (in samples_place_raw_wide)
# But also their province is LTS (in clinical_additional_features) 
# I think they are from an giang 
# Manually fixing this by keeping an giang in clinical_additional_features and removing LTS so now it matches with samples_place_raw
clinical_additional_features <- clinical_additional_features[-which(duplicated(clinical_additional_features$PID)),]

# ALL PATIENTS NOW MATCHING - NO DIFFERENCES
c(setdiff(samples_place_raw_wide$PID, clinical_additional_features$PID), setdiff(clinical_additional_features$PID, samples_place_raw_wide$PID))

## Combine clinical_additional_features and samples_place_raw (clinical + experimental info for 396 samples) -------------------------------------------------------
#make both files the same order so we can combine them
clinical_additional_features_longmatch <- clinical_additional_features[match(samples_place_raw$PID, clinical_additional_features$PID),]
clinical_additional_features_longmatch$PID == samples_place_raw$PID

clinical_excel2 <- cbind(samples_place_raw,clinical_additional_features_longmatch[,!colnames(clinical_additional_features_longmatch)=="PID"])

#### PAXGENES that are in the new clinical file but not the old one ------
# Most are the DECIMAL ONES - to fix
setdiff(clinical_excel2$PAXGENE,clinical_raw$PAXGENE)
nonmatching_PAXGENE_excel1vs2 <- clinical_excel2[match(setdiff(clinical_excel2$PAXGENE,clinical_raw$PAXGENE),clinical_excel2$PAXGENE), c("PID","PAXGENE")]


#Below line of code shows that there are 21 not in the old file
# BUT 14 of these are from patients that are in the old file, it' just the decimal paxgene that is the problem)
# 7 of them don't match because there are 6 PATIENTS (2 of the paxgenes are from same patient)that are in the new file but not the old file
# "10-11-0112"    "01-03-0138"    "09-13-0079"    "09-41-0011"    "07-19-0145"    "01-23-0039-03" "01-23-0039-03"
nonmatching_PAXGENE_excel1vs2$PID[!nonmatching_PAXGENE_excel1vs2$PID %in% clinical_raw$PID]

#this shows the 6 unmatched patients that are in new file and not the old file
setdiff(clinical_excel2$PID,clinical_raw$PID) 

### PAXGENES that are in the old clinical file but not the new one --------------------------------
#Note most of these are because of the decimals
nonmatching_PAXGENE_excel2vs1 <- clinical_raw[match(setdiff(clinical_raw$PAXGENE, clinical_excel2$PAXGENE),clinical_raw$PAXGENE), c("PID","PAXGENE")]
nonmatching_PAXGENE_excel2vs1$PID[!nonmatching_PAXGENE_excel2vs1$PID %in% clinical_excel2$PID]

#this shows 5 patients that are in the old file and not the new file
setdiff(clinical_raw$PID, clinical_excel2$PID) 



# ================================================================================== #
## 6.4a. MERGE THE OLD AND NEW CLINICAL FILES & REMOVE PAXGENES  WITH DECIMALS ======
# ================================================================================== #

# This contains all clinical samples in the original clinical files and in the new clinical files - even if not common between the two
clinical_raw2 <- merge(clinical_raw, clinical_excel2, by = colnames(clinical_raw)[1:7], all = TRUE) #all = TRUE means we keep all PIDs, not just those that are present in both objects

#07/02/24 - there are now 419 but that is because there are 14 samples with decimal points
#Remove the ones with the decimal points
clinical_raw2 <- clinical_raw2[-which(grepl("\\.", clinical_raw2$PAXGENE)),] #we have to use \\ to escape the dot, because a dot means 'any character' in R


length(unique(clinical_raw2$PID)) #205 unique


# ================================================================================== #
## 6.5. FIX UP THE AGE, SOME MISSING ===============================================
# ================================================================================== #

#from looking at data frame, the age is calculated is equal to the ages provided in additional features (except additional features has some not calculated yet so just remove it and keep my calculated ages)
# get ages for the ones that were only in new file (age.y), assign to age.x and then delete age.y
clinical_raw2[which(is.na(clinical_raw2$age.x)), "age.x"] <- round(clinical_raw2[which(is.na(clinical_raw2$age.x)), "age.y"], digits = 0)

#one of them is still the year so fix that
clinical_raw2[which(clinical_raw2$age.x == 1969), "age.x"] <- 2020-1969 #2020 was the date of collection for this patient

#delete age.y
clinical_raw2 <- clinical_raw2[,-which(colnames(clinical_raw2) =="age.y")]
colnames(clinical_raw2)[which(colnames(clinical_raw2) == "age.x")] <- "age"

any(duplicated(clinical_raw2$PAXGENE))
any(is.na(clinical_raw2$PAXGENE))


# ================================================================================== #
## 6.6. FIX UP SEX & SMOKING STATUS ================================================
# ================================================================================== #

unique(clinical_raw2$sex)
unique(clinical_raw2$smokedlastmonth)


#make sure any sample that does not have female or male is designated NA (not the character"NA")
clinical_raw2$sex <- tolower(clinical_raw2$sex)
clinical_raw2$sex <- ifelse(clinical_raw2$sex == "male", "male", 
                            ifelse(clinical_raw2$sex == "female", "female", NA))

#make sure smoking status only has three options
clinical_raw2$smoking_status <- tolower(clinical_raw2$smokedlastmonth)
clinical_raw2$smoking_status <- ifelse(clinical_raw2$smoking_status == "former", "former", 
                                       ifelse(clinical_raw2$smoking_status == "current", "current", "never"))



# ================================================================================== #
## 6.7. PIVOT WIDE TO VISUALISE UNIQUE PATIENTS ====================================
# ================================================================================== #
#dont need

# ================================================================================== #
## 6.8. FILL IN MISSING METADATA TO REMOVE DUPS ====================================
# ================================================================================== #

#Fill in the missing values
clinical_raw2 <- as.data.frame(
  clinical_raw2 %>%
    group_by(PID) %>% 
    fill(everything(), .direction = "downup") %>% #fill missing values for all columns within each PID (missing due to the decimal place discrepancy)
    ungroup() #ungroup data
) 

clinical_raw2_wide_nodecimal <- pivot_wider(clinical_raw2[,c(1:4,6:ncol(clinical_raw2))],  #col 5 is date-collection, so is not the same within each patient 
                                          names_from ="Timepoint", 
                                          values_from = "PAXGENE") 

duplicated(clinical_raw2_wide_nodecimal$PID)

row.names(clinical_raw2) <- clinical_raw2$PAXGENE

# ================================================================================== #
## 6.9. SAVE ALL VARIATIONS OF CLINICAL FILE (W AND W/O DECIMAL)  ==================
# ================================================================================== #
#dont need


# ================================================================================== #
# 7. MAKE EXPRESSION FILE  ==================
# ================================================================================== #

expression_raw <- norm_data

colnames(expression_raw)[duplicated(colnames(expression_raw), fromLast = TRUE)]
#3200266 this is a duplicate! There are two columns of data for this same paxgene
# Duplicate will get removed when we match up clinical_raw2 and expression_raw anyway, it will keep the first one. So don't need the line below
#fromLast = TRUE automatically takes the first one when matching in next steps(which should be the one from plate 3) 
# expression_raw <- expression_raw[, -duplicated(colnames(expression_raw), fromLast = TRUE)]

#The samples in clinical file that are not in expression file
data.frame(row.names(clinical_raw2)[(row.names(clinical_raw2) %in% colnames(expression_raw)) == FALSE])

#The samples in expression file that are not in clinical file
colnames(expression_raw)[(colnames(expression_raw)  %in%  row.names(clinical_raw2))== FALSE]



## Rearrange to same order as expression files -------------------------------------------------------------------
matches <- intersect(row.names(clinical_raw2), colnames(expression_raw))

clinical <- clinical_raw2[matches,]
expression <- expression_raw[,matches]

row.names(clinical) == colnames(expression)

# ================================================================================== #
# 8.  MORE ADJUSTMENTS TO CLINICAL FILE  ===========================================
# ================================================================================== #

#rename for ease
clinical[which(clinical$Disease == "Household contact"), "Disease"] <- "HC"
clinical[which(clinical$Disease == "Index patient"), "Disease"] <- "TB"

clinical$condition <- paste0(clinical$Disease, "_", clinical$Timepoint)
clinical$condition <- factor(clinical$condition)

#Include plate number to correct for batch effects
clinical$Plate <- "Plate"

for (i in 1:length(listofdata_beforenorm)){
  Plate <- listofdata_beforenorm[[i]] #get plates 1:7
  clinical[intersect(clinical$PAXGENE, colnames(Plate)), "Plate"] <- i #match the colnames from those plates with clinical paxgenes, get plate number
}

# ================================================================================== #
# 9.  SAVE CLINICAL AND EXPRESSION FILE  ===========================================
# ================================================================================== #
#dont need


# ================================================================================== #
# 10. GET PAIRED SAMPLES ONLY  ==========================================
# ================================================================================== #

#Exclude any patients in plate 1 and 2, which are those who ONLY had T0 measurements
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




# ================================================================================== #
# 2_TB_qc.R
# 11. REMOVE 2 OUTLIERS AND EDIT FILES ==============================================
# ================================================================================== #

# #REMOVE OUTLIERS THAT WERE IDENTIFIED IN PCA PLOT FROM CLINICAL AND EXPRESSION FILES
clinical <- clinical[!row.names(clinical) == "3204316",]
clinical <- clinical[!row.names(clinical) == "3201323",]

expression <- expression[,row.names(clinical)]

#Paired samples
unpaired_samples <- c(colnames(HC_Plate2_ct), colnames(TB_Plate1_ct))

#setdiff with clinical, which now has 2 outliers removed
clinical_paired <- clinical[c(setdiff(row.names(clinical),unpaired_samples)),] #setdiff is the opposite of intersect
clinical_paired_wide <- pivot_wider(clinical_paired[,-c(5,21,22)], names_from = Timepoint, values_from = PAXGENE)
#-22 because of the 3200266 patient who has data in different plates
duplicated(clinical_paired_wide$PID)

#These 4 patients are meant to have paired T0/T6 samples but are missing expression data for T0 or T6. Removed
PID_missingpair <- clinical_paired_wide[(rowSums(is.na(clinical_paired_wide[,c(18:21)])) >=3),"PID"]
PID_missingpair <- PID_missingpair$PID

clinical_paired <- clinical_paired[-c(match(PID_missingpair, clinical_paired$PID)),]
expression_paired <- expression[,row.names(clinical_paired)]







# LIMMA -----------------------------------------------------------------

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
    
  }
}



# SHINY ------------------------------------------------------------------------------------------------------------------
# UI ---------------------
ui <- dashboardPage(
  
  dashboardHeader(title = "TB vs HC Gene Expression"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("1) View input data", tabName = "upload", icon = icon("th")),
      menuItem("2) View normalised data", tabName = "viewdata", icon = icon("th")),
      menuItem("3) Principle Component Analysis", tabName = "pca", icon = icon("th")),
      menuItem("4) Differential Expression", tabName = "comparison", icon = icon("th")),
      menuItem("5) Gene Expression Boxplots", tabName = "boxplot", icon = icon("th")),
      menuItem("6) Signature Analysis", tabName = "signature", icon = icon("th")),
      menuItem("7) Validation of Signature", tabName = "signature", icon = icon("th"))

    ) #close sidebarMenu
  ), #closedashboardSidebar
  
  
  dashboardBody(
    title = "Household Contact vs TB Gene Expression",
    tabItems(
      
      
      ## UI 1) Upload Tab  ------------------------------------------------------------------------------
      tabItem("upload",
              tabsetPanel(
                tabPanel("Ct values",
                         
                         #Show PCR File -----------------#
                         h2("View PCR Data"),
                         selectInput("selectfile",
                                     label = "Select PCR data file to view",
                                     choices = names(listofdata_beforenorm)),
                         div(style = 'overflow-x: scroll', DTOutput("viewfile"))
    
                         
                ),#close tabPanel
                
                tabPanel("Clinical information",
                         
                         #Show Clinical file -----------------#
                         h2("View Clinical Info"),
                         selectInput("selectfile_clinical",
                                     label = "Select clinical info file to view",
                                     choices = c("TB_Plate1_info", 
                                     "HC_Plate2_info", 
                                     "HC_Plate3.4_Info",
                                     "TB_Plate5.6.7_Info")),
                         div(style = 'overflow-x: scroll',DTOutput("viewfile_clinical"))
                         
                )
              )
              
      ), #close TabItem
      
      
      
      ## UI 2) View normalised data Tab  ---------------------------------------------------------------------------------
      tabItem("viewdata",
              h2("View normalised expression data"),
              DTOutput("viewexpression"),
              h2("View clinical info master table"),
              div(style = 'overflow-x: scroll', DTOutput("viewclinical"))
      ), #close tabItem
      
      
      ## UI 3)
      tabItem("pca",
              h2("Principle Component Analysis"),
              selectInput("variable", label = "Variable of Interest",
                          choices = c("Plate", "Disease", "Timepoint", "condition")),
              plotOutput("pcaplot",  height="auto") #auto usually doesn't work due to how html/css is rendered. need to add to server 
              ),
              
      ## UI 4) Differential Expression tab  -------------------------------------------------------------------------------
      tabItem("comparison",

        selectInput("comparison", label = "Comparison of interest",
                    choices = c(nameconvert$contrast)),
      tabsetPanel(
      tabPanel("Differential Expression Results",
              h2("Volcano Plot"),
              plotOutput("volcano"),
              h2("Differential Expression Results Table"),
              downloadButton("downloadtable","Download Results Table"), #download table
              div(style = 'overflow-x: scroll',DT::DTOutput("deg")), #show table
              ),
      tabPanel("Heatmap",
               plotOutput("heatmap", height="1000px", width = "100%")
              )
      ) #close tabsetPanel
      ), #close tabItem
      
      ## UI 5) Boxplots tab  ---------------------------------------------------------------------------------------------
      tabItem("boxplot",
              selectizeInput("gene",
                             label = "Enter gene (HGNC symbol) of interest (case sensitive)",
                             choices = NULL,
                             width = "50%"),
              h4("Press backspace to search"),
              radioButtons("boxplotshow",
                           label = "Show:",
                           choices = c("TB and HC All timepoint", "Only T6 & T0", "Only TB T0 & HC T0", "Only TB")),
              
              plotOutput("boxplot",
                         width = "1000px",
                         height = "700px"),
              
              column(
                width = 8,
                fluidRow(
                  "Nominal P values from Mann-Whitney U test shown. Error bars show Mean + SE")
              ),
              
              
              column(width = 4,
                     downloadButton("downloadboxplot",
                                    "Download data for this gene")
              )
      ), #close tabItem
      
      ##UI 6) SIGNATURE GSVA -----------------------------------------------
      tabItem("signature",
              selectizeInput( "gene2",
                              label = "Enter gene (HGNC symbol) of interest (case sensitive)",
                              choices =  c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"),
                              multiple = TRUE,
                              width = "50%",
                              options = list(
                                'plugins' = list('remove_button')
                                # 'persist' = TRUE), 
                              )
                              
              ),
              plotOutput("boxplotgsva",
                         height = "500px")
              )

    ) #close tabItems
  )#close Dashboard body
) #close UI



#SERVER ----------------------------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  options(shiny.usecairo=TRUE)
  
  ## SERVER 1)  VIEW DATA ---------------------------------------------------------------------------------------------------------------
  
  #View PCR data file
  output$viewfile <-
    renderDT({
      datatable(listofdata_beforenorm[[input$selectfile]])
    })


  #View clinical file
  output$viewfile_clinical <-
    renderDT({
      datatable(listofpatientinfo[[input$selectfile_clinical]])
    })


  ## SERVER 2)   VIEW NORMALISED EXPRESSION AND MASTER CLINICAL ---------------------------------------------------------------------------------------------------------------

  # View expression & clinical file ------------#
   output$viewexpression <-  renderDT({
     expression
   })
  
  output$viewclinical <-  renderDT({
    clinical
  })
  
   
  ## SERVER 3) VIEW PCA PLOTS ---------------------------------------------------------------------------------------------------------------------
  output$pcaplot <- renderPlot({
    
    pca_res <- prcomp(t(as.data.frame(expression)), scale. = TRUE, center = TRUE) 
    
    clinical_pca <- cbind(pca_res$x, clinical)
    
    variable <- as.factor(clinical_pca[,input$variable])
    
    pcaplot <- pairs(clinical_pca[,c(1:5)],
          pch = 19, 
          cex = 1, 
          oma=c(3,3,3,15),
          col = variable)
    par(xpd=TRUE)
    legend("bottomright", fill = unique(variable), legend = c( levels(variable)))
    
    print(pcaplot)
  },
  height = function() {
    session$clientData$output_pcaplot_width
  })
  
   ## SERVER 4)  DIFFERENTIAL EXPRESSION AND VOLCANO PLOT --------------------------------------------------------------------------------------
   
  ### 4a: Volcano Plot --------------------------
   output$volcano<-renderPlot({
  
     c <- nameconvert[which(nameconvert$contrast == input$comparison),"label"]
     volcano <- ggplot(listoftT[[c]], aes(x = logFC, y = -log10(P.Value))) +
       geom_point(aes(color = legend)) +
       scale_color_manual(values = c("Downregulated" = "blue", "Not Sig" = "grey", "Upregulated" = "red"))+
       geom_hline(yintercept =-log10(max(as.data.frame(listoftT2[[c]][,"P.Value"]))),colour="black", linetype="dashed") +
       geom_vline(xintercept =-1,colour="black", linetype="dashed")+
       geom_vline(xintercept =1,colour="black", linetype="dashed")+
       geom_text_repel(data = as.data.frame(listoftT2[[c]]),
                       aes(label= genename),
                       size = 3, 
                       box.padding = unit(0.35, "lines"),
                       point.padding = unit(0.3, "lines") 
       ) +
       theme_bw(base_size = 12) +
       
       theme(legend.position = "bottom") +
       
       labs(title = input$comparison)
     
     print(volcano)
   })

   ### 4b: Show DEG table -------
   output$deg <- DT::renderDT({
     c <- nameconvert[which(nameconvert$contrast == input$comparison),"label"]
     datatable(listoftT[[c]], filter = list(position = "top", clear = TRUE))
   })
   
   ### 4c: Download Table --------
   output$downloadtable<-downloadHandler(
     filename = function() {
       paste(input$comparison, "_results.csv", sep = "")
     },
     content = function(file) {
       c <- nameconvert[which(nameconvert$contrast == input$comparison),"label"]
       write.csv(listoftT[[c]]
         , file, row.names = T)
     }
   )
  

  ## UPDATE SELECTIZEINPUT -----------------------------------------------------------------------------------------
  
  updateSelectizeInput(session, "gene", choices = row.names(expression), server = TRUE)
  updateSelectizeInput(session, "gene2", choices = row.names(expression), server = TRUE,
                       selected = c("IFITM1", "CD274", "TAP1", "GBP5", "GBP2", "S100A8", "FCGR1CP"))
  
  
  
  ## SERVER 5)  BOXPLOTS    ----------------------------------------------------------------------------------------------
  output$boxplot <-
    
    #Unpaired boxplots 
    renderPlot({
      
      boxplotdata_unpaired <- as.data.frame(t(expression))
      
      boxplotdata_unpaired <- boxplotdata_unpaired[row.names(clinical),]
      
      boxplotdata_unpaired <- cbind(boxplotdata_unpaired,
                                    PID = clinical$PID,
                                    Disease = clinical$Disease,
                                    Timepoint = clinical$Timepoint,
                                    Condition = clinical$condition)
      
      
      boxplot <- boxplotdata_unpaired[,c(input$gene,
                                         "PID",
                                         "Disease",
                                         "Timepoint",
                                         "Condition")]
      
      
      
      
      colnames(boxplot)[1] <- "gene"
      
      ### Get P-Values manually ------------------------------------------
      
      #create a list of length-2 vectors specifying the groups to be compared eg. list(c("A", "B"), c("B", "C"))
      
      comparisonlist <- list()
      for(i in 1:(nrow(nameconvert)-1)){
        x <- as.character(nameconvert[i,c("c1","c2")])
        comparisonlist[[i]] <- x
      }
      
      
      # boxplot$Condition <- factor(boxplot$Condition, levels = c("HC_T0", "TB_T0", "TB_T2", "TB_T4", "TB_T6", "HC_T6"))
      
      #### STEP 1) THIS IS  A TABLE OF THE COMPARISONS I ACTUALLY WANT ###
      stat.table <- boxplot %>%
        wilcox_test(gene ~ Condition,
               comparisons = comparisonlist)
      stat.table<- stat.table %>%
        add_xy_position(x = "Timepoint", dodge = 0.8)
      stat.table$contrast <- paste0(stat.table$group2, "-" ,stat.table$group1)
      

      #Non grouped 
      #### step 2) THIS GIVES US X POSITIONS BUT TOO MANY COMPARISONS ###
      stat.table2 <- boxplot %>%
        # group_by(Timepoint) %>%
        wilcox_test(gene ~ Condition)
      stat.table2<- stat.table2 %>%
        add_xy_position(x = "Condition", dodge = 0.8, scales = "free_y")
      stat.table2$contrast <- paste0(stat.table2$group2, "-" ,stat.table2$group1)
      
      
      
      #### filter STAT TABLE 2 TO ONLY INCLUDE COMPARISONS OF INTEREST
      stat.table3 <- stat.table2[match(stat.table$contrast, stat.table2$contrast),]
      
      #Check that these are in the same order
      gsub(" ", "", nameconvert$contrast[1:9]) == stat.table3$contrast
      # one is same comparison just written different order
      
      #### STEP 3) CREATE A TABLE WITH GROUPS NEEDED FOR STAT_P_VALUE MANUAL
      stat.table4 <- cbind(stat.table3, resultsname = nameconvert$label[1:9])
      
      ####  manually change y positions - this is for aesthetics after seeing boxplot, all p-values are staggered but some can be next to eachother
      stat.table4[which(stat.table4$resultsname == "contrast7"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast2"), "y.position"]
      stat.table4[which(stat.table4$resultsname == "contrast5"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"] + (stat.table4$y.position[1] - stat.table4$y.position[2])*2
      stat.table4[which(stat.table4$resultsname == "contrast9"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast2"), "y.position"]
      stat.table4[which(stat.table4$resultsname == "contrast8"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast2"), "y.position"]
      stat.table4[which(stat.table4$resultsname == "contrast4"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"] + stat.table4$y.position[1] - stat.table4$y.position[2]
      stat.table4[which(stat.table4$resultsname == "contrast3"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"]
      stat.table4[which(stat.table4$resultsname == "contrast6"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast1"), "y.position"] + (stat.table4$y.position[1] - stat.table4$y.position[2])*3
      

      #P value brackets are too far away from the highest point ?? why???
      #manually move down
      lowestbracket <- min(stat.table4$y.position)
      lowestbracketiwant <- max(boxplot$gene) + max(boxplot$gene)*0.02
      distancetomovedown <- lowestbracket - lowestbracketiwant
      stat.table4$y.position <- stat.table4$y.position - distancetomovedown
      
      #### STEP 4) REPLACE P VALUES FROM T-TEST WITH THOSE OUR DIFF EXP ANALYSIS  
      for (i in stat.table4$resultsname) {
        ifelse(input$gene %in% listoftT[[i]]$genename,
               stat.table4[which(stat.table4$resultsname == i),"p"] <- listoftT[[i]][input$gene, "P.Value"],
               stat.table4[which(stat.table4$resultsname == i),"p"] <- "NA")
      }
      stat.table4$p <- signif(as.numeric(stat.table4$p), digits = 3)
      stat.table4[which(is.na(stat.table4$p)), "p"] <- "NA" #the NA characters have become code NAs from the previous code, so reassign them to character NAs - surely there is a better way to do this ??, may as if !na statement in the ggplot line for showing p value
      
      
      ### CONDITIONAL - GENERATE BOXPLOT -----
      
      # universal theme
      # boxplot_theme <- theme(axis.title = element_text(size = 22),
      #                        axis.text = element_text(size = 22),
      #                        title = element_text(size = 18))


      
      if(input$boxplotshow == "Only TB T0 & HC T0"){
        boxplot <- boxplot[which(boxplot$Condition == "HC_T0" | boxplot$Condition == "TB_T0"),]
        stat.table4 <- stat.table4[which(stat.table4$contrast == "TB_T0-HC_T0"),]
        stat.table4$xmax <- 2
        stat.table4$y.position <- max(boxplot$gene)+max(boxplot$gene*0.02)

      plotdisplay <-  ggplot(boxplot, aes(
        x = Condition,
        y = gene)) +
        
        theme_bw(base_size = 20) +
        
        # geom_boxplot( aes(fill = Disease))+
        
      geom_dotplot(binaxis = "y",
                   stackdir='center',
                   position=position_dodge(),
                   stackratio = 0.3,
                   aes(fill = Disease),
                   dotsize = 0.55
                   )+

      # geom_jitter(position=position_jitter(0.2),
      #             alpha = 0.5,
      #             size = 3,
      #             aes(color = Disease))+

      
        # facet_wrap(~ cse, strip.position = "bottom") +
        # 
        stat_summary(fun.data = "mean_se",
                     geom = "errorbar",
                     width = 0.15,
                     linewidth = 0.7) +
        
        stat_summary(fun.y=mean,
                     geom="point", 
                     color="red",
                     size = 1.5) +
      
        labs(title = input$gene) +
        ylab (label = "Normalised Expression") +
        xlab (label = "Group") +
        # theme(legend.position = "None")
        # +
        
        stat_pvalue_manual(stat.table4,
                           label = "p",
                           tip.length = 0.02,
                           # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                           #  bracket.shorten = 0.1,
                           size = 4)
      
      print (plotdisplay)
      }
      
      
      
      if(input$boxplotshow == "Only T6 & T0"){
        boxplot <- boxplot[which(boxplot$Timepoint == "T6" | boxplot$Timepoint == "T0"),]
        boxplot$Condition <- factor(boxplot$Condition, levels = c("HC_T0", "TB_T0",  "HC_T6", "TB_T6"))
        stat.table4 <- stat.table4[which(stat.table4$contrast == "TB_T0-HC_T0" | stat.table4$contrast == "TB_T6-TB_T0" | stat.table4$contrast == "HC_T6-HC_T0" | stat.table4$contrast == "TB_T6-HC_T6"),]
        stat.table4$xmin <- c(1,1,2,3)
        stat.table4$xmax <- c(2,3,4,4)
        lowestbracket <- min(stat.table4$y.position)
        lowestbracketiwant <- max(boxplot$gene) + max(boxplot$gene)*0.02
        distancetomovedown <- lowestbracket - lowestbracketiwant
        stat.table4$y.position <- stat.table4$y.position - distancetomovedown
        
        
        plotdisplay <-  ggplot(boxplot, aes(
          x = Condition,
          y = gene)) +
          
          theme_bw(base_size = 20)+
          
          # geom_boxplot( aes(fill = Disease))+
          
         
          geom_dotplot(binaxis = "y",
                       stackdir='center',
                       position=position_dodge(),
                       stackratio = 0.3,
                       aes(fill = Disease),
                       dotsize = 0.55
          )+
          
          # geom_jitter(position=position_jitter(0.2),
          #             alpha = 0.5,
          #             size = 3,
          #             aes(color = Disease))+ 
          # facet_wrap(~ cse, strip.position = "bottom") +
          # 
          stat_summary(fun.data = "mean_se",
                       geom = "errorbar",
                       width = 0.15,
                       linewidth = 0.7) +
          
          stat_summary(fun.y=mean,
                       geom="point", 
                       color="red",
                       size = 1.5) +
          
          labs(title = input$gene) +
          ylab (label = "Normalised Expression") +
          xlab (label = "Group") +
          # theme(legend.position = "None")
          # +
          
          stat_pvalue_manual(stat.table4,
                             label = "p",
                             tip.length = 0.02,
                             # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                             #  bracket.shorten = 0.1,
                             size = 4)
        
        print (plotdisplay)
      }
      
      
      if(input$boxplotshow == "Only TB"){
        boxplot <- boxplot[which(boxplot$Disease == "TB"),]
        # boxplot$Condition <- factor(boxplot$Condition, levels = c("HC_T0", "TB_T0",  "HC_T6", "TB_T6"))
        stat.table4 <- stat.table4[which(stat.table4$contrast == "TB_T6-TB_T0" | stat.table4$contrast == "TB_T6-TB_T2" | stat.table4$contrast == "TB_T6-TB_T4" | stat.table4$contrast == "TB_T4-HC_T0" | stat.table4$contrast == "TB_T2-TB_T0" | stat.table4$contrast == "TB_T4-TB_T2"),]
        stat.table4$xmin <- c(1,2,1,3,2)
        stat.table4$xmax <- c(4,4,2,4,3)
        stat.table4[which(stat.table4$resultsname == "contrast5"), "y.position"] <- stat.table4[which(stat.table4$resultsname == "contrast3"), "y.position"] + (stat.table4$y.position[1] - stat.table4$y.position[3])
        
        plotdisplay <-  ggplot(boxplot, aes(
          x = Condition,
          y = gene)) +
          
          theme_bw(base_size = 20)+
          
          # geom_boxplot( aes(fill = Disease))+
          
          geom_dotplot(binaxis = "y",
                       stackdir='center',
                       position=position_dodge(),
                       stackratio = 0.3,
                       aes(fill = Disease),
                       dotsize = 0.55
          )+
          # geom_jitter(position=position_jitter(0.2),
          #             alpha = 0.5,
          #             size = 3,
          #             aes(color = Disease))+
          
          scale_fill_manual(values="#00BFC4") +
          
          # facet_wrap(~ cse, strip.position = "bottom") +
          # 
          stat_summary(fun.data = "mean_se",
                       geom = "errorbar",
                       width = 0.15,
                       linewidth = 0.7) +
          
          stat_summary(fun.y=mean,
                       geom="point", 
                       color="red",
                       size = 1.5) +
          
          labs(title = input$gene) +
          ylab (label = "Normalised Expression") +
          xlab (label = "Group") +
          # theme(legend.position = "None")
          # +
          
          stat_pvalue_manual(stat.table4,
                             label = "p",
                             tip.length = 0.02,
                             # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                             #  bracket.shorten = 0.1,
                             size = 4)
        
        print (plotdisplay)
      }
      
      
      else{
        plotdisplay <-  ggplot(boxplot, aes(
          x = Condition,
          y = gene)) +
          
          theme_bw(base_size = 20)+
          
          # geom_boxplot( aes(fill = Disease))+
          
          geom_dotplot(binaxis = "y",
                       stackdir='center',
                       position=position_dodge(),
                       stackratio = 0.3,
                       aes(fill = Disease),
                       dotsize = 0.55
          )+
          # geom_jitter(position=position_jitter(0.2),
          #             alpha = 0.5,
          #             size = 3,
          #             aes(color = Disease))+
          
          # facet_wrap(~ cse, strip.position = "bottom") +
          # 
          stat_summary(fun.data = "mean_se",
                       geom = "errorbar",
                       width = 0.15,
                       linewidth = 0.7) +
          
          stat_summary(fun.y=mean,
                       geom="point", 
                       color="red",
                       size = 1.5) +
          labs(title = input$gene) +
          ylab (label = "Normalised Expression") +
          xlab (label = "Group") +
          # theme(legend.position = "None")
          # +
          
          stat_pvalue_manual(stat.table4,
                             label = "p",
                             tip.length = 0.02,
                             # bracket.nudge.y =  c(0, 0.2, 0, 0.2),
                             #  bracket.shorten = 0.1,
                             size = 4)
        
        print (plotdisplay)
      }
      
    }, #Close render plot
    res = 120
    
    ) #Close render plot
      

  output$downloadboxplot<-downloadHandler(
    
    
    filename = function() {
      paste0(input$gene,"_data.csv")
    },
    content = function(file) { #'content' is a function that takes a single argument file that is a file path (string) of a nonexistent temp file, and writes the content to that file path. 
      
      
      
      boxplotdata <- as.data.frame(t(expression))
      
      plot <- cbind(boxplotdata,
                    disease = clinical$disease,
                    age = clinical$age)
      
      
      boxplot <- plot[,c(input$gene,
                         "disease",
                         "age")]
      
      write.csv(
        boxplot
        ,file, row.names = T)
    }
  )
  
  
  
  
  
  ## SERVER 6) GSVA Signature analysis BOXPLOT ----------------------------
  
  output$boxplotgsva <-   renderPlot({
    
    validate(need(input$gene2, 'Choose genes!')) #the message "Choose genes!" shows if no genes have been picked
    #req(input$gene) #can use this line too - will just show blank screen, no error message or chosen message shown
    #without these, error message shows up if no genes are picked
    
    
    gene_set_list <- list(c(input$gene2))
    
# gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))

gsvapar <- gsvaParam(as.matrix(expression), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway

all(row.names(gsva_res) == row.names(clinical))


boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = as.character(clinical$condition),
                                    PID = as.character(clinical$PID)))


gsva_theme <- theme(axis.title = element_text(size = 24),
                    axis.text = element_text(size = 24),
                    title = element_text(size = 20),
                    legend.position = "None") 


my_paired_comparisons <- list(
  c("HC_T0", "HC_T6"), 
  c("TB_T0", "TB_T2"), 
  c("TB_T0", "TB_T4"), 
  c("TB_T0", "TB_T6")
)

my_unpaired_comparisons <- list(
  c("HC_T0", "TB_T0"), 
  c("HC_T6", "TB_T6")
)


# The normalizeCyclicLoess() function in R (from the limma package) is a non-parametric normalization method that adjusts for systematic biases in microarray or high-throughput sequencing data. 
# It does not explicitly force the data to follow a normal distribution (Gaussian), but it equalizes intensity distributions across samples.
# Hence why our data may not be normally distributed still


## Normality tests ---------------------------------------------------------------------------------------
#dont need to show here

boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)
boxplot_gsva$group <- factor(boxplot_gsva$group)


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
stat.table.gsva <- map_dfr(my_paired_comparisons, function(groups) {
  #1) extract the timepoint name
  g1 <- groups[1]
  g2 <- groups[2]
  
  #2) make a data frame only containing our paired samples, making use of the function we made above
  df <- get_paired_data(data = boxplot_gsva, 
                        group1 = g1, 
                        group2 = g2)
  
df$group <- as.character(df$group) #don't know why this is needed. HC_T0 vs HC_T6 works without this. but the rest don't. something to do with unused factor levels?
  
#3) run wilcox_test on the paired samples
  if (n_distinct(df$PID) > 0) {
    wilcox_test(df, V1 ~ group, paired = TRUE) %>%
        add_xy_position(x = "group") %>% 
      mutate(group1 = g1, group2 = g2)
  } else {
    tibble()  # return empty if no paired data
  }
})


## Assign xmin and xmax positions manually  ---------------------------------------------------------------------------------------
stat.table.gsva$xmin = c(1,3,3,3)
stat.table.gsva$xmax = c(2,4,5,6)

# ================================================================================== #
##  UNPAIRED COMPARISONS  ===================================================================
# ================================================================================== #
stat.table.gsva2 <- wilcox_test(boxplot_gsva, V1 ~ group, 
                                paired = FALSE,
                                comparisons = my_unpaired_comparisons) %>%
        add_xy_position(x = "group")
stat.table.gsva2 <- stat.table.gsva2[,!colnames(stat.table.gsva2) %in%
                                          c("p.adj","p.adj.signif")]

stat.table.gsva.all <- rbind(stat.table.gsva, stat.table.gsva2)

# stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
lowest_bracket <- max(boxplot_gsva$V1) + 0.05*(max(boxplot_gsva$V1))
stat.table.gsva.all$y.position <- seq(lowest_bracket, by= 0.1, length.out = nrow(stat.table.gsva.all))


## Boxplot ggplot2 ---------------------------------------------------------------------------------------
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = factor(group),
  # x = factor(group),
  y = as.numeric(boxplot_gsva[,1]))) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(aes(color = group),position = position_dodge(1)) +
  
  #For unpaired boxplots
  # geom_jitter(aes(color = group),
  #             alpha = 0.5,
  #             size = 2.5,
  #             width = 0.3) +
  
  #For paired boxplots
  geom_point(aes(color = group))+
  geom_line(aes(group = PID), color = "black", alpha = 0.2) +
  
  stat_pvalue_manual(stat.table.gsva.all,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+
  
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  

  labs(title = paste0("Signature Analysis"),
       caption = paste0(paste("Signature:",paste(c(input$gene2), collapse = " "),"\n"),
                       "Wilcoxin rank-sum test performed for paired samples\n",
                       "Mann-Whitney-U performed for HC vs TB comparisons")) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition")
    
    print (boxplotfinal2)
  })
      



  output$heatmap <-

    renderPlot({
      
      contrastindex <- which(nameconvert$contrast == input$comparison)
      a <- nameconvert[contrastindex, "c1"]
      b <- nameconvert[contrastindex, "c2"]
      
      # heatmapsamples <- clinical$PAXGENE[which(clinical$condition == "TB_T0" | clinical$condition == "HC_T0")]
      heatmapsamples <- clinical$PAXGENE[which(clinical$condition == a | clinical$condition == b)]
      
      clinical_heatmap <- clinical[as.character(heatmapsamples),]
      expression_heatmap <- expression[,as.character(heatmapsamples)]
      
      clinical_heatmap_ordered <- clinical_heatmap[order(clinical_heatmap$Timepoint,clinical_heatmap$Disease),]
      
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
      
      labs = cbind(Timepoint,Disease)
      
      
      #reorder expression columns/samples to be the same as clinical_ordered
      heatmapdata <- as.matrix(expression_heatmap[,row.names(clinical_heatmap_ordered)])
      heatmapdata <- heatmapdata[-which(row.names(heatmapdata) == "B2M" | row.names(heatmapdata) == "GAPDH"),]
      

      
     heatmapplot <-  heatmap3(heatmapdata, 
               Colv=NA, 
               scale = "row",
               balanceColor=T,
               labCol=NA,
               showColDendro = F, 
               showRowDendro = F,
               margins=c(1,1),
               ColSideLabs = F, 
               ColSideColors =labs,
               cexRow=1.5,
               legendfun=function()
                 showLegend(legend=c("HC","TB","T0", "T2", "T4", "T6"),
                            col=c("seagreen","red3","lightgrey", "lightblue","skyblue3", "steelblue4"),
                            cex=1.5))
      print(heatmapplot)
    }) #close renderPlot

} #close Server


  
      

  

#RUN ----------------------------------------------------------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)



