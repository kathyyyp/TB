#mRNA response in TB patients
# ================================================================================== #
# OPTIONAL. LOAD WORKSPACES ========================================================
# ================================================================================== #

# load("workspaces/16_02_25_auc.Rdata")


# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

#Mac
# my_directory <- "/Volumes/One Touch/RBMB"
# setwd(paste0(my_directory,"/TB"))
# .libPaths("/Volumes/One Touch/rlibrary")


#Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
setwd(file.path(my_directory,"TB"))
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")


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

data.dir <- file.path(my_directory,"TB", "data")




#  Main aim; primary proposed strategy is to identify the potential biomarkers for monitoring Levofloxacin treatment responses of MDR-TB patients as compared with household contacts.

# (i)	Distinguish individuals with multipledrug-resistant active tuberculosis (MDR-TB) and household contacts with latent TB infection (LTBI)
# (ii)	Predict treatment responses and posttreatment (either levofloxacin or placebo) residual inflammation of individuals who are MDR-TB patients and household contacts at 2, 4 and 6 months

# 1: To determine the differences in biomarkers at initiation and completion of treatment in the levofloxacin group compared to the placebo group 
# 2: To determine if a biomarker signature can be identified for patients with microbiologically proven TB compared with contacts that do or do not develop infection. 

# Household contacts received either levofloxacin or placebo tablets once per day for 6 months. 
# TB patients had blood collected at 4 time points (at diagnosis-T0 and T2, T4 and T6-months post therapy commencement)
# Household contacts had blood collected at time of inclusion in study and at 6 months).

#By T2, can we predict treatment outcome?

# Samples = whole blood
# PCR

#Household contacts can be considered healthy controls (T0 = 100 and T6 = 50)



# ================================================================================== #
# 1. LOAD IN DATA ==================================================================
# ================================================================================== #

#File Legend
#Plate 1:TB
#Plate 2:Household Contact
#Plate 3-4:Household Contact T0 and T6
#Plate 5-7: TB patients with treatment 2,4 and 6 months

# LOAD IN DATA ------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd(file.path(my_directory,"TB", "data","raw"))

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
HC_More_Info <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "HC infos", skip = 1, col_names = TRUE)
HC_samples_place <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "HC samples place", col_names = FALSE)

TB_More_Info <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "index infos", skip = 1, col_names = TRUE)
TB_samples_place <- read_excel("Fluidigm assay 400 samples patients infos 030225.xlsx", sheet = "index patients samples place", col_names = FALSE)

## Plate layout ------------------------------------------------------------------------------------------------------------------------------------------------
plate_layout <- read_excel("Fluidugm data analysis-ct values raw data 231224-new.xlsx", sheet = "plate layout", col_names = TRUE)

setwd(file.path(my_directory, "TB"))




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

library(tidyverse)

TB_Plate1_info <- as.data.frame(TB_Plate1_info)
TB_Plate1_info <- TB_Plate1_info[,-1]
TB_Plate1_info <- TB_Plate1_info %>% select(-c(EDTA,
                                               Catalogue...9,
                                               Location...10,
                                               Catalogue...12,
                                               Location...13)) #dont need info in these columns
duplicated(TB_Plate1_info$PID)

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

HC_Plate3.4_Info <- as.data.frame(HC_Plate3.4_Info)
HC_Plate3.4_Info <- HC_Plate3.4_Info[,-1]
HC_Plate3.4_Info <- HC_Plate3.4_Info %>% select(-c(EDTA,
                                                   Catalogue...9,
                                                   Location...10,
                                                   Catalogue...12,
                                                   Location...13,
                                                   ...14,
                                                   ...15,
                                                   ...16)) #dont need info in these columns

HC_More_Info <- as.data.frame(HC_More_Info)
HC_More_Info <- HC_More_Info[,-1]

HC_samples_place<- as.data.frame(HC_samples_place)
HC_samples_place<- HC_samples_place[,-1]
HC_samples_place<- HC_samples_place%>% select(c(1:6,10)) #dont need info in these columns

TB_Plate5.6.7_Info <- as.data.frame(TB_Plate5.6.7_Info)
TB_Plate5.6.7_Info <- TB_Plate5.6.7_Info[,-1]
TB_Plate5.6.7_Info <- TB_Plate5.6.7_Info %>% select(-c(EDTA,
                                                       Catalogue...9,
                                                       Location...10,
                                                       Catalogue...12,
                                                       Location...13,
                                                       ...14,
                                                       ...15,
                                                       ...16)) #dont need info in these columns
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

# listofpatientinfo <- list(TB_Plate1_info = TB_Plate1_info, HC_Plate2_info = HC_Plate2_info, HC_Plate3.4_Info = HC_Plate3.4_Info, TB_Plate5.6.7_Info = TB_Plate5.6.7_Info)

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
# NEED TO FIX THE DECIMAL PAXGENES
# -- EDIT: THE ONES WITH DECIMALS WERE REMOVED
samples_place <- as.data.frame(rbind(TB_samples_place[,1:7], HC_samples_place))
samples_place$PAXGENE[duplicated(samples_place$PAXGENE)] #there are 5 duplicates
samples_place_raw <- samples_place[-which(duplicated(samples_place$PAXGENE)),]
samples_place_raw$PAXGENE <- as.character(samples_place_raw$PAXGENE) 

#200 patients total in samples_place_raw
# If this line comes up with a warning, it's because one of the patients has a mislabelled in the excel sheet, written as a different timepoint!!!!
# FIxed that manually in excel

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
## 6.4. MERGE THE OLD AND NEW CLINICAL FILES & REMOVE PAXGENES  WITH DECIMALS ======
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

#pivot wide to see the timepoint samples for each patient
clinical_raw2_wide <- as.data.frame(pivot_wider(clinical_raw2[,c(1:4,6:ncol(clinical_raw2))], 
                                                names_from ="Timepoint", 
                                                values_from = "PAXGENE"))

length(unique(clinical_raw2$PID))
clinical_raw2_wide[which(duplicated(clinical_raw2_wide$PID)),]

# ================================================================================== #
## 6.8. FILL IN MISSING METADATA TO REMOVE DUPS ====================================
# ================================================================================== #

# From the wide table above, there are duplicate PID rows
# WHY? Some of the rows are duplicated because their duplicate has missing metadata
# clinical_excel2 has every paxgene and all the extra metadata (duplicated for each patient)
# when clinical_excel2 was merged with clinical_raw, patients that are in clinical_raw but not in clinical_excel2 all have NA for the extra metadata (eg if a patient has 4 rows for the 4 TB timepoints, all 4 rows have NAs for the hiv, concomitant etc columns)
# however for the samples that have the decimal paxgenes, when merged, some patients had duplicate rows for the same timepoint because they were different paxgenes (eg Patient A has T0, T6 and another T6. I removed the duplicate row that came from clinical_excel2 because it had the decimal paxgene and trusted that the paxgene from clinical_raw is right. That means for one timepoint T0, they had extra metadata from clinical_excel2, but for the other timepoint they didn't (because it was from clinical_raw). We when we pivot, there is a discrepancy because all the rows no longer match for the same patient. so we use this command to group by PID and fill in any missing rows with the row above it (info from the same patient)
# But now, because we removed the decimal PIDs, the same samples from old clinical file which replaced them, dont have province and sex data

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

#clinical_raw2 has 405 samples
#expression_raw has 384 samples

clinical_raw_wide <- pivot_wider(clinical_raw[,c(1:4,6:8)], names_from ="Timepoint", values_from = "PAXGENE")
clinical_excel2_wide <- pivot_wider(clinical_excel2[,-5], names_from ="Timepoint", values_from = "PAXGENE") 

unique(clinical_excel2_wide$PID)


write.csv(clinical_raw2_wide_nodecimal, file.path(data.dir, "processed", "clinical_raw2_wide_nodecimal.csv")) #NO DECIMALS
write.csv(clinical_raw_wide, file.path(data.dir, "processed", "clinical_raw_wide.csv")) #WITH DECIMALS
write.csv(clinical_excel2_wide, file.path(data.dir, "processed", "clinical_excel2_wide.csv")) #WITH DECIMALS


setdiff(unique(clinical_excel2_wide$PID),unique(clinical_raw_wide$PID))



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

write.csv(clinical, "data/processed/pre-QC/clinical.csv")
write.csv(expression, "data/processed/pre-QC/expression.csv")





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

write.csv(clinical_paired, "data/processed/pre-QC/clinical_paired.csv")
write.csv(expression_paired, "data/processed/pre-QC/expression_paired.csv")



