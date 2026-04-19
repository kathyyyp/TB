# SCRIPT 2: Generating TB signature using ML 
# Method used to generate signature scores: GSVA
# ================================================================================== #
# A. SCRIPT SET UP =================================================================
# ================================================================================== #

# Windows
my_directory <- "C:/Users/165861_admin/OneDrive - UTS/Documents/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("C:/Users/165861_admin/OneDrive - UTS/rlibrary")

# Mac
my_directory <- "/Volumes/One Touch/RBMB"
main.dir <- file.path(my_directory, "TB")
setwd(file.path(main.dir))
.libPaths("/Volumes/One Touch/RBMB/rlibrary")

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
# ================================================================================== #
# B. SET UP DIRECTORY & OUTPUT PATHS ===============================================
# ================================================================================== #
output.dir <- file.path(main.dir, "output_v2")
signature.dir <- file.path(output.dir, "2_signature_generation")
if(!exists(signature.dir)) dir.create(signature.dir)

#For 2. Permutation Importance
current.roc.dir <- signature.dir
current.importance.dir <- file.path(current.roc.dir, "permutation_importance", "train_ML_16")
if(!exists(current.importance.dir)) dir.create(current.importance.dir, recursive = TRUE)

#For 3. GSVA
current.gsva.dir <- file.path(current.roc.dir, "GSVA", "ML_16")
if(!exists(current.gsva.dir)) dir.create(current.gsva.dir, recursive = TRUE)

# ================================================================================== #
# 1. LOAD IN DATA & GET UPREGULATED GENES ==========================================================
# ================================================================================== #
clinical <- read.csv(file.path(main.dir, "data/processed/post-QC/clinical.csv"), row.names = 1)
expression <- read.csv(file.path(main.dir, "data/processed/post-QC/expression.csv"), row.names = 1, check.names = FALSE)

# Add treatment outcome data #provided by Tess 10/03
TB_extraclinical <- read.csv(file.path(main.dir, "data", "raw","202603", "20250603_index patient biomarker with test set.csv"))
TB_validation_extraclinical <- read.csv(file.path(main.dir, "data", "raw","202603", "20250603_validation index patients with test set.csv"))
HC_extraclinical <- read_excel(file.path(main.dir, "data", "raw","202603", "20250608_hhc biomarkers for test set.xlsx"))
HC_rna_clinical <- read.csv(file.path(main.dir, "data", "raw","202603", "20251208_rna_household contacts.csv"))
TB_rna_clinical <- read.csv(file.path(main.dir, "data", "raw","202603", "20260108_rna_index data with tx outcomes.csv"))
TB_More_Info <- read_excel(file.path(main.dir, "data", "raw","Fluidigm assay 400 samples patients infos 030225.xlsx"), sheet = "index infos", skip = 1, col_names = TRUE)
TB_More_Info <- as.data.frame(TB_More_Info)
TB_More_Info <- TB_More_Info[-1,-1]

TB_extraclinical$pid %in% clinical$PID #41/105 TRUE
TB_validation_extraclinical$pid %in% clinical$PID #2/41 TRUE
HC_extraclinical$pid  %in% clinical$PID #33/99 TRUE
HC_rna_clinical$pid  %in% clinical$PID#1/100 TRUE
TB_rna_clinical$pid  %in% clinical$PID #ALL 113 FALSE

TB_More_Info$pid
colnames(TB_More_Info)[which(colnames(TB_More_Info) == "treatout_crf102")] <- "treatment_outcome"

clinical_TBonly <- clinical[which(clinical$Disease == "TB"),]

clinical_TBonly[which(!clinical_TBonly$PID %in% 
                        c( TB_extraclinical$pid, 
                           TB_validation_extraclinical$pid,
                           TB_More_Info$pid)),
                "PID"]
# [1] "07-14-0175" "10-11-0092" are in clinical but not the additional outcomes info



all_outcomes <- do.call(rbind, lapply(list(
  TB_extraclinical,
  TB_validation_extraclinical,
  TB_More_Info
), function(df) df[, c("pid", "treatment_outcome")]))

# Keep unique pid-outcome pairs
all_outcomes <- as.data.frame(all_outcomes)
row.names(all_outcomes) <- NULL
colnames(all_outcomes)[which(colnames(all_outcomes) == "pid")] <- "PID"
all_outcomes <- all_outcomes[match(unique(all_outcomes$PID), all_outcomes$PID),]
row.names(all_outcomes) <- NULL

# Merge into clinical_TBonly
clinical_TBonly <- merge(clinical_TBonly, all_outcomes, by = "PID", all.x = TRUE)

#Sent to tess to check if these patients have info
# clinical_nooutcomeinfo <- clinical_TBonly[which(is.na(clinical_TBonly$treatment_outcome)),]
# clinical_nooutcomeinfo <-clinical_nooutcomeinfo[match(unique(clinical_nooutcomeinfo$PID),clinical_nooutcomeinfo$PID),]
# write.csv(clinical_nooutcomeinfo, file.path(main.dir, "data", "processed", "TB_fluidigmdata_missing_tx_outcome.csv"))

table(clinical_TBonly$treatment_outcome)

clinical$treatment_outcome <- all_outcomes[match( clinical$PID, all_outcomes$PID), "treatment_outcome"]






clinical_wide <- as.data.frame(pivot_wider(clinical[,c("PID", "Disease", "Timepoint", "PAXGENE","treatment_outcome")], 
                                                names_from ="Timepoint", 
                                                values_from = "PAXGENE"))




### Get upregulated genes ############
# Get all 16 upregulated genes from our TB_T0 - HC_T0 differential expression analysis
diffexp_TBT0_vs_HCT0 <- read.csv(file.path(output.dir, "results", "tT2", "contrast1.csv"))
T0_upregulated <- diffexp_TBT0_vs_HCT0[which(diffexp_TBT0_vs_HCT0$legend == "Upregulated"), "genename"]

gene_set_of_interest <- T0_upregulated 
# Note: I also ran the ML model below to try and order all 46 genes by importance, but using the top 16 worked better

Disease <- clinical$Disease

#note 'expression' was already normalised by normalisecyclicloess
expression_auc <- as.data.frame(cbind(t(expression), Disease = Disease))#rows are samples instead of genes here

#Make Disease a factor
expression_auc[,"Disease"] <- factor(expression_auc[,"Disease"])

#Make gene expression numeric
expression_auc[,1:(ncol(expression_auc)-1)] <- sapply(expression_auc[,1:(ncol(expression_auc)-1)], as.numeric) 

#Subset
row.names(expression_auc) == row.names(clinical)
T0_samples <- row.names(clinical)[which(clinical$condition == "HC_T0" | clinical$condition == "TB_T0")]
expression_auc <- expression_auc[T0_samples,]



# ================================================================================== #
## Summary of ML ##
# use Repeated K-Fold Cross-Validation (5 folds, 5 repeats) to:
# Assess performance -  Outer CV loop
# split data into 5 folds, train with 4 folds and test with 1
#  Tune RF hyperparameters - inner CV loop
# For each training set in the outer loop, another 5-fold CV is run inside to find the best hyperparamters (mtry, splitrule and min.node.size)

#Then:
# Train RF on full training set
# Once the best hyperparameters are found, train RF on full training set (4 folds) with best hyperparameters
# Use trained RF to predict Disease on test fold	


# Shuffle data, repeat the whole process 5 times	

# ================================================================================== #




# ================================================================================== #
# 2. ORDER GENES BY PERMUTATION IMPORTANCE =========================================
# ================================================================================== #
## Define Repeated K-Fold CV Parameters 
# K-fold cross validation divides the data into k folds and trains the model k times, each time leaving out a different fold as the validation set (eg. If you have 50 samples, you can divide the data into 5 folds with each containing 10 samples. Then you can take 4 folds for training set and 1 fold for validation set (note validation set is not the same as testing set. both are used to evaluate a model's performanc but a validation set is used during the training process to fine-tune hyperparameters and select the best model, while a testing set is used after training to assess the final performance on completely unseen data, without making any further adjustments to the model. 
num_folds <- 5   # Number of CV folds
repeats <- 5      # Number of repetitions for performance estimation

## Define the hyperparameter grid and train control ------------------------------------
# Hyperparameters are the parameters that are explicitly defined to control the learning process before applying a machine-learning algorithm to a dataset
# These remain global since the grid and control methods are consistent across folds
rf_grid <- expand.grid(
  mtry = seq(1, 16), #This is a hyperparameter that controls the number of features considered for splitting a node when building each tree in the random forest (in a random forest, each tree is built with a random subset of the data and features) #Default is root of p, where p = number of features. Iris dataset has 4 features (sepal length, sepal width, petal length, petal width) (p=4), therefore the number of features considered at each node of the decision tree in random forest is 2.
  splitrule = "gini", #The Gini impurity is a metric used to measure how "impure" a node is in a decision tree (similar to entropy). A node with 100% of one class will have impurity score of 0, whereas a node with 50% of two classes (eg. 50% sepal width and 50% sepal length), will have a score of 0.5 because there is 50% chance of classifying a sample as one or the other
  min.node.size = c(1, 5, 10, 20, 50) # A hyperparameter that defines minimum number of samples (data points) required in a node before it can be split further. #Includes smaller values for finer splits and larger values to prevent overfitting
)

## For the inner tuning (specifiying parameters) on the training fold:
# Choose an appropriate number (e.g., 5 or 10) for inner CV folds.
inner_cv_control <- caret::trainControl(method = "cv", #method=cv means we choose k-fold cross validation.  
                                        number = 5, #number specifies number of resampling iterations or folds so 10 means data is split into 10 folds and will be trained 5 times 
                                        allowParallel = FALSE) #enables parallel computing during the resampling process which can significantly speed up the training process by utilizing multiple CPU cores.

# Initialize a list to store importance results from each fold
all_fold_importance <- list()

base_seed = 7

### Global Hyperparameter Tuning and Feature Importance Extraction 

# Example outer loop for repeats (we have 'repeats' defined)
overall_start_time = Sys.time()
for (r in seq_len(repeats)) {
  start_time = Sys.time()
  
  # Set a unique seed for each repeat
  set.seed(base_seed + r) #Note that set.seed only works for the next random draw? (i think), if you have another random draw, then you have to set seed again
  cat(sprintf("Repeat %d: Seed set to %d\n", r, base_seed + r))
  
  # Example creation of folds for this repetition
  # 'cv_folds' should be a list of indices defining each fold’s test set
  cv_folds <-
    caret::createFolds(expression_auc$Disease, k = num_folds, list = TRUE)
  
  for (fold_idx in seq_along(cv_folds)) {
    #There are 5 folds. Each fold has 38 or 39 samples (because we have 192 samples total and 192/5 = 38.5 samples per fold). Fold_idx takes one of these sets of samples
    
    cat(
      sprintf( "Repeat %d, Fold %d: Preparing training and test data...\n", r, fold_idx)
    )
    
    ### Define train/test indices for this fold
    test_indices <- cv_folds[[fold_idx]] #Indices for 38 or 39 samples (1 fold worth of samples)
    train_indices <- setdiff(seq_len(nrow(expression_auc)), test_indices)  #The remaining 154 or 153 samples that make up the training set (4 folds worth of samples)
    
    train_data <- expression_auc[train_indices,] #Subset expression data for training data indices Here, we have trained the model with ALL of the genes
    test_data <- expression_auc[test_indices,] #Subset testing data (unseen data for the model)
    
    ### Prepare the data frame for Random Forest tuning:
    # Here, we select the 16 genes. The model is trained on gene expression data from those 16 genes to predict Disease
    # Later we will assess permutation importance of those 16 genes 
    rf_train_data <-train_data[, c(gene_set_of_interest, "Disease")] 
    
    cat
    (sprintf("Repeat %d, Fold %d: Running hyperparameter tuning on training data only...\n",r,fold_idx)
    )
    
    ### Run caret::train multiple random forest (RF) models using different hyperparameters (from rf_grid) using only the training data of this fold
    rf_cv_fold <- caret::train(
      Disease ~ ., # This is a shorthand notation in R that means "all other variables" (which is all 16 genes) in the data frame, excluding the dependent variable. We want to use all genes to as features to predict disease
      data = rf_train_data,
      method = "ranger",
      tuneGrid = rf_grid,
      trControl = inner_cv_control,
      num.trees = 500 #Random forest parameter
    )
    
    ## Extract the best parameters (mtry, splitrule and min.node.size) for this fold (remember the 'best' is different for each fold and in each repetition)
    fold_best_params <- rf_cv_fold$bestTune
    
    cat(
      sprintf(
        "Repeat %d, Fold %d: Best Params found: mtry=%d, splitrule=%s, min.node.size=%d\n",
        r,
        fold_idx,
        fold_best_params$mtry,
        fold_best_params$splitrule,
        fold_best_params$min.node.size
      )
    )
    
    ##FINAL Define hyperparameter grid - using best parameters
    final_tuneGrid <- expand.grid(
      mtry = fold_best_params$mtry,
      splitrule = fold_best_params$splitrule,
      min.node.size = fold_best_params$min.node.size
    )
    
    cat( 
      sprintf("Repeat %d, Fold %d: Training final model with best parameters to extract importance...\n", r,fold_idx)
    )
    
    #Train final RF model on the training set for each fold using the best hyperparameters and features importance is extracted using permutation importance
    final_model_fold <- caret::train(
      Disease ~ .,
      data = rf_train_data,
      method = "ranger",
      trControl = caret::trainControl(method = "none", allowParallel = TRUE),
      # No CV here, just fit once on the full training set
      tuneGrid = final_tuneGrid,
      num.trees = 500,
      importance = "permutation" 
      ## Feature importance is assessed using permutation importance = assess how much the model's prediction performance (accuracy) would degrade if the expression values of each gene were randomly shuffled. if shuffling = significant decrease in model accuracy, then that gene must be improtant
    )
    
    # Extract feature importance for this fold #ie. how important is this feature in helping us differentiate HC and TB
    # 'Importance'? = Importance value of 0.1 means if we shuffled th expression values for that gene, the accuracy of the model would decrease by 10%, so the gene is somewhat important
    fold_importance <-
      final_model_fold$finalModel$variable.importance
    
    # Sort features by importance
    fold_importance <- sort(fold_importance, decreasing = TRUE) 
    
    # Create a rank vector: the first feature in fold_importance gets rank 1, next gets rank 2, etc.
    feature_ranks <- seq_along(fold_importance)
    
    # Predict on the test data
    test_predictions <-
      predict(final_model_fold, newdata = test_data) #note, even though test_data contains 40 genes, perm importance is still only calculated on 16 genes because the model was only trained on those features
    # Calculate overall accuracy on the test set
    test_accuracy <- mean(test_predictions == test_data$Disease) #this tells us how well our trained model predicts disease. test_predictions = the model's predicted classifications of "HB" and "TC"
    
    # Store fold-specific importance in a data frame with fold/repeat info
    fold_importance_df <- data.frame(
      Feature = names(fold_importance),
      Importance = as.numeric(fold_importance),
      Rank = feature_ranks, #remember feature_ranks = genes ordered by fold_importance (most important to least important)
      Test_Accuracy = test_accuracy,
      Repeat = r,
      Fold = fold_idx,
      stringsAsFactors = FALSE
    )
    
    # Append to the global list of all importance results
    all_fold_importance[[paste0("Repeat_", r, "_Fold_", fold_idx)]] <-
      fold_importance_df
    
    cat(
      sprintf("Repeat %d, Fold %d: Importance extracted and stored.\n",r,fold_idx )
    )
    
  } # end of fold loop
  end_time <- Sys.time()
  print(paste("Started at", overall_start_time, "Repeat", r, "Started at", start_time, "& finished at", end_time))
} # end of repeats loop
cat("All repeats and folds completed. Aggregating importance results...\n")

# Combine all fold-level importance data into one data frame
all_fold_importance_df <- do.call(rbind, all_fold_importance)

# Save the aggregated fold-level importance table if desired
write.csv(
  all_fold_importance_df,
  file = file.path(
    current.importance.dir,
    "02_fold_level_importance_results_genes.csv"
  )
)

cat(
  "Fold-level importance extraction complete. You can now proceed with frequency-based consensus or other methods to select top stable features.\n"
)


get_top_features <- function(feature_table, num_features) {
  # Check if requested number of features is within range
  max_features <- length(unique(feature_table$Feature)) #in this case length unique is 16 when feature_table = all_fold_importance_df
  if (num_features > max_features) { #step to check that the num_features wanted can actually be extracted. 
    message(paste("Requested number of features exceeds max (", max_features, "). Using max.", sep = ""))
    num_features <- max_features
  } else if (num_features < 1) {
    message("Requested < 1 feature. Using 1 feature.")
    num_features <- 1
  }
  
  # This is only relevant when we test 46 genes, when using the 16 genes, all of them are in the top20 so just use median
  # Calculate the frequency of each feature (gene) appearing in the top 20
  feature_frequency <- aggregate(Rank ~ Feature, data = feature_table[feature_table$Rank <= 20, ], FUN = length)
  colnames(feature_frequency)[2] <- "FrequencyInTop20" 
  
  # Calculate the median rank for each feature
  median_ranks <- aggregate(Rank ~ Feature, data = feature_table, FUN = median) #Remember, Rank in feature_table is the ranks for each gene accrding to importance. Here we take the median rank for each feature/gene
  colnames(median_ranks)[2] <- "MedianRank"
  
  # Merge frequency and median rank into a single data frame
  feature_stats <- merge(feature_frequency, median_ranks, by = "Feature", all.y = TRUE)

  # Replace na frequency (never in top 20) with 0
  feature_stats$FrequencyInTop10[is.na(feature_stats$FrequencyInTop20)] <- 0
  
  # Sort by frequency (descending) and then by median rank (ascending)
  feature_stats <- feature_stats[order(-feature_stats$FrequencyInTop20, feature_stats$MedianRank), ]
  
  # Select the top features based on the combined ranking
  top_features <- head(feature_stats$Feature, num_features) #manually saved this table into  "02_fold_level_importance_results_genes.csv" 05/03/26
  
  return(list(top_features = top_features,
              feature_stats = feature_stats))
  
}

num_features_to_use <- 16  # Specify how many features you want
top_features_cv <- get_top_features(feature_table = all_fold_importance_df, num_features = num_features_to_use)$top_features
feature_stats_cv <- get_top_features(feature_table = all_fold_importance_df,  num_features = num_features_to_use)$feature_stats

print(top_features_cv)

write.csv(top_features_cv, file = file.path(current.importance.dir, "top_features_16.csv"))
write.csv(feature_stats_cv, file = file.path(current.importance.dir, "summarised_feature_rankings.csv"))

#save note for myself
note <- paste0(
  "5 Repeat 5 CV RF Model trained on all 16 upregulated genes, assessed permutation importance of all 16 genes.\n",
  str_c(top_features_cv, collapse = " "))
write(note, file = file.path(current.importance.dir,"description_of_this_run.txt"))


#ML_16
top_features_16 <-read.csv( file.path(current.importance.dir, "top_features_16.csv"))$x


# ================================================================================== #
# 3. GSVA USING the 16 GENES ORDERED BY PERMUTATION IMPORTANCE ====================
#    TO DECIDE ON SIZE OF GENE SET 
# ================================================================================== #
# GSVA with cross-validation using 16-gene signature as ordered by above 5-fold CV with RF model 
gene_sets <- list()

# This is T0_upregulated reordered by importance
top_features_cv <- top_features_16 #ordered by importance in 16 genes #ML_16

print(current.roc.dir) 
if(!exists(file.path(current.roc.dir, "GSVA"))) dir.create(file.path(current.roc.dir, "GSVA"))

current.gsva.dir <- file.path(current.roc.dir, "GSVA", "ML_16")

if(!exists(current.gsva.dir)) dir.create(file.path(current.gsva.dir))

if(!exists(file.path(current.gsva.dir, "confusion_matrices"))) dir.create(file.path(current.gsva.dir, "confusion_matrices"))
if(!exists(file.path(current.gsva.dir, "roc_obj"))) dir.create(file.path(current.gsva.dir, "roc_obj"))

#save note for myself
note <- paste0(
  "5 Repeat 5 CV RF Model trained on 16 upregulated genes, assessed permutation importance on 16 genes.\n",
  str_c(top_features_cv, collapse = " "),
  "\nNow run GSVA with 16 sets")

write(note, file = file.path(current.gsva.dir,"description_of_this_run.txt"))

# Generate gene sets
for (i in seq_along(top_features_cv)) {
  gene_sets[[paste0("set_", i)]] <- top_features_cv[1:i]
  
}

# Example: Inspect the first few sets
print(gene_sets[3])

# Parameters
num_folds <- 5  # Number of folds for k-fold CV
repeats <- 5    # Number of repetitions for CV
nested_cv_results <- list()  # Store results for each gene set

library(plyr)
library("doParallel"); library("foreach")

# Iterate over all gene sets
for (set_name in names(gene_sets)) {
  start_time <- Sys.time()
  cat("Processing gene set:", set_name, "\n")
  
  # Current selected features for this gene set
  selected_features <- gene_sets[[set_name]]
  
  # The selected features form a single "gene set" for GSVA, GSVA needs list input
  gene_set_list <- list(current_set = selected_features)
  
  # Initialize storage for AUCs and confusion matrices
  repeat_aucs <- numeric()
  sensitivity <- numeric()
  specificity <- numeric()
  confusion_matrices <- list() #A confusion matrix represents the prediction summary in matrix form. It shows how many prediction are correct and incorrect per class. It helps in understanding the classes that are being confused by model as other class.
  roc_obj_list <- list()
  
  # Perform repeated k-fold CV
  #Open repeats loop
  results <- foreach(i = 1:repeats, .packages = c("caret", "GSVA", "pROC", "stats", "base")) %dopar% {
    # Set a unique seed for each repeat derived from the base seed
    set.seed(base_seed + i)
    
    # Create stratified folds for this repetition (based on disease factor)
    cv_folds <- caret::createFolds(expression_auc$Disease, k = num_folds, list = TRUE)
    
    # Manually validate/verify that in each fold, there is an even number of HC and TB
    # table(expression_auc[cv_folds[[5]], "Disease"])
    
    fold_results <- list()
    
    #Open fold loop
    for (fold in seq_along(cv_folds)) {
      # Split data into training and test sets
      test_indices <- cv_folds[[fold]]
      train_indices <- setdiff(seq_len(nrow(expression_auc)), test_indices)
      
      train_data <- expression_auc[train_indices, ]
      test_data <- expression_auc[test_indices, ]
      
      
      # Step 1: GSVA on training data (transform gene features into a composite score)
      gsvapar_train <- gsvaParam(t(train_data[, !colnames(train_data) == "Disease"]), 
                                 gene_set_list, 
                                 maxDiff = TRUE, 
                                 minSize = 1)
      gsva_train <- gsva(gsvapar_train)
      
      # Step 2: GSVA on test data
      gsvapar_test <- gsvaParam(t(test_data[, !colnames(test_data) == "Disease"]), 
                                gene_set_list, 
                                maxDiff = TRUE, 
                                minSize = 1)
      gsva_test <- gsva(gsvapar_test)
      
      # Prepare data for GLM
      glm_train_data <- data.frame(Score = gsva_train[1, ], Group = train_data$Disease)
      glm_test_data  <- data.frame(Score = gsva_test[1, ], Group = test_data$Disease)
      
      # THE GSVA SCORE IS WHAT IS USED FOR CLASSIFICATION. it is THE PREDICTOR. 
      # Since we have X amount of genes. the score should tell us if X amount is good for predicting
      
      # Step 3: Train a logistic GLM model
      # No hyperparameters tuned here; GLM is straightforward
      # Train a logistic GLM model (logistic regression is a type of GLM used to predict BINARY outcomes, whereas linear regression is used for CONTINUOUS outcomes)
      # We are testing if our predictor (GSVA score/Score) has an effect on the outcome (Disease/Group) - Since we have a gsva for for each amount of genes, this should tell us whether X amount of genes is a good enough predictor of disease
      glm_model <- glm(Group ~ Score, data = glm_train_data, family = binomial) #choose binomial because we have a binary outcome, HC or TB
      
      # Step 4: Predict probabilities on test set, based on the fitted logistic regression model above
      #For logistic regression, predict() gives us probabilities (between 0 and 1).
      #new_data is the data for which you want to make predictions, in this case it is the entire dataset
      #type = "response" for a logistic model gives you the predicted probabilities of the positive class (usually coded as 1) ie. probability of having disease
      test_probs <- predict(glm_model, newdata = glm_test_data, type = "response")
      
      
      # Step 5: Calculate ROC and determine optimal threshold
      # The ROC curve is a plot that shows the performance of a classification model across all possible threshold values. The two main metrics visualized on the ROC curve are:
      #True Positive Rate (TPR) or Sensitivity:y axis
      #False Positive Rate (FPR) or 1-Specificity: x axis
      #If the ROC curve is closer to the top-left corner (where FPR is 0 and TPR is 1), it indicates that the model is very good at distinguishing between positive and negative classes.
      roc_obj <- roc(glm_test_data$Group, test_probs) #if we access this object, it has sensitivity and specificity scores at every threshold
      

      optimal_threshold_coords <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
      optimal_threshold <- optimal_threshold_coords$threshold
      #x is the measure at which to extract coordinates = "best" the best argument refers to a method for selecting the "best" threshold on the ROC curve. The "best" threshold is often selected based on a trade-off between sensitivity (True Positive Rate) and specificity (True Negative Rate). The Youden’s Index: The threshold that maximizes sensitivity + specificity - 1.
      #  ret specifies the type of values to return, "threshold" is the cutoff threshold value.
      #threshold is the decision boundary. Above the threshold = that observation is positive, below = that observation is negative
      #if you set the threshold lower (eg 0.01) The model is more likely to predict "positive" for most instances but that = higher false positive rate and higher true positive rate. Higher threshold = lower false positive rate but also lower true positive rate 
      
      # Step 6: Convert probabilities to predictions using optimal threshold
      predicted_classes <- factor(ifelse(test_probs >= optimal_threshold, "TB", "HC"), #We assign TB to samples with a test probability more than the optimal threshold
                                  levels = c("HC", "TB"))
      actual_classes <- factor(glm_test_data$Group, levels = c("HC", "TB"))
      
      # Calculate confusion matrix for this model
      #used to evaluate the performance of a classification algorithm. It compares the actual labels (true outcomes) to the predicted labels, providing a detailed breakdown of the model's accuracy, precision, recall, F1-score, and other relevant metrics.
      confusion_matrix <- confusionMatrix(predicted_classes, actual_classes) #this contains sensitivity and specificity scoreat the chosen threshold
      
      # Store fold-level results
      fold_results[[fold]] <- list(
        Confusion_Matrix = confusion_matrix,
        AUC = auc(roc_obj), #auc calculates area under the roc_obj. 
        Fold = fold,
        Optimal_Threshold = optimal_threshold,
        ROC_obj = roc_obj,
        predicted_classes = predicted_classes,
        actual_classes = actual_classes
        
      )
      gc()
    } #close fold
    # Return the following object from the loop and store into 'results'
    fold_results
  } #close repeat
  #The AUC is the area under this curve, a larger area (closer to 1) means the model is better at distinguishing between the classes.
  #A straight line on a ROC curve = random classifier, the predictor model is no better than guessing
  
  
  # Consolidate results from all repeats and folds
  for (i in seq_along(results)) {
    for (fold in seq_along(results[[i]])) {
      confusion_matrices[[paste0("Repeat_", i, "_Fold_", fold)]] <- results[[i]][[fold]] # we made an empty list at the beginning so append to it
      repeat_aucs <- c(repeat_aucs, results[[i]][[fold]]$AUC) #we made an empty numeric at the beginning so append to it
      sensitivity <- c(sensitivity, results[[i]][[fold]]$Confusion_Matrix$byClass["Sensitivity"])
      specificity <- c(specificity, results[[i]][[fold]]$Confusion_Matrix$byClass["Specificity"])
    }
  }
  
  # Save confusion matrices for this gene set
  saveRDS(confusion_matrices, file = file.path(current.gsva.dir, "confusion_matrices", paste0("confusion_matrices_", set_name, ".RDS")))
  
  
  # Store average AUC for this gene set
  nested_cv_results[[set_name]] <- list(
    Mean_AUC = mean(repeat_aucs),
    AUCs = repeat_aucs,
    Sensitivity = sensitivity , #this is the sensitivity at the chosen threshsold 
    Mean_sensitivity = mean(sensitivity),
    Specificity = specificity, #this is the specificity at the chosen threshold 
    Mean_specificity = mean(specificity)
    
  )
  
  # Record end time and calculate elapsed time
  end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours <- floor((elapsed_seconds %% (24 * 3600)) / 3600)
  minutes <- floor((elapsed_seconds %% 3600) / 60)
  seconds <- round(elapsed_seconds %% 60)
  formatted_time <- sprintf("%02d:%02d:%02d", hours, minutes, seconds)
  
  cat(sprintf("Gene set: %s took %s (hh:mm:ss)\n", set_name, formatted_time))
} #close last gene set_16. loop done

# Save overall results
#nested_cv_results is a list containing 16 lists (set_1 to set_16). In each of these sets, there are 25x sensitivity, specificity and auc scores (because we split the data into 5 folds, used a fold for training and the rest for testing, then re-split the data 4 more times)
saveRDS(nested_cv_results, file = file.path(current.gsva.dir, "02_nested_cv_results_gene_sets.RDS"))
cat("Cross-validation and confusion matrix calculation complete.\n")

# Summarize results
summary_results <- data.frame(
  Gene_Set = names(nested_cv_results),
  Mean_AUC = sapply(nested_cv_results, function(x) x$Mean_AUC),
  t(sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3)))))
# Fold_AUCs = sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3), collapse = ", "))
# )

print(summary_results)

# summary_results <- cbind(summary_results[seq(1,nrow(summary_results),2),]) #remove single genes

plot_results <- summary_results[,-c(2)]
plot_results_long <- as.data.frame(pivot_longer(plot_results, cols= c(2:ncol(plot_results)), 
                                                names_to = "iteration", 
                                                values_to = "AUC"))
plot_results_long$num_genes <- as.numeric(gsub("set_", "", plot_results_long$Gene_Set))
plot_results_long$num_genes <- as.factor(plot_results_long$num_genes)


# ================================================================================== #
## 3.1. AUC CURVES FIGURES ==============================================================
# ================================================================================== #


auc_curves.dir <- file.path(current.gsva.dir, "figures")
if(!exists(auc_curves.dir)) dir.create(auc_curves.dir)


png(filename = file.path(auc_curves.dir, "auc_geneset_size_boxplot.png"), width = 530, height = 390, "px")
ggplot(data = plot_results_long,
       aes( x = as.factor(num_genes),
            y = as.numeric(AUC))) +
  geom_boxplot() +
  geom_point(alpha = .2)+
  scale_y_continuous( limits=c(0.5, 1)) +
  xlab("Genes") +  ylab("AUC") +
  ggtitle ("AUC Curve (5 x 5 CV)")

dev.off()

write.csv(plot_results_long, file.path(auc_curves.dir, "auc_geneset_size_boxplot_data.csv"))
plot_results_long$AUC

# ================================================================================== #
# 4. 3 DIFFERENT ORDERING METHODS  =================================================
# ================================================================================== #

# We have the 16 genes orderEd by permutation importance (top_features_16)
# Alternatively, we also have 16 genes ordered by FDR from the differential expression (T0_upregulated)
# We don't know which method of ordering is better, so will try
# 1) Ordered by Permutation importance
# 2) Ordered by FDR
# 3) Ordered by Average of the ranks from both 
intersect(top_features_16, T0_upregulated) #same 16
cbind(top_features_16, T0_upregulated)

genes_16 <- as.data.frame(rbind(cbind(top_features_16, seq(1:16)), cbind(T0_upregulated, seq(1:16))))

detach(package:plyr)

colnames(genes_16) <- c("Gene", "Rank")
genes_16$Gene <- as.factor(genes_16$Gene)
genes_16$Rank <- as.numeric(genes_16$Rank)
genes_16_mean <- group_by(genes_16, Gene) %>% summarise(mean = mean(Rank))
genes_16_mean <- as.data.frame(genes_16_mean[order(genes_16_mean$mean),])
genes_16_mean <- as.character(genes_16_mean$Gene)




# ================================================================================== #
# 5. GSVA WITH VARIOUS SIGNATURE SETS =========================================
# ================================================================================== #

# GSVA 5X5 CV using 6/7/8 gene signatureS ------------------------------------------------------------#

library(plyr)
library(dplyr)
library("doParallel"); library("foreach")

signature_set.dir <- file.path(current.roc.dir, "GSVA", "different_gene_sets")
if(!exists(signature_set.dir)) dir.create(signature_set.dir)
# GSVA gene sets

# Tess' choices
genesig_A_6 <- c("IFITM1", "P2RY14", "TAP1", "CD274", "GBP1", "GBP2")
genesig_B_6 <- c("IFITM1", "P2RY14", "S100A8", "CD274", "GBP1", "GBP2")

# Before sex correction
genesig_C_6 <- c("IFITM1", "CD274", "TAP1", "GBP5", "S100A8", "FCGR1CP") 

# With sex correction
# Top genes when ordered by significance of 16-gene signature
genesig_D_6 <- c(top_features_16[1:6]) 
genesig_D_7 <- c(top_features_16[1:7])
genesig_D_8 <- c(top_features_16[1:8]) # "IFITM1"  "CD274"   "TAP1"    "GBP5"    "GBP2"    "S100A8"  "FCGR1CP" "IFITM3" 
# genesig_D_4 = c("TAP1", "GBP5", "GBP2", "FCGR1CP") #chosen 10/03

# Top genes when ordered by Pval
genesig_E_6 <- c(T0_upregulated[1:6])
genesig_E_7 <- c(T0_upregulated[1:7])
genesig_E_8 <- c(T0_upregulated[1:8]) #"IFITM1"  "S100A8"  "FCGR1CP" "GBP2"    "IFITM3"  "P2RY14"  "GBP1"    "TAP1"

# Top genes overlapping between the two
genesig_F_6 <- genes_16_mean[1:6]
genesig_F_7 <- genes_16_mean[1:7]
genesig_F_8 <- genes_16_mean[1:8]

library(tidyverse)

#Define gene_sets
gene_sets <- list(genesig_A_6, genesig_B_6, 
                  genesig_D_6, genesig_D_7, genesig_D_8,
                  genesig_E_6, genesig_E_7, genesig_E_8,
                  genesig_F_6, genesig_F_7, genesig_F_8)


names(gene_sets) <- c("genesig_A_6", "genesig_B_6", "genesig_D_6", "genesig_D_7", "genesig_D_8",
                      "genesig_E_6", "genesig_E_7", "genesig_E_8", "genesig_F_6", "genesig_F_7", "genesig_F_8")



roc_func <- function(type, case, control){
  
  if(type == "Disease"){
    
    ### For DISEASE ----
    #Make Disease a factor
    Disease <- clinical2$Disease
    expression_auc <- as.data.frame(cbind(t(expression2), Disease = Disease))#rows are samples instead of genes here
    
    expression_auc[,"Disease"] <- factor(expression_auc[,"Disease"])
    
    #Make gene expression numeric
    expression_auc[,1:(ncol(expression_auc)-1)] <- sapply(expression_auc[,1:(ncol(expression_auc)-1)], as.numeric)
    
    #Subset
    row.names(expression_auc) == row.names(clinical2)
    T0_samples <- row.names(clinical2)[which(clinical2$condition == "HC_T0" | clinical2$condition == "TB_T0")]
    
    
    expression_auc <- expression_auc[T0_samples,]
    
    case = "TB"
    control = "HC"
  }
  
  
  if(type == "Timepoint"){
    
    ### For TIMEPOINT -----
    Timepoint <- clinical2$Timepoint
    expression_auc <- as.data.frame(cbind(t(expression2), Timepoint = Timepoint))#rows are samples instead of genes here
    
    
    #Make gene expression numeric
    expression_auc[,1:(ncol(expression_auc)-1)] <- sapply(expression_auc[,1:(ncol(expression_auc)-1)], as.numeric)
    
    row.names(expression_auc) == row.names(clinical2)
    
    
    #Subset
    subset_samples<- row.names(clinical2)[which(clinical2$condition == case | clinical2$condition == control)]
    expression_auc <- expression_auc[subset_samples,]
    
    #Make Timepoint a factor
    expression_auc[,"Timepoint"] <- factor(expression_auc[,"Timepoint"])
  }
  
  # Define parameters for 5x5 CV
  num_folds <- 5  # Number of folds for k-fold CV
  repeats <- 5    # Number of repetitions for CV
  
  
  nested_cv_results <- list()  # Store results for each gene set
  all_results <- list()
  
  base_seed = 7
  # Iterate over all gene sets
  for (set_name in names(gene_sets)) {
    
    start_time <- Sys.time()
    cat("Processing gene set:", set_name, "\n")
    
    # Current selected features for this gene set
    selected_features <- gene_sets[[set_name]]
    
    # The selected features form a single "gene set" for GSVA, GSVA needs list input
    gene_set_list <- list(current_set = selected_features)
    
    # Initialize storage for AUCs and confusion matrices
    repeat_aucs <- numeric()
    sensitivity <- numeric()
    specificity <- numeric()
    confusion_matrices <- list() #A confusion matrix represents the prediction summary in matrix form. It shows how many prediction are correct and incorrect per class. It helps in understanding the classes that are being confused by model as other class.
    
    
    # Perform repeated k-fold CV
    #foreach differs from a for loop in that its return is a list of values, whereas a for loop has no value and uses side effects to convey its result.
    results <- foreach(i = 1:repeats, .packages = c("caret", "GSVA", "pROC", "stats", "base")) %dopar% {
      # Set a unique seed for each repeat derived from the base seed
      set.seed(base_seed + i)
      
      # Create stratified folds for this repetition (based on  factor)
      cv_folds <- caret::createFolds(expression_auc[,type], k = num_folds, list = TRUE) #note, using expression instead of expression_auc$Timepoint here would have done the same thing. Because create folds just gives you indices for each fold
      
      # Manually validate/verify that in each fold, there is an even number of HC and TB
      # table(expression_auc[cv_folds[[5]], "Timepoint"])
      
      fold_results <- list()
      
      for (fold in seq_along(cv_folds)) {
        # Split data into training and test sets
        test_indices <- cv_folds[[fold]]
        train_indices <- setdiff(seq_len(nrow(expression_auc)), test_indices)
        
        train_data <- expression_auc[train_indices, ]
        test_data <- expression_auc[test_indices, ]
        
        
        # Step 1: GSVA on training data (transform gene features into a composite score)
        gsvapar_train <- gsvaParam(t(train_data[, !colnames(train_data) == type]), 
                                   gene_set_list, 
                                   maxDiff = TRUE, 
                                   minSize = 1)
        gsva_train <- gsva(gsvapar_train)
        
        # Step 2: GSVA on test data
        gsvapar_test <- gsvaParam(t(test_data[, !colnames(test_data) == type]), 
                                  gene_set_list, 
                                  maxDiff = TRUE, 
                                  minSize = 1)
        gsva_test <- gsva(gsvapar_test)
        
        # Prepare training data for GLM
        glm_train_data <- data.frame(Score = gsva_train[1, ], Group = train_data[,type])
        
        # The GLM model we make will be used to predict Timepoint based on the gsva scores in the test set
        glm_test_data  <- data.frame(Score = gsva_test[1, ], Group = test_data[,type]) 
        
        # THE GSVA SCORE IS WHAT IS USED FOR CLASSIFICATION. it is THE PREDICTOR. The gsva score is used by the model to predict Timepoint. We then apply this model to the test set to predict Timepoint. We calculate the ROC (false pos vs false neg rate) to see what threshold is a good classifier
        
        # Step 3: Train a logistic GLM model (with train)
        # No hyperparameters tuned here; GLM is straightforward
        # Train a logistic GLM model (logistic regression is a type of GLM used to predict BINARY outcomes, whereas linear regression is used for CONTINUOUS outcomes)
        # We are testing if our predictor (GSVA score/Score) has an effect on the outcome (Timepoint/Group) - Since we have a gsva for for each amount of genes, this should tell us whether X amount of genes is a good enough predictor of Timepoint
        glm_model <- glm(Group ~ Score, data = glm_train_data, family = binomial) #choose binomial because we have a binary outcome, HC or TB
        
        # Step 4: Predict probabilities on test set, based on the fitted logistic regression model above
        #For logistic regression, predict() gives us probabilities (between 0 and 1).
        #new_data is the data for which you want to make predictions, in this case it is the entire dataset
        #type = "response" for a logistic model gives you the predicted probabilities of the positive class (usually coded as 1) ie. probability of having Timepoint
        test_probs <- predict(glm_model, newdata = glm_test_data, type = "response") #this gives us numerical values for each sample, but then we have to decide on our threshold so that value can be assigned either HC or TB
        
        
        # Step 5: Calculate ROC and determine optimal threshold
        # The ROC curve is a plot that shows the performance of a classification model across all possible threshold values. The two main metrics visualized on the ROC curve are:
        #True Positive Rate (TPR) or Sensitivity:y axis
        #False Positive Rate (FPR) or 1-Specificity: x axis
        #If the ROC curve is closer to the top-left corner (where FPR is 0 and TPR is 1), it indicates that the model is very good at distinguishing between positive and negative classes.
        roc_obj <- roc(glm_test_data$Group, test_probs)
        
        # roc_obj_results[[fold]] <- roc_obj 
        
        # There are 29 different sensitivty and specificity scores, associated with 29 thresholds
        # There are 29 samples in this fold of test_prob, hence why 29 thresholds
        
        optimal_threshold_coords <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
        
        
        optimal_threshold <- optimal_threshold_coords$threshold
        
        # coords_at_optimal_threshold <- coords(roc_obj, x = optimal_threshold, ret = c("sensitivity", "specificity"))
        
        
        #x is the measure at which to extract coordinates = "best" the best argument refers to a method for selecting the "best" threshold on the ROC curve. The "best" threshold is often selected based on a trade-off between sensitivity (True Positive Rate) and specificity (True Negative Rate). The Youden’s Index: The threshold that maximizes sensitivity + specificity - 1.
        #  ret specifies the type of values to return, "threshold" is the cutoff threshold value.
        #threshold is the decision boundary. Above the threshold = that observation is positive, below = that observation is negative
        #if you set the threshold lower (eg 0.01) The model is more likely to predict "positive" for most instances but that = higher false positive rate and higher true positive rate. Higher threshold = lower false positive rate but also lower true positive rate 
        
        # Step 6: Convert probabilities to predictions/classifications using optimal threshold
        predicted_classes <- factor(ifelse(test_probs >= optimal_threshold, case, control), #We assign TB to samples with a test probability more than the optimal threshold
                                    levels = c(control, case))
        actual_classes <- factor(glm_test_data$Group, levels = c(control, case))
        
        # Calculate confusion matrix for this model
        #used to evaluate the performance of a classification algorithm. It compares the actual labels (true outcomes) to the predicted labels, providing a detailed breakdown of the model's accuracy, precision, recall, F1-score, and other relevant metrics.
        confusion_matrix <- confusionMatrix(predicted_classes, actual_classes) #this has the senstivity and specificity score at the optimal threshold, for repeat i fold i
        
        # Store fold-level results (for this repeat and fold)
        fold_results[[fold]] <- list(
          Confusion_Matrix = confusion_matrix,
          AUC = auc(roc_obj), #auc calculates area under the roc_obj. 
          Fold = fold,
          Optimal_Threshold = optimal_threshold,
          ROC_obj = roc_obj,
          predicted_classes = predicted_classes,
          actual_classes = actual_classes
          
        )
        
        gc()
      } #close fold
      # Return the following object from the loop and store into 'results'
      fold_results 
    } # close repeat
    #The AUC is the area under this curve, a larger area (closer to 1) means the model is better at distinguishing between the classes.
    #A straight line on a ROC curve = random classifier, the predictor model is no better than guessing
    
    
    # Consolidate results from all repeats and folds
    for (i in seq_along(results)) {
      for (fold in seq_along(results[[i]])) { #for each fold in this repeat
        confusion_matrices[[paste0("Repeat_", i, "_Fold_", fold)]] <- results[[i]][[fold]] #
        repeat_aucs <- c(repeat_aucs, results[[i]][[fold]]$AUC)
        sensitivity <- c(sensitivity, results[[i]][[fold]]$Confusion_Matrix$byClass["Sensitivity"])
        specificity <- c(specificity, results[[i]][[fold]]$Confusion_Matrix$byClass["Specificity"])
        
      }
    }
    
    
    # Store average AUC for this gene set
    nested_cv_results[[set_name]] <- list(
      Mean_AUC = mean(repeat_aucs),
      AUCs = repeat_aucs,
      Sensitivity = sensitivity , #this is the sensitivity at the chosen threshsold
      Specificity = specificity #this is the specificity at the chosen threshold
    )
    
    
    # Consolidate results from all repeats and folds
    for (i in seq_along(results)) {
      for (fold in seq_along(results[[i]])) { #for each fold in this repeat
        all_results[[paste0("Repeat_", i, "_Fold_", fold)]] <- results[[i]][[fold]] #
      }
    }
    
    #at the T0 back to the name for Disease 
    if(type == "Disease"){
      case = "TB_T0"
      control = "HC_T0"}
    
    if(!exists(file.path(signature_set.dir, set_name)))dir.create(file.path(signature_set.dir, set_name))
    saveRDS(all_results,file = file.path(signature_set.dir, set_name, paste0(set_name, outcome,"_",case,"vs",control,"_all_results.RDS")))
    
    # 'results' contains the specified results above (confusion matrix, AUC, senstivity and specificity) from every repeat and fold now
    
    # Record end time and calculate elapsed time
    end_time <- Sys.time()
    elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
    hours <- floor((elapsed_seconds %% (24 * 3600)) / 3600)
    minutes <- floor((elapsed_seconds %% 3600) / 60)
    seconds <- round(elapsed_seconds %% 60)
    formatted_time <- sprintf("%02d:%02d:%02d", hours, minutes, seconds)
    
    cat(sprintf("Gene set: %s took %s (hh:mm:ss)\n", set_name, formatted_time))
  }
  
  
  # Save overall results
  saveRDS(nested_cv_results, file = file.path(signature_set.dir, paste0(outcome, "_",case,"vs",control,"_nested_cv_results_different_gene_sets.RDS")))
  
  cat("Cross-validation and confusion matrix calculation complete.\n")
  
  # Summarize results
  summary_results  <- data.frame(
    Gene_Set = names(nested_cv_results),
    Mean_AUC = sapply(nested_cv_results, function(x) x$Mean_AUC),
    t(sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3)))))
  # Fold_AUCs = sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3), collapse = ", "))
  # )
  
  print(summary_results)
  
  
}






#Run the roc function with either all patients or only patients that were cured
for (outcome in c("Cured", "Any_tx_outcome")){
  
  
  #outcome loop around function calls at the bottom -------#
  #Subset for those who were cured
  if(outcome == "Cured"){
    
    clinical_cured <- clinical[which(clinical$treatment_outcome == "Cured" | clinical$Disease == "HC"),]
    expression_cured <- expression[,row.names(clinical_cured)]
    expression2 <-  expression_cured
    clinical2 <- clinical_cured
  }
  
  #Keep all patients regardless of treatment outcome
  if(outcome == "Any_tx_outcome"){
    expression2 <-  expression
    clinical2 <- clinical
  }
  
  
  roc_func(type = "Disease", case = "TB_T0", control = "HC_T0")
  roc_func(type = "Timepoint", case = "TB_T2", control = "TB_T0")
  roc_func(type = "Timepoint", case = "TB_T4", control = "TB_T0")
  roc_func(type = "Timepoint", case = "TB_T6", control = "TB_T0")
  roc_func(type = "Timepoint", case = "HC_T6", control = "HC_T0")
  
} #close outcome loop


for (outcome in c("Cured", "Any_tx_outcome")){
  #make empty table to save summarised results
  
  
  
  listofcasecontrol <- list(c("TB_T0", "HC_T0"),
                            c("TB_T2", "TB_T0"),
                            c("TB_T4", "TB_T0"),
                            c("TB_T6", "TB_T0"),
                            c("HC_T6", "HC_T0"))
  
  combined_summary_results <- data.frame(Gene_set = names(gene_sets),
                                         Genes = unlist(lapply(gene_sets, function (x) 
                                           paste(x,collapse = ",")
                                         )))
  
  for (i in 1:length(listofcasecontrol)){
    nested_cv_results <- readRDS(file.path(signature_set.dir, paste0(outcome, "_",listofcasecontrol[[i]][1],"vs",listofcasecontrol[[i]][2],"_nested_cv_results_different_gene_sets.RDS")))
    
    # Summarize results
    summary_results  <- data.frame(
      Gene_Set = names(nested_cv_results),
      Mean_AUC = sapply(nested_cv_results, function(x) x$Mean_AUC),
      t(sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3)))))
    # Fold_AUCs = sapply(nested_cv_results, function(x) paste(round(x$AUCs, 3), collapse = ", "))
    # )
    
    combined_summary_results <- cbind(combined_summary_results,  summary_results$Mean_AUC)
    colnames(combined_summary_results)[ncol(combined_summary_results)] <- paste0(listofcasecontrol[[i]][1],"vs",listofcasecontrol[[i]][2])
    
    
    #save the table
    write.csv(combined_summary_results, file = file.path(signature_set.dir, paste0(outcome, "_all_different_gene_sets_auc_results.csv")))
    
  }
}



# ================================================================================== #
## 5.1 BOXPLOTS COMPARING THE PERFORMANCES OF VARIOUS SIGNATURES  ==================
# ================================================================================== #
# Seperate boxplots for each timepoint/disease comparisons
# X axis is all the different signatures and y axis is the AUCs

boxplot_theme <- theme(axis.title = element_text(size = 24),
                       axis.text = element_text(size = 24),
                       title = element_text(size = 20),
                       legend.position = "None") 

contrasts <- c(
  "TB_T0vsHC_T0",
  "TB_T2vsTB_T0",
  "TB_T4vsTB_T0",
  "TB_T6vsTB_T0",
  "HC_T6vsHC_T0")

for(outcome in c("Cured", "Any_tx_outcome")){
  for (contrast in contrasts){
    nested_res <- readRDS(file.path(signature_set.dir, 
                                    paste0(outcome, "_", contrast, "_nested_cv_results_different_gene_sets.RDS"))) 
    
    nested_res_long <- do.call(rbind, lapply(names(nested_res), function(signature) {
      
      auc_values <- nested_res[[signature]]$AUCs
      data.frame(
        name =  rep(signature, length(auc_values)),
        aucs = auc_values)
    } ))
    
    
    
    boxplot <- ggplot(nested_res_long, aes(
      x = factor(name),
      y = as.numeric(aucs),
      colour = name)) +
      
      theme_bw()+
      
      boxplot_theme +
      
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      geom_boxplot(position = position_dodge(1)) +
      
      geom_jitter(aes(color = name),
                  alpha = 0.5,
                  size = 2,
                  width = 0.3) +
      stat_summary(fun.y = mean, fill = "red",
                   geom = "point", shape = 21, size =4,
                   show.legend = TRUE) +
      # # scale_x_discrete(labels= c("Control" = "Control", "Mild.moderate.COPD" = "mCOPD", "Severe.COPD" = "sCOPD"))+
      # scale_y_continuous(expand = c(0.07, 0, 0.07, 0)) +
      
      labs(title = paste("TB Signature Performance:", contrast ))+
      ylab (label = "AUC Score") +
      xlab (label = "Signature")
    
    #save boxplot
    if(!exists(file.path(signature_set.dir, "boxplot"))) dir.create(file.path(signature_set.dir, "boxplot"))
    ggsave(boxplot, filename = file.path(signature_set.dir, "boxplot", paste0(outcome, "_different_gene_sets_" ,contrast, "_boxplot.png")),
           width = 3500,
           height = 3200,
           units = "px" )
    
    #save data for the boxplot
    if(!exists(file.path(signature_set.dir, "boxplot_data"))) dir.create(file.path(signature_set.dir, "boxplot_data"))
    write.csv(nested_res_long, file.path(signature_set.dir, "boxplot_data", paste0(outcome, "_different_gene_sets_" ,contrast, "_boxplot_data.csv")))
    
    
  } #close contrast loop
}# close outcome loop





# ================================================================================== #
# GENESIG_D_7 PERFORMANCE AS A SIGNATURE ============================================
# ================================================================================== #
genesig_D_7_performance.dir <- file.path(signature_set.dir, "genesig_D_7")
if(!exists(genesig_D_7_performance.dir)) dir.create(genesig_D_7_performance.dir)

genesig_D_7_performance.figures.dir <- file.path(genesig_D_7_performance.dir, "figures")
if(!exists(genesig_D_7_performance.figures.dir)) dir.create(genesig_D_7_performance.figures.dir)





# ================================================================================== #
# 6. PLOT MIN/MAX ROC CURVES FOR GENESIG_D_7 =======================================
# ================================================================================== #
#each point of the curve represents a different threshold  (high threshold = most sample will be negative = high true negative rate (high specificity) and low true positibe rate (low sensitivity) because eveything is getting classified as a negative)

for (outcome in c("Cured", "Any_tx_outcome")){
  genesig_D_7_TBT0vsHCT0 <- readRDS(file.path(genesig_D_7_performance.dir, paste0("genesig_D_7", outcome, "_TB_T0vsHC_T0", "_all_results.RDS"))) ##change to csv
  genesig_D_7_TBT2vsTBT0 <- readRDS(file.path(genesig_D_7_performance.dir, paste0("genesig_D_7", outcome, "_TB_T2vsTB_T0", "_all_results.RDS")))
  genesig_D_7_TBT4vsTBT0 <- readRDS(file.path(genesig_D_7_performance.dir, paste0("genesig_D_7", outcome, "_TB_T4vsTB_T0", "_all_results.RDS")))
  genesig_D_7_TBT6vsTBT0 <- readRDS(file.path(genesig_D_7_performance.dir, paste0("genesig_D_7", outcome, "_TB_T6vsTB_T0", "_all_results.RDS")))
  genesig_D_7_HCT6vsHCT0 <- readRDS(file.path(genesig_D_7_performance.dir, paste0("genesig_D_7", outcome, "_HC_T6vsHC_T0", "_all_results.RDS")))
  
  listof_D7_results <- list(TBT0vsHCT0 = genesig_D_7_TBT0vsHCT0, 
                            TBT2vsTBT0 = genesig_D_7_TBT2vsTBT0, 
                            TBT4vsTBT0 = genesig_D_7_TBT4vsTBT0, 
                            TBT6vsTBT0 = genesig_D_7_TBT6vsTBT0, 
                            HCT6vsHCT0 = genesig_D_7_HCT6vsHCT0)
  
  
  
  
  #listofresults only containing roc_obj results
  roc_combined_list <- list()
  for (comparison in names(listof_D7_results)){
    comparison_list <- listof_D7_results[[comparison]]
    
    roc_list <- list()
    for (repeat_fold in names(comparison_list)) {
      roc_obj <- comparison_list[[repeat_fold]]$ROC_obj
      roc_list[[repeat_fold]] <- roc_obj
    }
    roc_combined_list[[comparison]] <- roc_list
  }
  
  
  
  #calculate fpr and tpr and store in roc_results_list
  roc_results_list <-list()
  
  for (comparison in names(roc_combined_list)){
    this_comparison <- roc_combined_list[[comparison]]
    
    for (repeat_fold in names(this_comparison)) {
      
      roc_obj <- this_comparison[[repeat_fold]]
      
      tpr_values <- this_comparison[[repeat_fold]]$sensitivities
      fpr_values <- 1 - this_comparison[[repeat_fold]]$specificities
      
      roc_obj$tpr <- tpr_values
      roc_obj$fpr <- fpr_values
      
      this_comparison[[repeat_fold]] <- roc_obj
    }
    
    roc_results_list[[comparison]] <- this_comparison
  }
  
  
  
  
  
  #Get min and max AUCs for every comparison
  roc_to_plot <- list()
  for (comparison in names(roc_results_list)){
    aucs <- as.numeric((lapply(roc_results_list[[comparison]], function(x) x$auc)))
    print(min(aucs))
    print(max(aucs))
    
    min_index <- which(aucs == min(aucs))[1] 
    max_index <- which(aucs == max(aucs))[1] #Some have the same max auc?????? like TBT4vsTBT0 has the same max auc at index 7 and 17? Just took the first one for now 
    
    fold_to_plot_min <- roc_results_list[[comparison]][[min_index]]
    fold_to_plot_max <- roc_results_list[[comparison]][[max_index]]
    
    roc_to_plot[[paste0(comparison, "_min")]] <- fold_to_plot_min
    roc_to_plot[[paste0(comparison, "_max")]] <- fold_to_plot_max
    
  }
  
  minmax_res <- do.call(rbind, lapply(roc_to_plot, function(x) {
    data.frame(auc = x$auc)}))
  
  minmax_res$group <- gsub("_.*", "", row.names(minmax_res))
  minmax_res$best_worst<- gsub(".*_", "", row.names(minmax_res))
  
  minmax_res <-pivot_wider(minmax_res, names_from = "best_worst", values_from = "auc")
  
  this_mean <- data.frame()
  roc_means <- data.frame()
  for (comparison in names(roc_results_list)){
    aucs <- as.numeric((lapply(roc_results_list[[comparison]], function(x) x$auc)))
    
    this_mean <- cbind(comparison,mean(aucs))
    roc_means <- rbind(roc_means,this_mean)}
  
  minmax_res <- cbind(comparison = minmax_res$group, 
                      mean_auc = roc_means$V2,
                      best_auc = minmax_res$max,
                      worst_auc = minmax_res$min)
  
  write.csv(minmax_res, file.path(genesig_D_7_performance.dir,paste0("genesig_D_7", outcome,"minmax_res_summary.csv")))
  
  
  
  ## MIN/MAX PLOTS GGPLOT2 ------------------------------------------------
  
  roc_to_plot[[paste0(comparison, "_min")]] <- fold_to_plot_min
  roc_to_plot[[paste0(comparison, "_max")]] <- fold_to_plot_max
  
  
  roc_data <- do.call(rbind, lapply(names(roc_to_plot), function(comparison) {
    data.frame(
      TPR = rev(roc_to_plot[[comparison]]$sensitivities),  # True Positive Rate
      FPR = rev(1 - roc_to_plot[[comparison]]$specificities),  # False Positive Rate
      auc = rev(roc_to_plot[[comparison]]$auc),
      comparison = comparison
      
    )
  }))
  
  
  roc_data$min_or_max <- gsub(".*_", "", roc_data$comparison)
  roc_data[which(roc_data$min_or_max == "min"), "min_or_max"] <- "Worst"
  roc_data[which(roc_data$min_or_max == "max"), "min_or_max"] <- "Best"
  
  roc_data$group <- gsub("_.*", "", roc_data$comparison)
  
  rocplot <- ggplot(roc_data, aes(x = FPR,
                                  y = TPR, 
                                  color = group,
                                  linetype = min_or_max)) +
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
      title = "Signature Performance",
      x = "FPR (1 - Specificity)",
      y = "TPR (Sensitivity)",
      caption = "Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP") 
  
  ggsave(rocplot, 
         filename = file.path(genesig_D_7_performance.figures.dir, paste0("genesig_D_7", outcome,"ROC_best_and_worst.png")), 
         width = 3000, 
         height = 3200, 
         units = "px")
  
} #close outcome loop

# ================================================================================== #
# 7. GSVA PLOT BOXPLOTS FOR GENESIG_D_7 =================================================
# ================================================================================== #


gene_set_list <- list(c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP"))



#Run the roc function with either all patients or only patients that were cured
for (outcome in c("Cured", "Any_tx_outcome")){
  
  
  #outcome loop around function calls at the bottom -------#
  #Subset for those who were cured
  if(outcome == "Cured"){
    
    clinical_cured <- clinical[which(clinical$treatment_outcome == "Cured" | clinical$Disease == "HC"),]
    expression_cured <- expression[,row.names(clinical_cured)]
    expression2 <-  expression_cured
    clinical2 <- clinical_cured
  }
  
  #Keep all patients regardless of treatment outcome
  if(outcome == "Any_tx_outcome"){
    expression2 <-  expression
    clinical2 <- clinical
  }

gsvapar <- gsvaParam(as.matrix(expression2), 
                     gene_set_list, 
                     maxDiff = TRUE, 
                     minSize = 1)

gsva_res <- gsva(gsvapar) #dont need to transpose because next line takes row 1 anyway


all(row.names(gsva_res) == row.names(clinical2))



boxplot_gsva <- as.data.frame(cbind(gsva = t(gsva_res),
                                    group = as.character(clinical2$condition),
                                    PID = as.character(clinical2$PID)))



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
# library(ggplot2)
# 
# plot(density(expression2)) #quite a long tail on the right - not normal
# ks.test(expression2, "pnorm", mean(expression2), sd(expression2))
# 
# library("nortest")
# ad.test(expression2)
# 
# ks.test(expression2, "pnorm")
# 
# qqnorm(expression2)
# qqline(expression2, col="red")
# 

# ================================================================================== #
##  PAIRED COMPARISONS  ===================================================================
# ================================================================================== #


boxplot_gsva$V1 <- as.numeric(boxplot_gsva$V1)
boxplot_gsva$group <- factor(boxplot_gsva$group)

# View(pivot_wider(boxplot_gsva, names_from = group, values_from = PID))

# We need to do per-comparison filtering
# We can't just filter out the unpaired samples and run a wilcox_test on that because a paired wilcoxin rank sum test will expect one value from each group 
# Say we are comparing TB_T0 and TB_T6. There may be patients with TB_T0 and TB_T2 but not TB_T6. They’ll still be present in the data and wilcox_test will extract the TB_T0 for that patient (since the comparison specificed TB_T0 and TB_T6. However, since the test is paired it won't work because wilcox_test needs one measurement from Group1 (TB_T0) and one from Group2 (TB_T6)
# When comparison = TB_T0 vs TB_T6, wilcoxin rank sum test will expect one value from Group 1 (TB_T0) and one value from Group 2 (TB_T6)
# Tried it paired and got the error below #'x' and 'y' must have the same length (because both groups don't have the same number of samples, since some patients don't have paired data)
# Unpaired will work but its not what we want here (it works because wilcox_test will rank all values from both groups and compare group ranks (not per patient) =

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
boxplotfinal <- ggplot(boxplot_gsva, aes(
  x = factor(group),
  y = as.numeric(boxplot_gsva[,1],
                 group = group))) +
  
  theme_bw()+
  
  gsva_theme +
  
  geom_boxplot(position = position_dodge(1)) +
  
  
  stat_pvalue_manual(stat.table.gsva.all,
                     label = "p",
                     tip.length = 0.01,
                     size = 6)+
  
  stat_summary(fun.y = mean, fill = "red",
               geom = "point", shape = 21, size =4,
               show.legend = TRUE) +
  
  labs(title = paste0("Signature Analysis"),
       caption = paste("Signature: IFITM1, CD274, TAP1, GBP5, GBP2, S100A8, FCGR1CP\n",
                       "Wilcoxin rank-sum test performed for paired samples\n",
                       "Mann-Whitney-U performed for HC vs TB comparisons\n",
                       "Enrichment Score: GSVA score for the 7 genes")) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition")

boxplotfinal2 <- boxplotfinal +  
  geom_jitter(aes(color = group),
              alpha = 0.5,
              size = 2.5,
              width = 0.3) 

ggsave(boxplotfinal2, filename = file.path(genesig_D_7_performance.figures.dir, paste0("genesig_D_7_", outcome, "_gsva_paired_boxplot.png")), 
       width = 3500, 
       height = 3200, 
       units = "px" )

boxplotfinal_lines <- boxplotfinal +     geom_point(aes(color = group))+
 geom_line(aes(group = PID), color = "black", alpha = 0.2) 
ggsave(boxplotfinal_lines, filename = file.path(genesig_D_7_performance.figures.dir, paste0("genesig_D_7_", outcome, "_gsva_paired_boxplot_lines.png")), 
       width = 3500, 
       height = 3200, 
       units = "px" )

write.csv(boxplot_gsva, file.path(genesig_D_7_performance.figures.dir, paste0("genesig_D_7_", outcome, "_boxplot_gsva_paired_boxplot_data.csv")))


}


# ================================================================================== #
# GENESIG_D_7 MORE PLOTS ============================================
# ================================================================================== #
genesig_D_7.dir <- file.path(output.dir, "2.1_genesig_D_7", "GSVA")
if(!exists(genesig_D_7.dir)) dir.create(genesig_D_7.dir, recursive = TRUE)
# ================================================================================== #
## HEATMAP ON OWN DATA  (ALL GENES) =====================================================
# ================================================================================== #
#X axis: Samples
#Y-axis: genes
#data: expression

clinical <- read.csv(file.path(main.dir,"data/processed/post-QC/clinical.csv"), row.names = 1, check.names = FALSE)
expression <- read.csv(file.path(main.dir,"data/processed/post-QC/expression.csv"), row.names = 1, check.names = FALSE)
tT_contrast1 <- read.csv(file.path(output.dir,"1_differential_expression/tT/contrast1.csv"), row.names = 1, check.names = FALSE)
top_features_16 <- read.csv(file.path(current.importance.dir,"top_features_16.csv"), row.names = 1, check.names = FALSE)
top_features_16 <- top_features_16$x
# top_features_4 <- c("TAP1", "GBP5", "GBP2", "FCGR1CP")

gene_order <- c(top_features_16, setdiff(row.names(expression), top_features_16))

# heatmapsamples <- clinical$PAXGENE[which(clinical$condition == "TB_T0" | clinical$condition == "HC_T0")] #Only T0
# heatmapsamples <- clinical$PAXGENE[which(clinical$Disease == "TB")]

heatmapsamples <- clinical$PAXGENE

clinical_heatmap <- clinical[as.character(heatmapsamples),]
expression_heatmap <- expression[,as.character(heatmapsamples)]

clinical_heatmap_ordered <- clinical_heatmap[order(clinical_heatmap$condition),] #all timepoints


#reorder expression columns/samples to be the same as clinical_ordered
heatmapdata <- as.matrix(expression_heatmap[,row.names(clinical_heatmap_ordered)])
heatmapdata <- heatmapdata[-which(row.names(heatmapdata) == "B2M" | row.names(heatmapdata) == "GAPDH"),]
# heatmapdata <-rbind(heatmapdata,condition = clinical_heatmap_ordered$condition)


sample_annot <- data.frame(Sample = clinical_heatmap_ordered$PAXGENE, 
                           condition = clinical_heatmap_ordered$condition)
sample_annot$Sample <- as.character(sample_annot$Sample)

library(ggtext)
p <- as.data.frame(heatmapdata) %>% 
  rownames_to_column("Gene") %>% 
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "expression")

p <- left_join(p, sample_annot, by = "Sample") 
p <- p %>% 
  # group_by(Gene) %>%
  mutate(expression_z = as.vector(scale(expression, scale = TRUE, center = TRUE))) #%>%
# ungroup()

p <- mutate(p, Gene_label = ifelse(p$Gene %in% top_features_16[1:7],
                                   paste0("**", p$Gene, "**"),  # Bolded label
                                   p$Gene))


p$Gene_label <- factor(p$Gene, levels = rev(gene_order))
for (i in top_features_16[1:7]){
  
  levels(p$Gene_label)[levels(p$Gene_label) == i] <-  paste0("**", i, "**")
  
}




heatmapplot <- ggplot(p, aes(x=Sample, y=Gene_label, fill=expression_z))+
  geom_tile(colour="white", linewidth=0.1)+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 14),
        axis.text.y = element_markdown(size = 10),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 14),
        plot.caption = element_text(size = 14))+
  guides(fill=guide_colorbar(title="Z-Scored \nNormalised Expression"))+
  facet_grid(. ~ condition, scales = "free_x", space = "free_x")+
  scale_fill_gradient2(
    low = "blue",       
    mid = "white",     
    high = "red",    
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  labs(title = "Gene expression at all timepoints",
       caption = "Fluidigm data normalised using Cyclic Loess Normalisation. 
       Z-Scored normalised data shown.")

print(heatmapplot)

this.figure.dir <- file.path(genesig_D_7.dir, "heatmap")
if(!exists(this.figure.dir))dir.create(this.figure.dir, recursive = TRUE)


ggsave(heatmapplot, filename = file.path(this.figure.dir, "heatmap_alltimepoints.png"), 
       width = 3200, 
       height = 2800, 
       units = "px")






# ================================================================================== #
# 8. CORRELATE CLINICAL FEATURES WITH GENESIG_D_7 ==================================
# ================================================================================== #
clinical_correlation.dir <- file.path(genesig_D_7.dir, "clinical_correlation")
if(!exists(clinical_correlation.dir)) dir.create(clinical_correlation.dir)




## Patient demographics--------
current_clinical <- clinical[-c(which(is.na(clinical$sex))),]
current_clinical[which(is.na(current_clinical$sex)), "sex"] <- "NA"
current_clinical[which(is.na(current_clinical$smoking_status)), "smoking_status"] <- "NA"

current_clinical$age_category <- cut(as.numeric(current_clinical$age),
                                     breaks=c(15, 25, 35, 45, 55, 65, Inf),
                                     labels = c("15-24", "25-34", "35-44", "45-54", "55-64", "65+"),
                                     right = FALSE) #range will include lower bound but not the upper bound

detach("package:plyr", unload = TRUE)

library("dplyr")
clinical_demographics <- t(current_clinical %>% 
                            dplyr::group_by(condition) %>% 
                            dplyr::summarise(
                              total_patients = n(),  # Count total patients in each disease group
                              male_patients = sum(sex == "male"),  # Count male patients in each disease group
                              female_patients = sum(sex == "female"),  # Count male patients in each disease group
                              no_sex_data = sum(sex == "NA"),
                              neversmoker = sum (smoking_status == "never"),
                              currentsmoker = sum (smoking_status == "current"),
                              former = sum (smoking_status == "former"),
                              no_smoke_data = sum(smoking_status == "NA"),
                              age15_24 = sum (age_category == "15-24"),
                              age25_34 = sum (age_category == "25-34"),
                              age35_44 = sum (age_category == "35-44"),
                              age45_54 = sum (age_category == "45-54"),
                              age55_64 = sum (age_category == "55-64"),
                              age65plus= sum (age_category == "65+")
                            ))

write.csv(clinical_demographics,file.path(clinical_correlation.dir, "counts_per_clinicalfeature.csv"))


# ================================================================================== #
## 8.1 GSVA TB T0 vs HC T0 =========================================================
# ================================================================================== #
if(!exists(file.path(clinical_correlation.dir, "disease"))) dir.create(file.path(clinical_correlation.dir, "disease"))
  
genesig_D_7 <- c("IFITM1","CD274","TAP1","GBP5","GBP2","S100A8","FCGR1CP")
listofgroups <- list(genesig_D_7)

T0_samples <- row.names(clinical)[which(clinical$condition == "HC_T0" | clinical$condition == "TB_T0")]
expression_T0 <- expression[,T0_samples]
clinical_T0 <- clinical[T0_samples,]

library(GSVA)
gsvapar <- gsvaParam(as.matrix(expression_T0), listofgroups, maxDiff = TRUE)
gsva_res = gsva(gsvapar)

gsva_res=t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd 


### SEX BOXPLOT -------------------------------------------------------------
boxplot=cbind(Score = gsva_res,
              Disease = clinical_T0$Disease, 
              Condition = as.character(clinical_T0$condition), 
              Age = clinical_T0$age, 
              Sex = clinical_T0$sex, 
              Smoking_status =clinical_T0$smoking_status)

boxplot <- as.data.frame(boxplot)
colnames(boxplot)[1] <- "Score"
boxplot$Score <- as.numeric(boxplot$Score)


#Remove NAs
boxplot2 <- boxplot[-which(is.na(boxplot$Sex)),]

  
library(rstatix)


#Get P-values for WITHIN disease group comparisons (male vs female within each disease)
stat.table <- boxplot2  %>%
  group_by(Disease) %>% 
  wilcox_test(Score ~ Sex,
              paired = FALSE) %>%
  add_xy_position(x = "Disease")
stat.table <- stat.table[which(stat.table$p < 0.05),]
stat.table$y.position <- max(boxplot2$Score) + 0.05*(max(boxplot2$Score))
# stat.table$y.position <- as.numeric(stat.table$y.position)
# stat.table <- stat.table[which(stat.table$p < 0.05),]

#Pairwise comparison BETWEEN disease groups
stat.table2 <-  boxplot2  %>%
  group_by(Sex) %>%
  wilcox_test(Score ~ Disease,
              paired = FALSE) %>%
  add_xy_position(x = "Disease", group = "Sex")
#Current y positions of brackets are too close to eachother. double it
new.bracket.distance <- 0.1
stat.table2 <- stat.table2[which(stat.table2$p < 0.05),]

for (i in 1:length(stat.table2$y.position)){
  stat.table2$y.position[i] <- max(stat.table$y.position) + new.bracket.distance*i
}

library(ggpubr)
#Split into two plots, one for HC and one for TB and facet? or ggarrange 
boxplotfinal2 <- ggplot(boxplot2, aes(
  x = as.factor(Disease),
  y = Score)) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5),
        legend.position = "bottom") +
  
  geom_boxplot(aes(fill = Sex)) +
  
  
  labs(title = label,
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  #within groups
  stat_pvalue_manual(stat.table,
                     label = "p",
                     tip.length = 0.01,
                     size = 4) +
  #between groups
  stat_pvalue_manual(stat.table2,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)

ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "disease", paste0("sex",score,".png")),  width = 2000,
       height = 2000,
       units = "px" )



### Smoking BOXPLOT -----------------------------------------------------------------------------------------------

boxplot_gsva=cbind(gsva = gsva_res, Disease = clinical[T0_samples,]$Disease, 
                   Condition = as.character(clinical[T0_samples,]$condition), 
                   Age = clinical[T0_samples,]$age, 
                   Sex = clinical[T0_samples,]$sex, 
                   Smoking_status = clinical[T0_samples,]$smoking_status)

boxplot_gsva <- as.data.frame(boxplot_gsva)
boxplot_gsva$genesig_D_6 <- as.numeric(boxplot_gsva$genesig_D_6)
boxplot_gsva$genesig_D_7 <- as.numeric(boxplot_gsva$genesig_D_7)

#remove NAs
boxplot_gsva <- boxplot_gsva[-which(is.na(boxplot_gsva$Smoking_status)),]

library(rstatix)

#Get P-values for WITHIN disease group comparisons (male vs female within each disease)
stat.table.gsva <- boxplot_gsva  %>%
  group_by(Disease) %>% 
  wilcox_test(genesig_D_7 ~ Smoking_status,
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "Disease")
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
stat.table.gsva$y.position <- max(boxplot_gsva$genesig_D_7) + 0.05*(max(boxplot_gsva$genesig_D_7))
new.bracket.distance <- 0.1
# stat.table.gsva$y.position[c(2)] <- max(boxplot_gsva$genesig_D_7) + new.bracket.distance

#Pairwise comparison BETWEEN disease groups
stat.table.gsva2 <-  boxplot_gsva  %>%
  group_by(Smoking_status) %>%
  wilcox_test(genesig_D_7 ~ Disease,
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "Disease", group = "Smoking_status")
stat.table.gsva2 <- stat.table.gsva2[which(stat.table.gsva2$p < 0.05),]
for (i in 1:length(stat.table.gsva2$y.position)){
  stat.table.gsva2$y.position[i] <- max(stat.table.gsva$y.position) + new.bracket.distance*i
}

library(ggpubr)
#Split into two plots, one for HC and one for TB and facet? or ggarrange 
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(Disease),
  y = as.numeric(genesig_D_7))) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5),
        legend.position = "bottom") +
  
  geom_boxplot(aes(fill = Smoking_status)) +
  
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  #within groups
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 4) +
  #between groups
  stat_pvalue_manual(stat.table.gsva2,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)

ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "disease", "smoking_gsva_genesig_D_7_sig.png"),  width = 2000,
       height = 2000,
       units = "px" )



### Age BOXPLOT ----------------------------------------------------------------------------------------------------------  

boxplot_gsva=cbind(gsva = gsva_res, Disease = clinical[T0_samples,]$Disease, 
                   Condition = as.character(clinical[T0_samples,]$condition), 
                   Age = clinical[T0_samples,]$age, 
                   Sex = clinical[T0_samples,]$sex, 
                   Smoking_status = clinical[T0_samples,]$smoking_status)

boxplot_gsva <- as.data.frame(boxplot_gsva)
boxplot_gsva$genesig_D_6 <- as.numeric(boxplot_gsva$genesig_D_6)
boxplot_gsva$genesig_D_7 <- as.numeric(boxplot_gsva$genesig_D_7)


# table(clinical$sex)
boxplot_gsva$age_category <- cut(as.numeric(boxplot_gsva$Age),
                                 breaks=c(15, 25, 35, 45, 55, 65, Inf),
                                 labels = c("15-24", "25-34", "35-44", "45-54", "55-64", "65+"),
                                 right = FALSE) #range will include lower bound but not the upper bound


#### Group by Age (X AXIS = DISEASE)  -----------------------------------------------
#Get P-values for WITHIN disease group comparisons (male vs female within each disease)
stat.table.gsva <- boxplot_gsva  %>%
  group_by(Disease) %>% 
  wilcox_test(genesig_D_7 ~ age_category,
              paired = FALSE) %>%
  add_xy_position(x = "Disease")
lowest.bracket <- max(boxplot_gsva$genesig_D_7) 
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
new.bracket.distance <- 0.15
for (i in 1:length(stat.table.gsva$y.position)){
  stat.table.gsva$y.position[i] <-lowest.bracket + new.bracket.distance*i
}

# Comment this out to hide the disease comparison p-values
#Pairwise comparison BETWEEN disease groups
stat.table.gsva2 <-  boxplot_gsva  %>%
  group_by(age_category) %>%
  wilcox_test(genesig_D_7 ~ Disease,
              paired = FALSE) %>%
  add_xy_position(x = "Disease", group = "age_category")

stat.table.gsva2 <- stat.table.gsva2[which(stat.table.gsva2$p < 0.05),]

for (i in 1:length(stat.table.gsva2$y.position)){
  stat.table.gsva2$y.position[i] <- max(stat.table.gsva$y.position) + new.bracket.distance*i
}



boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(Disease),
  y = as.numeric(genesig_D_7))
) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5)) +
  
  geom_boxplot(aes(fill = age_category)) +
  
  scale_fill_manual(values = c("15-24" = "yellow", 
                               "25-34" = "skyblue", 
                               "35-44" = "purple", 
                               "45-54" = "pink", 
                               "55-64" = "orange", 
                               "65+" = "red"))+
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  
  #within groups
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 4) +
  # #between groups # can combine all stats into one graph but too messy
  stat_pvalue_manual(stat.table.gsva2,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)


ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "disease", "age_gsva_genesig_D_7_combined.png"),  
       width = 2500,
       height = 2000,
       units = "px" )




#### Group by Disease (X AXIS = AGE) -----------------------------------

#Pairwise comparison BETWEEN disease groups
stat.table.gsva2 <-  boxplot_gsva  %>%
  group_by(age_category) %>%
  wilcox_test(genesig_D_7 ~ Disease,
              paired = FALSE) %>%
  add_xy_position(x = "age_category", group = "Disease")

lowest.bracket <- max(boxplot_gsva$genesig_D_7) + 0.01*(max(boxplot_gsva$genesig_D_7))
stat.table.gsva2 <- stat.table.gsva2[which(stat.table.gsva2$p < 0.05),]


boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(age_category),
  y = as.numeric(genesig_D_7))
) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5)) +
  
  geom_boxplot(aes(fill = Disease)) +
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Age") +
  
  
  #between groups
  stat_pvalue_manual(stat.table.gsva2,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)


ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "disease","age_gsva_genesig_D_7_comparedisease.png"),  
       width = 2500,
       height = 2000,
       units = "px" )





# ================================================================================== #
## 8.2 GSVA TB T0 vs T2, vs T4 vs T6 ===============================================
# ================================================================================== #
if(!exists(file.path(clinical_correlation.dir, "timepoint"))) dir.create(file.path(clinical_correlation.dir, "timepoint"))

listofgroups <- list(genesig_D_6, genesig_D_7)

TB_samples <- row.names(clinical)[which(clinical$Disease == "TB")]
expression_TB <- expression[,TB_samples]
clinical_TB <- clinical[TB_samples,]

#if p<0.05, then not normally distributed. only female.TB appears normal
shapiro.test(expression_TB[which(clinical_TB$sex == "male" & clinical_TB$Disease == "TB" )])
shapiro.test(expression_TB[which(clinical_TB$smoking_status == "never" & clinical_TB$Disease == "TB" )])
# replace shaprio.test with hist() or qqnorm() to visualise normality

library(GSVA)
gsvapar <- gsvaParam(as.matrix(expression_TB), listofgroups, maxDiff = TRUE)
gsva_res = gsva(gsvapar)

gsva_res=t(gsva_res) #the results tell us how much the set of genes was represented in each sample. ie. enrichment score of 0.9 is high- meaning the genes of interest showed up alot in sample X - now when we group the samples by copd and non copd, we can see whether certain genes are enriched in samples with or without copd 
colnames(gsva_res)=c("genesig_D_6", "genesig_D_7")










### Sex BOXPLOT -----------------------------------------------------------------------------------------------
boxplot_gsva=cbind(gsva = gsva_res, Disease = clinical_TB$Disease, 
                   Condition = as.character(clinical_TB$condition), 
                   Age = clinical_TB$age, 
                   Sex = clinical_TB$sex, 
                   Smoking_status =clinical_TB$smoking_status)

boxplot_gsva <- as.data.frame(boxplot_gsva)
boxplot_gsva$genesig_D_6 <- as.numeric(boxplot_gsva$genesig_D_6)
boxplot_gsva$genesig_D_7 <- as.numeric(boxplot_gsva$genesig_D_7)


#Remove NAs
boxplot_gsva <- boxplot_gsva[-which(is.na(boxplot_gsva$Sex)),]

row.names(clinical) == row.names(gsva_res)

table(clinical$sex)


library(rstatix)

#Get P-values for WITHIN disease group comparisons (male vs female within each disease)
stat.table.gsva <- boxplot_gsva  %>%
  group_by(Condition) %>% 
  wilcox_test(genesig_D_7 ~ Sex,
              paired = FALSE) %>%
  add_xy_position(x = "Condition")
stat.table.gsva$y.position <- max(boxplot_gsva$genesig_D_7) + 0.05*(max(boxplot_gsva$genesig_D_7))
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]

#Pairwise comparison BETWEEN Condition groups
stat.table.gsva2 <-  boxplot_gsva  %>%
  group_by(Sex) %>%
  wilcox_test(genesig_D_7 ~ Condition,
              paired = FALSE) %>%
  add_xy_position(x = "Condition", group = "Sex")
#Current y positions of brackets are too close to eachother. double it
new.bracket.distance <- 0.1
stat.table.gsva2 <- stat.table.gsva2[which(stat.table.gsva2$p < 0.05),]

for (i in 1:length(stat.table.gsva2$y.position)){
  stat.table.gsva2$y.position[i] <- max(stat.table.gsva$y.position) + new.bracket.distance*i
}

library(ggpubr)
#Split into two plots, one for HC and one for TB and facet? or ggarrange 
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(Condition),
  y = as.numeric(genesig_D_7))) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5),
        legend.position = "bottom") +
  
  geom_boxplot(aes(fill = Sex)) +
  
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  #within groups
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 4) +
  #between groups
  stat_pvalue_manual(stat.table.gsva2,
                     label = "p",
                     tip.length = 0.01,
                     size = 4)

ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "timepoint","sex_gsva_genesig_D_7.png"),  width = 2000,
       height = 2000,
       units = "px" )



### Smoking BOXPLOT -----------------------------------------------------------------------------------------------

boxplot_gsva=cbind(gsva = gsva_res, Disease = clinical_TB$Disease, 
                   Condition = as.character(clinical_TB$condition), 
                   Age = clinical_TB$age, 
                   Sex = clinical_TB$sex, 
                   Smoking_status =clinical_TB$smoking_status)

boxplot_gsva <- as.data.frame(boxplot_gsva)
boxplot_gsva$genesig_D_6 <- as.numeric(boxplot_gsva$genesig_D_6)
boxplot_gsva$genesig_D_7 <- as.numeric(boxplot_gsva$genesig_D_7)


#remove NAs
boxplot_gsva <- boxplot_gsva[-which(is.na(boxplot_gsva$Smoking_status)),]

library(rstatix)

#Get P-values for WITHIN Condition group comparisons (male vs female within each Condition)
stat.table.gsva <- boxplot_gsva  %>%
  group_by(Condition) %>% 
  wilcox_test(genesig_D_7 ~ Smoking_status,
              paired = FALSE,
              p.adjust.method = "BH") %>%
  add_xy_position(x = "Condition")
stat.table.gsva$y.position <- max(boxplot_gsva$genesig_D_7) + 0.05*(max(boxplot_gsva$genesig_D_7))
new.bracket.distance <- 0.2
stat.table.gsva$y.position[c(2,5,8,11)] <- max(boxplot_gsva$genesig_D_7) + new.bracket.distance
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]

# library(scales)
# stat.table.gsva$p[3] <- scientific(stat.table.gsva$p[3], digits = 3)


library(ggpubr)
#Split into two plots, one for HC and one for TB and facet? or ggarrange 
boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(Condition),
  y = as.numeric(genesig_D_7))) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5),
        legend.position = "bottom") +
  
  geom_boxplot(aes(fill = Smoking_status)) +
  
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  #within groups
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 3.8) 

ggsave(boxplotfinal2, 
       file = file.path(clinical_correlation.dir, "timepoint", "smoking_gsva_genesig_D_7.png"),  
       width = 2800,
       height = 2000,
       units = "px" )

table(clinical_TB$condition[clinical$smoking_status == "former"])

### Age BOXPLOT ----------------------------------------------------------------------------------------------------------  

boxplot_gsva=cbind(gsva = gsva_res, Disease = clinical_TB$Disease, 
                   Condition = as.character(clinical_TB$condition), 
                   Age = clinical_TB$age, 
                   Sex = clinical_TB$sex, 
                   Smoking_status =clinical_TB$smoking_status)

boxplot_gsva <- as.data.frame(boxplot_gsva)
boxplot_gsva$genesig_D_6 <- as.numeric(boxplot_gsva$genesig_D_6)
boxplot_gsva$genesig_D_7 <- as.numeric(boxplot_gsva$genesig_D_7)



# table(clinical$sex)
boxplot_gsva$age_category <- cut(as.numeric(boxplot_gsva$Age),
                                 breaks=c(15, 25, 35, 45, 55, 65, Inf),
                                 labels = c("15-24", "25-34", "35-44", "45-54", "55-64", "65+"),
                                 right = FALSE) #range will include lower bound but not the upper bound


#### Group by Age (X AXIS = Condition)  -----------------------------------------------
#Get P-values for WITHIN Condition group comparisons (male vs female within each Condition)
stat.table.gsva <- boxplot_gsva  %>%
  group_by(Condition) %>% 
  wilcox_test(genesig_D_7 ~ age_category,
              paired = FALSE) %>%
  add_xy_position(x = "Condition")
lowest.bracket <- max(boxplot_gsva$genesig_D_7) 
stat.table.gsva <- stat.table.gsva[which(stat.table.gsva$p < 0.05),]
new.bracket.distance <- 0.15
for (i in 1:length(stat.table.gsva$y.position)){
  stat.table.gsva$y.position[i] <-lowest.bracket + new.bracket.distance*i
}

# stat.table.gsva$y.position[4] <- stat.table.gsva$y.position[1]
# stat.table.gsva$y.position[5] <- stat.table.gsva$y.position[2]



boxplotfinal2 <- ggplot(boxplot_gsva, aes(
  x = as.factor(Condition),
  y = as.numeric(genesig_D_7))
) +
  
  theme_bw()+
  
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        title = element_text(size = 11.5)) +
  
  geom_boxplot(aes(fill = age_category)) +
  
  scale_fill_manual(values = c("15-24" = "yellow", 
                               "25-34" = "skyblue", 
                               "35-44" = "purple", 
                               "45-54" = "pink", 
                               "55-64" = "orange", 
                               "65+" = "red"))+
  
  labs(title = "TB 7-Gene Signature",
       caption = str_wrap(paste("Signature:", paste(genesig_D_7, collapse = " ")))) +
  ylab (label = "Enrichment Score") +
  xlab (label = "Condition") +
  
  
  #within groups
  stat_pvalue_manual(stat.table.gsva,
                     label = "p",
                     tip.length = 0.01,
                     size = 4) 

ggsave(boxplotfinal2, file = file.path(clinical_correlation.dir, "timepoint", "age_gsva_genesig_D_7_compareage.png"),  
       width = 3000,
       height = 2000,
       units = "px" )
