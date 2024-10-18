# load area above curve (aac) data
ctrpv2_aac <- readRDS(file = 'model_data/aac_data/CTRPv2.rds')
prism_aac <- readRDS(file = 'model_data/aac_data/PRISM.rds')
gdsc1_aac <- readRDS(file = 'model_data/aac_data/GDSC1.rds')
gdsc2_aac <- readRDS(file = 'model_data/aac_data/GDSC2.rds')

# load omics data
ccle_rnaseq_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RNAseq.rds')
ccle_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RPPA.rds')
ccle_ms_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE MS.rds')
mclp_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/MCLP RPPA.rds')

library(dplyr)   # for data manipulation
library(tidyr)   # for data wrangling
library(purrr)   # for functional programming
library(tibble)  # for tibble manipulation

aac_formatter <- function(aac) {
  aac_long <- aac %>%
    pivot_longer(cols = -c(cell_line, lineage, sublineage), 
                 names_to = "treatment", values_to = "sensitivity")
  
  # retain only iap inhibitors of interest
  aac_long <- aac_long[aac_long$treatment %in% c("AT406", "BIRINAPANT", "LCL161", "GDC0152", "AZD5582", "IAP5620", "IAP7638", "EMBELIN"), ]
  
  # remove rows in aac data where the 'sensitivity' column is NA
  aac_variable <- aac_long[!is.na(aac_long$sensitivity), ]
  
  return(aac_variable)
}

ctrpv2_aac <- aac_formatter(ctrpv2_aac)
prism_aac <- aac_formatter(prism_aac)
gdsc1_aac <- aac_formatter(gdsc1_aac)
gdsc2_aac <- aac_formatter(gdsc2_aac)

# helper function to identify common model context
common <- function(datasets) {
  Reduce(intersect, lapply(datasets, function(dataset) dataset$cell_line))
}

common_model_context <- common(list(ctrpv2_aac, prism_aac, gdsc1_aac, gdsc2_aac)) # n = 30

library(mRMRe)
library(caret)
library(ggplot2)

feature_count_selection <- function(omics, aac, common_model_context, max_features) {
  
  # retrieve unique treatments
  treatments <- unique(aac$treatment)
  
  # store accuracy results 
  accuracy_results <- list()
  
  for (treatment in treatments) {
    aac_subset <- aac[aac$treatment ==  treatment, ]
    if (nrow(aac_subset) < 5) {
      next
    }
    
    common_samples <- intersect(colnames(omics), aac_subset$cell_line)
    if (length(common_samples) == 0) {
      stop("no common samples between omics and sensitivity datasets for", treatment)
    }
    
    # subset omics data and sensitivity data to common samples
    omics_subset <- omics[, common_samples]
    aac_subset <- aac_subset[aac_subset$cell_line %in% common_samples, ]
    
    # remove rows in omics subset with 0 variance
    omics_subset <- omics_subset[apply(omics_subset, 1, var) != 0, ]
    
    # prepare feature matrix (X) and continuous outcome (Y)
    X <- t(omics_subset)
    Y <- as.numeric(aac_subset$sensitivity)
    
    if (nrow(X) < 5) {
      next
    }
    
    # store accuracy results for treatment
    accuracy_result <- data.frame(n_features = integer(), mse = numeric())
    
    common_contexts <- aac_subset$cell_line %in% common_model_context
    common_train_indices <- which(common_contexts)
    
    # separate non-common samples
    different_indices <- which(!common_contexts)
    total_indices <- length(common_train_indices) + length(different_indices)
    
    # determine number of non-common samples for training and testing (0.8:0.2 rule)
    train_size <- round(0.8 * total_indices)  # 80% of total samples for training
    test_size <- total_indices - train_size   # 20% of total samples for testing
    
    train_balance <- train_size - length(common_train_indices)
    
    # Check if there are enough non-common samples for training
    if (train_balance > length(different_indices)) {
      warning("Not enough non-common samples for training, adjusting the split.")
      train_balance <- length(different_indices)  # Set to maximum available non-common samples
    }
    
    set.seed(42) # for reproducibility
    if (train_balance > 0) {
      train_different_contexts <- sample(different_indices, train_balance)
    } else {
      train_different_contexts <- integer(0)  # No non-common samples for training
    }
    
    different_train_indices <- which(different_indices %in% train_different_contexts)
    
    test_contexts <- setdiff(different_indices, train_different_contexts)
    test_indices <- which(different_indices %in% test_contexts)
    
    # final training and test indices
    train_index <- c(common_train_indices, different_train_indices)
    test_index <- test_indices
    
    # subset X and Y to selected features
    X_train <- X[train_index, ]
    X_test <- X[test_index, ]
    Y_train <- Y[train_index]
    Y_test <- Y[test_index]
    
    # iterate over number of features to select (i.e. 2 to max)
    for (n_features in 2:max_features) {
      
      # perform mRMR feature selection
      feature_selection <- mRMR.data(data = data.frame(X, Y))
      selected_features <- mRMR.classic(data = feature_selection, target_indices = ncol(X) + 1, feature_count = n_features)
      
      # subset X to the selected features
      selected_feature_indices <- as.numeric(selected_features@filters[[1]])
      X_selected_train <- X_train[, selected_feature_indices]
      X_selected_test <- X_test[, selected_feature_indices]
      
      # add try-catch block to handle errors in model fitting
      try({
        # train a regression tree model (continuous AAC outcome)
        model <- train(X_selected_train, Y_train, method = "rpart", 
                       trControl = trainControl(method = "cv", number = 5), 
                       tuneGrid = data.frame(cp = 0.01))  # set cp manually
        
        # predict on the test set
        predictions <- predict(model, X_selected_test)
        
        # calculate Mean Squared Error
        mse <- mean((predictions - Y_test)^2)
        
        accuracy_result <- rbind(accuracy_result, data.frame(n_features = n_features, mse = mse))
      }, silent = TRUE)  # continue on error
    }
    
    # store results for current treatment
    accuracy_results[[treatment]] <- accuracy_result
    
    # plot MSE vs number of features for the current treatment
    mse_feature_plot <- ggplot(accuracy_result, aes(x = n_features, y = mse)) +
      geom_line() +
      geom_point() +
      labs(title = paste("MSE vs Number of Features for", treatment), 
           x = "Number of Features", y = "Mean Squared Error") +
      theme_grey()
    
    print(mse_feature_plot)
  }
  
  # return overall accuracy results
  return(accuracy_results)
}


ccle_rnaseq_ctrpv2_fcs <- feature_count_selection(ccle_rnaseq_variables, ctrpv2_aac, common_model_context, 100)
ccle_rnaseq_prism_fcs <- feature_count_selection(ccle_rnaseq_variables, prism_aac, common_model_context, 100)
ccle_rnaseq_gdsc1_fcs <- feature_count_selection(ccle_rnaseq_variables, gdsc1_aac, common_model_context, 100)
ccle_rnaseq_gdsc2_fcs <- feature_count_selection(ccle_rnaseq_variables, gdsc2_aac, common_model_context, 100)

ccle_rppa_ctrpv2_fcs <- feature_count_selection(ccle_rppa_variables, ctrpv2_aac, common_model_context, 100)
ccle_rppa_prism_fcs <- feature_count_selection(ccle_rppa_variables, prism_aac, common_model_context, 100)
ccle_rppa_gdsc1_fcs <- feature_count_selection(ccle_rppa_variables, gdsc1_aac, common_model_context, 100)
ccle_rppa_gdsc2_fcs <- feature_count_selection(ccle_rppa_variables, gdsc2_aac, common_model_context, 100)

ccle_ms_ctrpv2_fcs <- feature_count_selection(ccle_ms_variables, ctrpv2_aac, common_model_context, 100)
ccle_ms_prism_fcs <- feature_count_selection(ccle_ms_variables, prism_aac, common_model_context, 100)
ccle_ms_gdsc1_fcs <- feature_count_selection(ccle_ms_variables, gdsc1_aac, common_model_context, 100)
ccle_ms_gdsc2_fcs <- feature_count_selection(ccle_ms_variables, gdsc2_aac, common_model_context, 100)

mclp_rppa_ctrpv2_fcs <- feature_count_selection(mclp_rppa_variables, ctrpv2_aac, common_model_context, 100)
mclp_rppa_prism_fcs <- feature_count_selection(mclp_rppa_variables, prism_aac, common_model_context, 100)
mclp_rppa_gdsc1_fcs <- feature_count_selection(mclp_rppa_variables, gdsc1_aac, common_model_context, 100)
mclp_rppa_gdsc2_fcs <- feature_count_selection(mclp_rppa_variables, gdsc2_aac, common_model_context, 100)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'mrmr', 'feature_count', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_ctrpv2_fcs, file = 'model_data/mrmr/feature_count/ccle_rnaseq/CTRPv2.rds')
saveRDS(ccle_rnaseq_prism_fcs, file = 'model_data/mrmr/feature_count/ccle_rnaseq/PRISM.rds')
saveRDS(ccle_rnaseq_gdsc1_fcs, file = 'model_data/mrmr/feature_count/ccle_rnaseq/GDSC1.rds')
saveRDS(ccle_rnaseq_gdsc2_fcs, file = 'model_data/mrmr/feature_count/ccle_rnaseq/GDSC2.rds')

saveRDS(ccle_rppa_ctrpv2_fcs, file = 'model_data/mrmr/feature_count/ccle_rppa/CTRPv2.rds')
saveRDS(ccle_rppa_prism_fcs, file = 'model_data/mrmr/feature_count/ccle_rppa/PRISM.rds')
saveRDS(ccle_rppa_gdsc1_fcs, file = 'model_data/mrmr/feature_count/ccle_rppa/GDSC1.rds')
saveRDS(ccle_rppa_gdsc2_fcs, file = 'model_data/mrmr/feature_count/ccle_rppa/GDSC2.rds')

saveRDS(ccle_ms_ctrpv2_fcs, file = 'model_data/mrmr/feature_count/ccle_ms/CTRPv2.rds')
saveRDS(ccle_ms_prism_fcs, file = 'model_data/mrmr/feature_count/ccle_ms/PRISM.rds')
saveRDS(ccle_ms_gdsc1_fcs, file = 'model_data/mrmr/feature_count/ccle_ms/GDSC1.rds')
saveRDS(ccle_ms_gdsc2_fcs, file = 'model_data/mrmr/feature_count/ccle_ms/GDSC2.rds')

saveRDS(mclp_rppa_ctrpv2_fcs, file = 'model_data/mrmr/feature_count/mclp_rppa/CTRPv2.rds')
saveRDS(mclp_rppa_prism_fcs, file = 'model_data/mrmr/feature_count/mclp_rppa/PRISM.rds')
saveRDS(mclp_rppa_gdsc1_fcs, file = 'model_data/mrmr/feature_count/mclp_rppa/GDSC1.rds')
saveRDS(mclp_rppa_gdsc2_fcs, file = 'model_data/mrmr/feature_count/mclp_rppa/GDSC2.rds')

# function to determine optimal feature count
mrmr_selection <- function(fcs, omics, aac, common_model_context) {
  
  selected_features_list <- list()
  
  for (treatment in names(fcs)) {
    treatment_data <- fcs[[treatment]]
    
    # sort by the number of features to ensure the correct order
    treatment_data <- treatment_data[order(treatment_data$n_features), ]
    
    # find the optimal number of features based on the MSE criteria
    optimal_n_features <- treatment_data$n_features[1]  # start with the smallest n_features
    optimal_mse <- treatment_data$mse[1]  # start with the corresponding MSE for the smallest n_features
    
    for (i in 2:nrow(treatment_data)) {
      current_n_features <- treatment_data$n_features[i]
      current_mse <- treatment_data$mse[i]
      
      # calculate the percentage change in MSE
      mse_reduction <- (optimal_mse - current_mse) / optimal_mse
      
      # if MSE reduction is 5% or more, update optimal_n_features
      if (mse_reduction >= 0.05) {
        optimal_n_features <- current_n_features
        optimal_mse <- current_mse
      }
    }
    
    aac_subset <- aac[aac$treatment == treatment, ]
    if (nrow(aac_subset) < 5) {
      next
    }
    
    common_samples <- intersect(colnames(omics), aac_subset$cell_line)
    if (length(common_samples) == 0) {
      stop("no common samples between omics and sensitivity datasets for", treatment)
    }
    
    # subset omics and AAC data to common samples
    omics_subset <- omics[, common_samples]
    aac_subset <- aac_subset[aac_subset$cell_line %in% common_samples, ]
    
    # remove rows in omics subset with 0 variance
    omics_subset <- omics_subset[apply(omics_subset, 1, var) != 0, ]
    
    # prepare feature matrix (X) and continuous outcome (Y)
    X <- t(omics_subset)
    Y <- as.numeric(aac_subset$sensitivity)
    
    if (nrow(X) < 5) {
      next
    }
    
    # perform mRMR feature selection using the optimal_n_features
    feature_selection <- mRMR.data(data = data.frame(X, Y))
    selected_features <- mRMR.classic(data = feature_selection, target_indices = ncol(X) + 1, feature_count = optimal_n_features)
    
    # extract selected feature indices
    selected_feature_indices <- as.numeric(selected_features@filters[[1]])
    
    # extract names for selected feature indices
    selected_feature_names <- rownames(omics_subset)[selected_feature_indices]
    
    # subset X to the selected features
    X_selected <- X[, selected_feature_indices]
    
    # store the treatment, selected features, and additional parameters in the new list
    selected_features_list[[treatment]] <- list(
      optimal_n_features = optimal_n_features,
      optimal_mse = optimal_mse, 
      selected_features = selected_feature_names,
      feature_selection_object = feature_selection  # Store the feature selection object for further investigation
    )
  }
  
  # return new list with the selected features and statistics for each treatment
  return(selected_features_list)
}

ccle_rnaseq_ctrpv2_fs <- mrmr_selection(ccle_rnaseq_ctrpv2_fcs, ccle_rnaseq_variables, ctrpv2_aac, common_model_context)
ccle_rnaseq_prism_fs <- mrmr_selection(ccle_rnaseq_prism_fcs, ccle_rnaseq_variables, prism_aac, common_model_context)
ccle_rnaseq_gdsc1_fs <- mrmr_selection(ccle_rnaseq_gdsc1_fcs, ccle_rnaseq_variables, gdsc1_aac, common_model_context)
ccle_rnaseq_gdsc2_fs <- mrmr_selection(ccle_rnaseq_gdsc2_fcs, ccle_rnaseq_variables, gdsc2_aac, common_model_context)

ccle_rppa_ctrpv2_fs <- mrmr_selection(ccle_rppa_ctrpv2_fcs, ccle_rppa_variables, ctrpv2_aac, common_model_context)
ccle_rppa_prism_fs <- mrmr_selection(ccle_rppa_prism_fcs, ccle_rppa_variables, prism_aac, common_model_context)
ccle_rppa_gdsc1_fs <- mrmr_selection(ccle_rppa_gdsc1_fcs, ccle_rppa_variables, gdsc1_aac, common_model_context)
ccle_rppa_gdsc2_fs <- mrmr_selection(ccle_rppa_gdsc2_fcs, ccle_rppa_variables, gdsc2_aac, common_model_context)

ccle_ms_ctrpv2_fs <- mrmr_selection(ccle_ms_ctrpv2_fcs, ccle_ms_variables, ctrpv2_aac, common_model_context)
ccle_ms_prism_fs <- mrmr_selection(ccle_ms_prism_fcs, ccle_ms_variables, prism_aac, common_model_context)
ccle_ms_gdsc1_fs <- mrmr_selection(ccle_ms_gdsc1_fcs, ccle_ms_variables, gdsc1_aac, common_model_context)
ccle_ms_gdsc2_fs <- mrmr_selection(ccle_ms_gdsc2_fcs, ccle_ms_variables, gdsc2_aac, common_model_context)

mclp_rppa_ctrpv2_fs <- mrmr_selection(mclp_rppa_ctrpv2_fcs, mclp_rppa_variables, ctrpv2_aac, common_model_context)
mclp_rppa_prism_fs <- mrmr_selection(mclp_rppa_prism_fcs, mclp_rppa_variables, prism_aac, common_model_context)
mclp_rppa_gdsc1_fs <- mrmr_selection(mclp_rppa_gdsc1_fcs, mclp_rppa_variables, gdsc1_aac, common_model_context)
mclp_rppa_gdsc2_fs <- mrmr_selection(mclp_rppa_gdsc2_fcs, mclp_rppa_variables, gdsc2_aac, common_model_context)

# load omics features metadata
ccle_rnaseq_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RNAseq.csv')
ccle_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RPPA.csv')
ccle_ms_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE MS.csv')
mclp_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/MCLP RPPA.csv')

# function to map assay_id to gene names based on metadata
gene_mapper <- function(selected_features_list, metadata) {
  
  for (treatment in names(selected_features_list)) {
    
    # get the selected features (assay IDs) for the current treatment
    selected_features <- selected_features_list[[treatment]][["selected_features"]]
    
    # find the matching gene names from metadata where assay_id matches selected_features
    genes <- metadata$gene[match(selected_features, metadata$assay_id)]
    
    # store the corresponding gene names in selected_features_gene
    selected_features_list[[treatment]][["selected_features_gene"]] <- genes
  }
  
  return(selected_features_list)
}

ccle_rnaseq_ctrpv2_fs <- gene_mapper(ccle_rnaseq_ctrpv2_fs, ccle_rnaseq_metadata)
ccle_rnaseq_prism_fs <- gene_mapper(ccle_rnaseq_prism_fs, ccle_rnaseq_metadata)
ccle_rnaseq_gdsc1_fs <- gene_mapper(ccle_rnaseq_gdsc1_fs, ccle_rnaseq_metadata)
ccle_rnaseq_gdsc2_fs <- gene_mapper(ccle_rnaseq_gdsc2_fs, ccle_rnaseq_metadata)

ccle_rppa_ctrpv2_fs <- gene_mapper(ccle_rppa_ctrpv2_fs, ccle_rppa_metadata)
ccle_rppa_prism_fs <- gene_mapper(ccle_rppa_prism_fs, ccle_rppa_metadata)
ccle_rppa_gdsc1_fs <- gene_mapper(ccle_rppa_gdsc1_fs, ccle_rppa_metadata)
ccle_rppa_gdsc2_fs <- gene_mapper(ccle_rppa_gdsc2_fs, ccle_rppa_metadata)

ccle_ms_ctrpv2_fs <- gene_mapper(ccle_ms_ctrpv2_fs, ccle_ms_metadata)
ccle_ms_prism_fs <- gene_mapper(ccle_ms_prism_fs, ccle_ms_metadata)
ccle_ms_gdsc1_fs <- gene_mapper(ccle_ms_gdsc1_fs, ccle_ms_metadata)
ccle_ms_gdsc2_fs <- gene_mapper(ccle_ms_gdsc2_fs, ccle_ms_metadata)

mclp_rppa_ctrpv2_fs <- gene_mapper(mclp_rppa_ctrpv2_fs, mclp_rppa_metadata)
mclp_rppa_prism_fs <- gene_mapper(mclp_rppa_prism_fs, mclp_rppa_metadata)
mclp_rppa_gdsc1_fs <- gene_mapper(mclp_rppa_gdsc1_fs, mclp_rppa_metadata)
mclp_rppa_gdsc2_fs <- gene_mapper(mclp_rppa_gdsc2_fs, mclp_rppa_metadata)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'mrmr', 'feature_selection', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_ctrpv2_fs, file = 'model_data/mrmr/feature_selection/ccle_rnaseq/CTRPv2.rds')
saveRDS(ccle_rnaseq_prism_fs, file = 'model_data/mrmr/feature_selection/ccle_rnaseq/PRISM.rds')
saveRDS(ccle_rnaseq_gdsc1_fs, file = 'model_data/mrmr/feature_selection/ccle_rnaseq/GDSC1.rds')
saveRDS(ccle_rnaseq_gdsc2_fs, file = 'model_data/mrmr/feature_selection/ccle_rnaseq/GDSC2.rds')

saveRDS(ccle_rppa_ctrpv2_fs, file = 'model_data/mrmr/feature_selection/ccle_rppa/CTRPv2.rds')
saveRDS(ccle_rppa_prism_fs, file = 'model_data/mrmr/feature_selection/ccle_rppa/PRISM.rds')
saveRDS(ccle_rppa_gdsc1_fs, file = 'model_data/mrmr/feature_selection/ccle_rppa/GDSC1.rds')
saveRDS(ccle_rppa_gdsc2_fs, file = 'model_data/mrmr/feature_selection/ccle_rppa/GDSC2.rds')

saveRDS(ccle_ms_ctrpv2_fs, file = 'model_data/mrmr/feature_selection/ccle_ms/CTRPv2.rds')
saveRDS(ccle_ms_prism_fs, file = 'model_data/mrmr/feature_selection/ccle_ms/PRISM.rds')
saveRDS(ccle_ms_gdsc1_fs, file = 'model_data/mrmr/feature_selection/ccle_ms/GDSC1.rds')
saveRDS(ccle_ms_gdsc2_fs, file = 'model_data/mrmr/feature_selection/ccle_ms/GDSC2.rds')

saveRDS(mclp_rppa_ctrpv2_fs, file = 'model_data/mrmr/feature_selection/mclp_rppa/CTRPv2.rds')
saveRDS(mclp_rppa_prism_fs, file = 'model_data/mrmr/feature_selection/mclp_rppa/PRISM.rds')
saveRDS(mclp_rppa_gdsc1_fs, file = 'model_data/mrmr/feature_selection/mclp_rppa/GDSC1.rds')
saveRDS(mclp_rppa_gdsc2_fs, file = 'model_data/mrmr/feature_selection/mclp_rppa/GDSC2.rds')

# function to merge different models for the same omics dataset in a list
outputs_merger <- function(...) {
  
  # capture the input lists and their names
  input_lists <- list(...)
  list_names <- as.character(match.call())[-1]  # get the names of the input lists
  
  # initialize an empty list to store the merged results
  merged_list <- list()
  
  # loop over each input list and its corresponding name
  for (i in seq_along(input_lists)) {
    current_list <- input_lists[[i]]  # current feature selection list
    source_name <- list_names[i]  # name of the current list
    
    # extract omics_source, omics, and aac_source from the name
    source_parts <- unlist(strsplit(source_name, "_"))
    omics_source <- source_parts[1]
    omics <- source_parts[2]
    aac_source <- source_parts[3]
    
    # loop over each treatment in the current list
    for (treatment in names(current_list)) {
      
      # get the treatment data and add the source information
      treatment_data <- current_list[[treatment]]
      treatment_data$omics_source <- omics_source
      treatment_data$omics <- omics
      treatment_data$aac_source <- aac_source
      
      # create a unique treatment name by appending the aac_source
      unique_treatment_name <- paste(treatment, "(", aac_source, ")", sep = "")
      
      # add this treatment to the merged list
      merged_list[[unique_treatment_name]] <- treatment_data
    }
  }
  
  # return the merged list
  return(merged_list)
}

ccle_rnaseq_fs_list <- outputs_merger(ccle_rnaseq_ctrpv2_fs, ccle_rnaseq_prism_fs, ccle_rnaseq_gdsc1_fs, ccle_rnaseq_gdsc2_fs)
ccle_rppa_fs_list <- outputs_merger(ccle_rppa_ctrpv2_fs, ccle_rppa_prism_fs, ccle_rppa_gdsc1_fs, ccle_rppa_gdsc2_fs)
ccle_ms_fs_list <- outputs_merger(ccle_ms_ctrpv2_fs, ccle_ms_prism_fs, ccle_ms_gdsc1_fs, ccle_ms_gdsc2_fs)
mclp_rppa_fs_list <- outputs_merger(mclp_rppa_ctrpv2_fs, mclp_rppa_prism_fs, mclp_rppa_gdsc1_fs, mclp_rppa_gdsc2_fs)