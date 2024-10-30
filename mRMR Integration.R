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

library(randomForest)
library(xgboost)
library(caret)

# function to quantify feature importance with tree-based algorithms:
  # 1: random forest (rf) - bagging technique (majority vote)
  # 2: XGboost (xgb)* - boosting technique, performs better if there is class imbalance
feature_importance <- function(omics, aac, common_model_context) {
  
  # initialize empty list to store results
  results <- list(RandomForest = list(), XGBoost = list())
  
  # retrieve unique treatments
  treatments <- unique(aac$treatment)

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
   
    # random forest model feature importance
    rf_model <- randomForest(X, Y, importance = TRUE)
    rf_importance <- importance(rf_model, 
                                type = 1) # calculate mean decrease in accuracy
    rf_importance <- data.frame(Feature = rownames(rf_importance), Importance = rf_importance[, 1])
    results[["RandomForest"]][[treatment]] <- rf_importance
    
    # XGBoost model feature importance
    xgb_data <- xgb.DMatrix(data = as.matrix(X), label = Y)
    
    # determine optimal nrounds using cross-validation with early stopping
    cv_results <- xgb.cv(
      data = xgb_data,
      max_depth = 6,
      eta = 0.1,
      nrounds = 1000, # large initial nrounds to allow early stopping
      objective = "reg:squarederror",
      nfold = 5, # 5-fold cross-validation
      early_stopping_rounds = 10, # stop if no improvement over 10 rounds
      verbose = 0
    )
    
    best_nrounds <- cv_results$best_iteration
    
    xgb_model <- xgboost(data = xgb_data, 
                         max_depth = 6, 
                         eta = 0.1, # smaller value of eta to reduce risk of overfitting (default=0.3)
                         nrounds = best_nrounds, 
                         objective = "reg:squarederror", # regression with squared loss
                         verbose = 0)
    xgb_importance <- xgb.importance(feature_names = colnames(X), model = xgb_model)
    xgb_importance <- xgb_importance %>% select(Feature = Feature, Importance = Gain)
    results[["XGBoost"]][[treatment]] <- xgb_importance
    
  }
  
  return(results)
}

ctrpv2_ccle_rnaseq_ftimp <- feature_importance(ccle_rnaseq_variables, ctrpv2_aac, common_model_context)
prism_ccle_rnaseq_ftimp <- feature_importance(ccle_rnaseq_variables, prism_aac, common_model_context)
gdsc1_ccle_rnaseq_ftimp <- feature_importance(ccle_rnaseq_variables, gdsc1_aac, common_model_context)
gdsc2_ccle_rnaseq_ftimp <- feature_importance(ccle_rnaseq_variables, gdsc2_aac, common_model_context)

omdirs <- c('ccle_rnaseq')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'mrmre', 'feature_importance', dir), recursive = TRUE))
saveRDS(ctrpv2_ccle_rnaseq_ftimp, file = 'model_data/mrmre/feature_importance/ccle_rnaseq/CTRPv2.rds')
saveRDS(prism_ccle_rnaseq_ftimp, file = 'model_data/mrmre/feature_importance/ccle_rnaseq/PRISM.rds')
saveRDS(gdsc1_ccle_rnaseq_ftimp, file = 'model_data/mrmre/feature_importance/ccle_rnaseq/GDSC1.rds')
saveRDS(gdsc2_ccle_rnaseq_ftimp, file = 'model_data/mrmre/feature_importance/ccle_rnaseq/GDSC2.rds')

# function to select optimal feature count based on feature importance selection
  # rf - 1000 features with the highest feature importance (ensuring no -ve and 0 vals)
  # xg - all

feature_count_selection <- function(omics, aac, common_model_context, feature_importance) {
  
  # store accuracy results for each model and treatment
  accuracy_results <- list(RandomForest = list(), XGBoost = list())
  
  # retrieve unique treatments
  treatments <- unique(aac$treatment)
  
  for (treatment in treatments) {
    aac_subset <- aac[aac$treatment == treatment, ]
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
    
    # determine training and testing indices as per previous implementation
    common_contexts <- aac_subset$cell_line %in% common_model_context
    common_train_indices <- which(common_contexts)
    different_indices <- which(!common_contexts)
    
    # split into train and test indices
    set.seed(42)
    train_index <- c(common_train_indices, sample(different_indices, length(different_indices) * 0.8))
    test_index <- setdiff(different_indices, train_index)
    
    X_train <- X[train_index, ]
    X_test <- X[test_index, ]
    Y_train <- Y[train_index]
    Y_test <- Y[test_index]
    
    # rf feature selection
    if (!is.null(feature_importance[["RandomForest"]][[treatment]])) {
      rf_features <- feature_importance[["RandomForest"]][[treatment]]
      
      # select top 1000 features with positive importance
      rf_features <- rf_features[rf_features$Importance > 0, ]
      rf_features <- head(rf_features[order(-rf_features$Importance), "Feature"], 1000)
      rf_indices <- match(rf_features, colnames(X_train))
      rf_indices <- rf_indices[!is.na(rf_indices)]  # remove any NAs in case features are missing in the subset
      
      # evaluate performance for rf-selected features
      accuracy_result_rf <- data.frame(n_features = integer(), mse = numeric())
      for (n_features in 2:length(rf_indices)) {
        X_selected_train <- X_train[, rf_indices[1:n_features]]
        X_selected_test <- X_test[, rf_indices[1:n_features]]
        
        try({
          # train regression tree model on selected features
          model <- train(X_selected_train, Y_train, method = "rpart", 
                         trControl = trainControl(method = "cv", number = 5), 
                         tuneGrid = data.frame(cp = 0.01))
          predictions <- predict(model, X_selected_test)
          mse <- mean((predictions - Y_test)^2)
          
          accuracy_result_rf <- rbind(accuracy_result_rf, data.frame(n_features = n_features, mse = mse))
        }, silent = TRUE)
      }
      
      accuracy_results[["RandomForest"]][[treatment]] <- accuracy_result_rf
      
      # plot for rf-selected features
      rf_plot <- ggplot(accuracy_result_rf, aes(x = n_features, y = mse)) +
        geom_line() + geom_point() +
        labs(title = paste("MSE vs Number of Features for", treatment, "(Random Forest)"), 
             x = "Number of Features", y = "Mean Squared Error") +
        theme_grey()
      print(rf_plot)
    }
    
    # xgb feature selection
    if (!is.null(feature_importance[["XGBoost"]][[treatment]])) {
      xgb_features <- feature_importance[["XGBoost"]][[treatment]]$Feature
      xgb_indices <- match(xgb_features, colnames(X_train))
      xgb_indices <- xgb_indices[!is.na(xgb_indices)]
      
      # evaluate performance for XGBoost-selected features
      accuracy_result_xgb <- data.frame(n_features = integer(), mse = numeric())
      for (n_features in 2:length(xgb_indices)) {
        X_selected_train <- X_train[, xgb_indices[1:n_features]]
        X_selected_test <- X_test[, xgb_indices[1:n_features]]
        
        try({
          # train regression tree model on selected features
          model <- train(X_selected_train, Y_train, method = "rpart", 
                         trControl = trainControl(method = "cv", number = 5), 
                         tuneGrid = data.frame(cp = 0.01))
          predictions <- predict(model, X_selected_test)
          mse <- mean((predictions - Y_test)^2)
          
          accuracy_result_xgb <- rbind(accuracy_result_xgb, data.frame(n_features = n_features, mse = mse))
        }, silent = TRUE)
      }
      
      accuracy_results[["XGBoost"]][[treatment]] <- accuracy_result_xgb
      
      # plot for XGBoost-selected features
      xgb_plot <- ggplot(accuracy_result_xgb, aes(x = n_features, y = mse)) +
        geom_line() + geom_point() +
        labs(title = paste("MSE vs Number of Features for", treatment, "(XGBoost)"), 
             x = "Number of Features", y = "Mean Squared Error") +
        theme_grey()
      print(xgb_plot)
    }
  }
  
  # return overall accuracy results for RandomForest and XGBoost features
  return(accuracy_results)
}

ctrpv2_ccle_rnaseq_fcs <- feature_count_selection(ccle_rnaseq_variables, ctrpv2_aac, common_model_context, ctrpv2_ccle_rnaseq_ftimp)
prism_ccle_rnaseq_fcs <- feature_count_selection(ccle_rnaseq_variables, prism_aac, common_model_context, prism_ccle_rnaseq_ftimp)
gdsc1_ccle_rnaseq_fcs <- feature_count_selection(ccle_rnaseq_variables, gdsc1_aac, common_model_context, gdsc1_ccle_rnaseq_ftimp)
gdsc2_ccle_rnaseq_fcs <- feature_count_selection(ccle_rnaseq_variables, gdsc2_aac, common_model_context, gdsc2_ccle_rnaseq_ftimp)

library(mRMRe)

mrmr_selection <- function(fcs, omics, aac, ftimp) {
  selected_features_list <- list()
  
  for (model in names(fcs)) {
    selected_features_list[[model]] <- list()

    for (treatment in names(fcs[[model]])) {
      treatment_data <- fcs[[model]][[treatment]]
      
      treatment_data <- treatment_data[order(treatment_data$n_features), ]
      
      # determine the optimal number of features based on the MSE criteria
      optimal_n_features <- treatment_data$n_features[1]  # start with the smallest n_features
      optimal_mse <- treatment_data$mse[1]  # start with the corresponding MSE for the smallest n_features
      
      # check if feature importance data exists for the current model and treatment
      if (!is.null(ftimp[[model]][[treatment]])) {
        ftimp_data <- ftimp[[model]][[treatment]]
        
        # select features based on model type
        if (model == "RandomForest") {
          # for RandomForest, select top 1000 features by importance
          top_features <- head(ftimp_data$Feature[order(-ftimp_data$Importance)], 1000)
        } else if (model == "XGBoost") {
          # for XGBoost, include all features
          top_features <- ftimp_data$Feature
        }
        
        # subset AAC data for the current treatment
        aac_subset <- aac[aac$treatment == treatment, ]
        if (nrow(aac_subset) < 5) {
          next  # skip if there are fewer than 5 samples
        }
        
        # find common samples between omics and AAC data for this treatment
        common_samples <- intersect(colnames(omics), aac_subset$cell_line)
        if (length(common_samples) == 0) {
          stop("no common samples between omics and sensitivity datasets for treatment: ", treatment)
        }
        
        # subset omics and AAC data to common samples
        omics_subset <- omics[, common_samples]
        aac_subset <- aac_subset[aac_subset$cell_line %in% common_samples, ]
        
        # filter omics data by selected top features
        omics_subset <- omics_subset[top_features, , drop = FALSE]
        
        # prepare feature matrix (X) and continuous outcome (Y)
        X <- t(omics_subset)
        Y <- as.numeric(aac_subset$sensitivity)
        
        if (nrow(X) < 5) {
          next  # skip if there are fewer than 5 samples after subsetting
        }
        
        # loop through rows to find optimal n_features with MSE reduction ≥ 5%
        for (i in 2:nrow(treatment_data)) {
          current_n_features <- treatment_data$n_features[i]
          current_mse <- treatment_data$mse[i]
          
          # calculate the percentage change in MSE
          mse_reduction <- (optimal_mse - current_mse) / optimal_mse
          
          # update optimal_n_features if MSE reduction is 5% or more
          if (mse_reduction >= 0.05) {
            optimal_n_features <- current_n_features
            optimal_mse <- current_mse
          }
        }
        
        # store results for the current model-treatment
        selected_features_list[[model]][[treatment]] <- list(
          optimal_n_features = optimal_n_features,
          optimal_mse = optimal_mse,
          selected_features = top_features[1:optimal_n_features]
        )
      }
    }
  }
  
  # return the nested list with selected features and statistics for each model-treatment
  return(selected_features_list)
}

ctrpv2_ccle_rnaseq_fs <- mrmr_selection(ctrpv2_ccle_rnaseq_fcs, ccle_rnaseq_variables, ctrpv2_aac, ctrpv2_ccle_rnaseq_ftimp)
prism_ccle_rnaseq_fs <- mrmr_selection(prism_ccle_rnaseq_fcs, ccle_rnaseq_variables, prism_aac, prism_ccle_rnaseq_ftimp)
gdsc1_ccle_rnaseq_fs <- mrmr_selection(gdsc1_ccle_rnaseq_fcs, ccle_rnaseq_variables, gdsc1_aac, gdsc1_ccle_rnaseq_ftimp)
gdsc2_ccle_rnaseq_fs <- mrmr_selection(gdsc2_ccle_rnaseq_fcs, ccle_rnaseq_variables, gdsc2_aac, gdsc2_ccle_rnaseq_ftimp)

omdirs <- c('ccle_rnaseq')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'mrmre', 'feature_selection', dir), recursive = TRUE))
saveRDS(ctrpv2_ccle_rnaseq_fs, file = 'model_data/mrmre/feature_selection/ccle_rnaseq/CTRPv2.rds')
saveRDS(prism_ccle_rnaseq_fs, file = 'model_data/mrmre/feature_selection/ccle_rnaseq/PRISM.rds')
saveRDS(gdsc1_ccle_rnaseq_fs, file = 'model_data/mrmre/feature_selection/ccle_rnaseq/GDSC1.rds')
saveRDS(gdsc2_ccle_rnaseq_fs, file = 'model_data/mrmre/feature_selection/ccle_rnaseq/GDSC2.rds')

# load omics features metadata
ccle_rnaseq_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RNAseq.csv')

# function to map assay_id to gene names based on metadata
gene_mapper <- function(selected_features_list, metadata) {
  
  for (model in names(selected_features_list)) {

    for (treatment in names(selected_features_list[[model]])) {
      
      # retrieve selected features (assay IDs) for the current model and treatment
      selected_features <- selected_features_list[[model]][[treatment]][["selected_features"]]
      
      # find the matching gene names from metadata where assay_id matches selected_features
      genes <- metadata$gene[match(selected_features, metadata$assay_id)]
      
      # store the corresponding gene names in selected_features_gene within the nested list
      selected_features_list[[model]][[treatment]][["selected_features_gene"]] <- genes
    }
  }
  
  return(selected_features_list)
}

ctrpv2_ccle_rnaseq_fs <- gene_mapper(ctrpv2_ccle_rnaseq_fs, ccle_rnaseq_metadata)
prism_ccle_rnaseq_fs <- gene_mapper(prism_ccle_rnaseq_fs, ccle_rnaseq_metadata)
gdsc1_ccle_rnaseq_fs <- gene_mapper(gdsc1_ccle_rnaseq_fs, ccle_rnaseq_metadata)
gdsc2_ccle_rnaseq_fs <- gene_mapper(gdsc2_ccle_rnaseq_fs, ccle_rnaseq_metadata)

# function to merge different models for the same omics dataset in a list
outputs_merger <- function(...) {
  
  # capture the input lists and their names
  input_lists <- list(...)
  list_names <- as.character(match.call())[-1]  # get the names of the input lists
  
  # initialize an empty list to store the merged results with model-based structure
  merged_list <- list()
  
  # loop over each input list and its corresponding name
  for (i in seq_along(input_lists)) {
    current_list <- input_lists[[i]]  # current model-specific feature selection list
    source_name <- list_names[i]  # name of the current list
    
    # extract omics_source, omics, and aac_source from the name
    source_parts <- unlist(strsplit(source_name, "_"))
    omics_source <- source_parts[2]
    omics <- source_parts[3]
    aac_source <- source_parts[1]
    
    # loop over each model in the current list
    for (model in names(current_list)) {
      
      # ensure the model entry exists in the merged list
      if (!model %in% names(merged_list)) {
        merged_list[[model]] <- list()
      }
      
      # loop over each treatment within the current model
      for (treatment in names(current_list[[model]])) {
        
        # get the treatment data and add the source and model information
        treatment_data <- current_list[[model]][[treatment]]
        treatment_data$omics_source <- omics_source
        treatment_data$omics <- omics
        treatment_data$aac_source <- aac_source
        treatment_data$model <- model  # Add model information
        
        # create a unique treatment name by appending the aac_source
        unique_treatment_name <- paste(treatment, "(", aac_source, ")", sep = "")
        
        # add this treatment to the merged list under the correct model
        merged_list[[model]][[unique_treatment_name]] <- treatment_data
      }
    }
  }
  
  # Return the merged list in the format [[model]][[unique_treatment_name]]
  return(merged_list)
}

ccle_rnaseq_fs_list <- outputs_merger(ctrpv2_ccle_rnaseq_fs, prism_ccle_rnaseq_fs, gdsc1_ccle_rnaseq_fs, gdsc2_ccle_rnaseq_fs)

library(ggplot2)
library(dplyr)

performance_comparator <- function(fs_list) {
  
  # prepare an empty data frame to collect plot data
  plot_data <- data.frame(
    unique_treatment_name = character(),
    model = character(),
    optimal_mse = numeric(),
    stringsAsFactors = FALSE
  )
  
  # reshape list into a data frame suitable for plotting
  for (model in names(fs_list)) {
    for (unique_treatment_name in names(fs_list[[model]])) {
      treatment_data <- fs_list[[model]][[unique_treatment_name]]
      plot_data <- rbind(plot_data, data.frame(
        unique_treatment_name = unique_treatment_name,
        model = model,
        optimal_mse = treatment_data[["optimal_mse"]]
      ))
    }
  }
  
  # round optimal_mse values for annotation
  plot_data <- plot_data %>%
    mutate(optimal_mse_label = round(optimal_mse, 5))  # create labels rounded to 5 s.f.
  
  # create the bar plot with ggplot2
  mse_plot <- ggplot(plot_data, aes(x = unique_treatment_name, y = optimal_mse, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(
      title = "MSE by Treatment and Model",
      x = "Treatment (AAC Source)",
      y = "MSE"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # rotate x-axis labels for readability
    scale_fill_viridis_d(option = "D") +  # use the viridis color palette
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +  # add dotted line at y = 0.05 to denote cut-off
    geom_text(aes(label = optimal_mse_label), vjust = -0.3, 
              position = position_dodge(0.9), size = 3)  # annotate values on top of bars
  
  # display the plot
  print(mse_plot)
}

performance_comparator(ccle_rnaseq_fs_list)
