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
  
  return(aac_long)
}

ctrpv2_aac <- aac_formatter(ctrpv2_aac)
prism_aac <- aac_formatter(prism_aac)
gdsc1_aac <- aac_formatter(gdsc1_aac)
gdsc2_aac <- aac_formatter(gdsc2_aac)

library(data.table)
aac_omics_combiner <- function(aac, omics) {
  aac_dt <- as.data.table(aac)
  aac_wide <- dcast(aac_dt, cell_line + lineage + sublineage ~ treatment, value.var = "sensitivity")
  aac_wide <- aac_wide[cell_line %in% colnames(omics)]
  
  omics_sub <- omics[, colnames(omics) %in% aac_wide$cell_line, drop = FALSE]
  omics_dt <- as.data.table(t(omics_sub), keep.rownames = "cell_line")
  combined <- merge(aac_wide, omics_dt, by = "cell_line", all.x = TRUE)
  
  return(combined)
}

ctrpv2_ccle_rnaseq_combination <- aac_omics_combiner(ctrpv2_aac, ccle_rnaseq_variables)
prism_ccle_rnaseq_combination <- aac_omics_combiner(prism_aac, ccle_rnaseq_variables)
gdsc1_ccle_rnaseq_combination <- aac_omics_combiner(gdsc1_aac, ccle_rnaseq_variables)
gdsc2_ccle_rnaseq_combination <- aac_omics_combiner(gdsc2_aac, ccle_rnaseq_variables)

ctrpv2_ccle_rppa_combination <- aac_omics_combiner(ctrpv2_aac, ccle_rppa_variables)
prism_ccle_rppa_combination <- aac_omics_combiner(prism_aac, ccle_rppa_variables)
gdsc1_ccle_rppa_combination <- aac_omics_combiner(gdsc1_aac, ccle_rppa_variables)
gdsc2_ccle_rppa_combination <- aac_omics_combiner(gdsc2_aac, ccle_rppa_variables)

ctrpv2_ccle_ms_combination <- aac_omics_combiner(ctrpv2_aac, ccle_ms_variables)
prism_ccle_ms_combination <- aac_omics_combiner(prism_aac, ccle_ms_variables)
gdsc1_ccle_ms_combination <- aac_omics_combiner(gdsc1_aac, ccle_ms_variables)
gdsc2_ccle_ms_combination <- aac_omics_combiner(gdsc2_aac, ccle_ms_variables)

ctrpv2_mclp_rppa_combination <- aac_omics_combiner(ctrpv2_aac, mclp_rppa_variables)
prism_mclp_rppa_combination <- aac_omics_combiner(prism_aac, mclp_rppa_variables)
gdsc1_mclp_rppa_combination <- aac_omics_combiner(gdsc1_aac, mclp_rppa_variables)
gdsc2_mclp_rppa_combination <- aac_omics_combiner(gdsc2_aac, mclp_rppa_variables)

ccle_rnaseq_combinations <- list(ctrpv2 = ctrpv2_ccle_rnaseq_combination,
                                 prism = prism_ccle_rnaseq_combination,
                                 gdsc1 = gdsc1_ccle_rnaseq_combination,
                                 gdsc2 = gdsc2_ccle_rnaseq_combination)

ccle_rppa_combinations <- list(ctrpv2 = ctrpv2_ccle_rppa_combination,
                          prism = prism_ccle_rppa_combination,
                          gdsc1 = gdsc1_ccle_rppa_combination,
                          gdsc2 = gdsc2_ccle_rppa_combination)

ccle_ms_combinations <- list(ctrpv2 = ctrpv2_ccle_ms_combination,
                             prism = prism_ccle_ms_combination,
                             gdsc1 = gdsc1_ccle_ms_combination,
                             gdsc2 = gdsc2_ccle_ms_combination)

mclp_rppa_combinations <- list(ctrpv2 = ctrpv2_mclp_rppa_combination,
                               prism = prism_mclp_rppa_combination,
                               gdsc1 = gdsc1_mclp_rppa_combination,
                               gdsc2 = gdsc2_mclp_rppa_combination)

# function to identify most appropriate lambda
library(glmnet)
library(ggplot2)

cross_validation <- function(combinations, omics_variables, plot_cv = TRUE) {
  set.seed(42)  # for reproducibility
  model_results <- list()
  
  for (name in names(combinations)) {
    cat("processing:", name, "\n")
    
    data <- combinations[[name]]
    
    # identify feature and drug columns
    feature_columns <- intersect(colnames(data), rownames(omics_variables))
    drug_columns <- setdiff(colnames(data), c("cell_line", "lineage", "sublineage", rownames(omics_variables)))
    
    dataset_results <- list()
    
    for (drug in drug_columns) {
      cat("processing:", drug, "\n")
      
      # 1: extract omics features and response variable (AAC values)
      omics <- data[, ..feature_columns, with = FALSE]
      aac <- data[[drug]]
      
      # 2: remove rows with NA in AAC values
      non_na_indices <- which(!is.na(aac))
      if (length(non_na_indices) == 0) {
        cat("skipping:", drug, "due to all NA AAC values.\n")
        next
      }
      
      omics <- omics[non_na_indices, , drop = FALSE]  
      aac <- aac[non_na_indices]  
      
      # check if there are enough unique non-NA AAC values for modeling
      if (length(unique(aac)) == 1) {
        cat("skipping:", drug, "due to constant AAC values.\n")
        next
      }
      
      cat("valid AAC data points:", length(non_na_indices), "\n")
      cat("omics data dimensions:", dim(omics), "\n")
      
      # check if there's enough data for training
      if (nrow(omics) < 2) {
        cat("skipping:", drug, "due to insufficient valid omics features for non-NA AAC values.\n")
        next
      }
      
      # 3: prepare data for glmnet (LASSO regression)
      x_train <- as.matrix(omics)  # predictor matrix
      y_train <- aac  # response variable
      
      # 4: perform 5-fold cross-validation using glmnet
      cv_lasso <- tryCatch({
        cv.glmnet(x_train, y_train, alpha = 1, type.measure = "mse", nfolds = 5)  # 5 fold cross-validation
      }, error = function(e) {
        cat("skipping:", drug, "due to glmnet error:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(cv_lasso)) {
        next
      }
      
      best_lambda <- cv_lasso$lambda.min  # optimal lambda from cross-validation
      lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
      
      # 5: store model and performance for the drug
      dataset_results[[drug]] <- list(
        model = lasso_model,
        best_lambda = best_lambda,
        cv_lasso = cv_lasso  # Store the cv.glmnet object for later use (e.g., plotting)
      )
      
      # 6: optionally plot the cross-validation curve
      if (plot_cv) {
        cat("plotting cross-validation curve for", drug, "\n")
        
        # plot the cross-validation curve using base R plot
        plot(cv_lasso)
      }
    }
    
    model_results[[name]] <- dataset_results
  }
  
  return(model_results)
}

ccle_rnaseq_model_cv <- cross_validation(ccle_rnaseq_combinations, ccle_rnaseq_variables)
ccle_rppa_model_cv <- cross_validation(ccle_rppa_combinations, ccle_rppa_variables)
ccle_ms_model_cv <- cross_validation(ccle_ms_combinations, ccle_ms_variables)
mclp_rppa_model_cv <- cross_validation(mclp_rppa_combinations, mclp_rppa_variables)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'lasso', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_model_cv, file = 'model_data/lasso/ccle_rnaseq/cross validation.rds')
saveRDS(ccle_rppa_model_cv, file = 'model_data/lasso/ccle_rppa/cross validation.rds')
saveRDS(ccle_ms_model_cv, file = 'model_data/lasso/ccle_ms/cross validation.rds')
saveRDS(mclp_rppa_model_cv, file = 'model_data/lasso/mclp_rppa/cross validation.rds')

split_data <- function(combinations, split_ratio = 0.8) {
  set.seed(42)  # for reproducibility
  
  # 1: identify common cell lines across all datasets
  cell_lines_list <- lapply(combinations, function(data) data$cell_line)
  common_cell_lines <- Reduce(intersect, cell_lines_list)
  
  split_results <- list()
  
  # define all possible drug columns
  drug_list <- c("AT406", "BIRINAPANT", "GDC0152", "EMBELIN", "LCL161", "AZD5582", "IAP5620", "IAP7638")
  
  for (name in names(combinations)) {
    aac_data <- combinations[[name]]
    
    # 2: find which drug columns are available in the current dataset
    drug_columns <- intersect(drug_list, colnames(aac_data))
    feature_columns <- setdiff(colnames(aac_data), c("cell_line", "lineage", "sublineage", drug_columns))
    
    if (length(drug_columns) == 0) {
      cat("no drug columns found in dataset:", name, "\n")
      next
    }
    
    # 3: split data into common and non-common cell lines
    common_data <- aac_data[aac_data$cell_line %in% common_cell_lines, ]
    non_common_data <- aac_data[!aac_data$cell_line %in% common_cell_lines, ]
    
    drug_split <- list()  # list to store split results for each drug
    
    # process each available drug
    for (drug in drug_columns) {
      cat("processing:", drug, "in", toupper(name), "\n")
      
      # remove rows where AAC values are NA for the specific drug
      valid_data <- aac_data[!is.na(aac_data[[drug]]), ]
      
      # split valid data into common and non-common cell lines
      common_train <- valid_data[valid_data$cell_line %in% common_cell_lines, ]
      non_common_data <- valid_data[!valid_data$cell_line %in% common_cell_lines, ]
      
      num_common <- nrow(common_train)
      
      # calculate total number of cells required for training set (80%)
      total_cells <- nrow(valid_data)
      required_train_size <- round(total_cells * split_ratio)
      
      # calculate non-common cell lines required to fill training set quota
      non_common_cells_needed_for_train <- required_train_size - num_common
      
      # ! checkpoint if common cells fulfill the ratio -> no negative count
      if (non_common_cells_needed_for_train < 0) {
        non_common_cells_needed_for_train <- 0
      }
      
      # split non-common cell lines by lineage and fill the training set
      train_indices <- c()
      test_indices <- c()
      
      for (lineage_type in unique(non_common_data$lineage)) {
        lineage_data <- non_common_data[non_common_data$lineage == lineage_type, ]
        num_rows <- nrow(lineage_data)
        
        # proportionally select how many non-common cells from each lineage are needed for training
        train_size_for_lineage <- min(round(num_rows * (non_common_cells_needed_for_train / nrow(non_common_data))), num_rows)
        
        # sample indices for training and test sets
        lineage_train_indices <- sample(seq_len(num_rows), size = train_size_for_lineage)
        
        # add indices to train and test sets
        train_indices <- c(train_indices, rownames(lineage_data)[lineage_train_indices])
        test_indices <- c(test_indices, rownames(lineage_data)[-lineage_train_indices])
      }
      
      # create train and test sets using calculated indices for non-common data
      train_set <- non_common_data[as.numeric(train_indices), ]
      test_set <- non_common_data[as.numeric(test_indices), ]
      
      # combine common and non-common train sets
      final_train_set <- rbind(train_set, common_train)
      
      # drop all other drug columns except for the current drug in the train and test sets
      final_train_set <- final_train_set[, c("cell_line", "lineage", "sublineage", feature_columns, drug), with = FALSE]
      final_test_set <- test_set[, c("cell_line", "lineage", "sublineage", feature_columns, drug), with = FALSE]
      
      # store split results for the drug, with only the relevant AAC column
      drug_split[[drug]] <- list(train = final_train_set, test = final_test_set)
    }
    
    split_results[[name]] <- drug_split
  }
  
  return(split_results)
}

split_ccle_rnaseq_combinations <- split_data(ccle_rnaseq_combinations)
split_ccle_rppa_combinations <- split_data(ccle_rppa_combinations)
split_ccle_ms_combinations <- split_data(ccle_ms_combinations)
split_mclp_rppa_combinations <- split_data(mclp_rppa_combinations)

# function to fit GLM with lasso regularization
lasso_regression <- function(split_combinations, cv, omics_variables) {
  set.seed(42)  # for reproducibility
  final_models <- list()

  for (dataset_name in names(split_combinations)) {
    cat("Processing dataset:", dataset_name, "\n")
    
    dataset_results <- list()

    for (drug in names(split_combinations[[dataset_name]])) {
      cat("Processing drug:", drug, "\n")
      
      # extract the train and test data for the current drug
      train_data <- split_combinations[[dataset_name]][[drug]]$train
      test_data <- split_combinations[[dataset_name]][[drug]]$test
      
      # identify the omics feature columns (which are shared between data and omics variables)
      feature_columns <- intersect(colnames(train_data), rownames(omics_variables))
      
      # if no cross-validation data for the drug, skip it
      if (!is.null(cv[[dataset_name]][[drug]])) {
        best_lambda <- cv[[dataset_name]][[drug]][["best_lambda"]]
        
        # extract omics features and response (AAC values) for training and test sets
        train_omics <- train_data[, ..feature_columns, with = FALSE]
        test_omics <- test_data[, ..feature_columns, with = FALSE]
        
        train_aac <- train_data[[drug]]
        test_aac <- test_data[[drug]]
        
        # remove rows with NA AAC values in the training data
        valid_train_indices <- which(!is.na(train_aac))
        x_train <- as.matrix(train_omics[valid_train_indices, , drop = FALSE])
        y_train <- train_aac[valid_train_indices]
        
        # ensure there are enough valid rows in the training set
        if (length(valid_train_indices) < 2) {
          cat("Skipping drug:", drug, "due to insufficient valid training data.\n")
          next
        }
        
        # fit the Lasso model using the best lambda from cross-validation
        lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
        
        # remove rows with NA AAC values in the test data
        valid_test_indices <- which(!is.na(test_aac))
        x_test <- as.matrix(test_omics[valid_test_indices, , drop = FALSE])
        actual_test_aac <- test_aac[valid_test_indices]
        
        # ensure there are enough valid rows in the test set
        if (length(valid_test_indices) == 0) {
          cat("Skipping predictions for drug:", drug, "due to no valid test data.\n")
          next
        }
        
        # make predictions on the valid test set
        predictions <- predict(lasso_model, newx = x_test)
        
        # calculate Mean Squared Error (MSE)
        mse <- mean((actual_test_aac - predictions)^2, na.rm = TRUE)
        
        # calculate the variance of the actual values in the test set
        variance_actual <- var(actual_test_aac, na.rm = TRUE)
        
        # calculate the normalized MSE (nMSE)
        nMSE <- if (!is.na(variance_actual) && variance_actual > 0) mse / variance_actual else NA
        
        # store the model, predictions, and performance metrics for the drug
        dataset_results[[drug]] <- list(
          model = lasso_model,
          mse = mse,
          nMSE = nMSE,
          predictions = predictions,
          actual = actual_test_aac,
          best_lambda = best_lambda
        )
      } else {
        cat("Skipping drug:", drug, "due to missing cross-validation results.\n")
      }
    }
    
    # store the results for this dataset
    final_models[[dataset_name]] <- dataset_results
  }
  
  return(final_models)
}
    
   
ccle_rnaseq_model_lasso <- lasso_regression(split_ccle_rnaseq_combinations, ccle_rnaseq_model_cv, ccle_rnaseq_variables)
ccle_rppa_model_lasso <- lasso_regression(split_ccle_rppa_combinations, ccle_rppa_model_cv, ccle_rppa_variables)
ccle_ms_model_lasso <- lasso_regression(split_ccle_ms_combinations, ccle_ms_model_cv, ccle_ms_variables)
mclp_rppa_model_lasso <- lasso_regression(split_mclp_rppa_combinations, mclp_rppa_model_cv, mclp_rppa_variables)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'lasso', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_model_lasso, file = 'model_data/lasso/ccle_rnaseq/model.rds')
saveRDS(ccle_rppa_model_lasso, file = 'model_data/lasso/ccle_rppa/model.rds')
saveRDS(ccle_ms_model_lasso, file = 'model_data/lasso/ccle_ms/model.rds')
saveRDS(mclp_rppa_model_lasso, file = 'model_data/lasso/mclp_rppa/model.rds')

# function to generate residual & qq plots
residuals_plot <- function(lasso_model){
  for (aac_source in names(lasso_model)) {
    for (drug in names(lasso_model[[aac_source]])) {
      predictions <- lasso_model[[aac_source]][[drug]][["predictions"]]
      actual <- lasso_model[[aac_source]][[drug]][["actual"]]
      
      residuals <- actual - predictions
      plot_data <- data.frame(Predictions = as.vector(predictions), Residuals = as.vector(residuals))
      
      residual_plot <- ggplot(data = plot_data, aes(x = Predictions, y = Residuals)) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        labs(title = paste("Residual Plot for", drug, "in", toupper(aac_source)),
             x = "Predicted Values",
             y = "Residuals") +
        theme_grey()
      
      print(residual_plot)
      
      qq_plot <- ggplot(data = plot_data, aes(sample = Residuals)) +
        stat_qq() +
        stat_qq_line(color = 'red') +
        labs(title = paste("Q-Q Plot for", drug, "in", toupper(aac_source)),
             x = "Theoretical Quantiles",
             y = "Sample Quantiles") +
        theme_grey()
      
      print(qq_plot)
    }
  }
}

residuals_plot(ccle_rnaseq_model_lasso)
residuals_plot(ccle_rppa_model_lasso)
residuals_plot(ccle_ms_model_lasso)
residuals_plot(mclp_rppa_model_lasso)

# load omics features metadata
ccle_rnaseq_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RNAseq.csv')
ccle_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RPPA.csv')
ccle_ms_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE MS.csv')
mclp_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/MCLP RPPA.csv')

parameters_extractor <- function(lasso_model, metadata) {
  parameters_list <- list()
  
  for (aac_source in names(lasso_model)) {
    for (drug in names(lasso_model[[aac_source]])) {
      beta_object <- lasso_model[[aac_source]][[drug]][["model"]][["beta"]]
      beta_coefficients <- beta_object@x  # non-zero beta coefficients
      non_zero_indices <- which(beta_object != 0)  # indices of non-zero features
      dimnames_full <- beta_object@Dimnames[[1]]  # full list of feature names
      non_zero_features <- dimnames_full[non_zero_indices]  # fitted features
     
      # map each fitted feature (assay_id) to the corresponding gene from the metadata
      mapped_genes <- sapply(non_zero_features, function(assay_id) {
        gene_name <- metadata$gene[metadata$assay_id == assay_id]
        if (length(gene_name) > 0) {
          return(gene_name)
        } else {
          return(assay_id) # if no match is found, keep the original assay_id
        }
      })
      
      # extract other key parameters
      alpha <- lasso_model[[aac_source]][[drug]][["model"]][["a0"]][["s0"]]  
      lambda <- lasso_model[[aac_source]][[drug]][["model"]][["lambda"]]  
      dev_ratio <- lasso_model[[aac_source]][[drug]][["model"]][["dev.ratio"]] 
      mse <- lasso_model[[aac_source]][[drug]][["mse"]]  
      nMSE <- lasso_model[[aac_source]][[drug]][["nMSE"]]
      predictions <- lasso_model[[aac_source]][[drug]][["predictions"]] 
      actual <- lasso_model[[aac_source]][[drug]][["actual"]]  
      best_lambda <- lasso_model[[aac_source]][[drug]][["best_lambda"]]  
      
      parameters_list[[aac_source]][[drug]] <- list(
        aac_source = aac_source,
        drug = drug,
        beta_coefficients = beta_coefficients,
        non_zero_features = mapped_genes,
        alpha = alpha,
        lambda = lambda,
        dev_ratio = dev_ratio,
        mse = mse,
        nMSE = nMSE,
        predictions = predictions,
        actual = actual,
        best_lambda = best_lambda
      )
    }
  }
  return(parameters_list) 
}

ccle_rnaseq_model_parameters <- parameters_extractor(ccle_rnaseq_model_lasso, ccle_rnaseq_metadata)
ccle_rppa_model_parameters <- parameters_extractor(ccle_rppa_model_lasso, ccle_rppa_metadata)
ccle_ms_model_parameters <- parameters_extractor(ccle_ms_model_lasso, ccle_ms_metadata)
mclp_rppa_model_parameters <- parameters_extractor(mclp_rppa_model_lasso, mclp_rppa_metadata)

library(ggplot2)
library(viridis)

# function to generate coefficient plots
coefficient_plots <- function(parameters_list) {
  
  for (aac_source in names(parameters_list)) {
    for (drug in names(parameters_list[[aac_source]])) {
      
      # extract relevant information
      non_zero_features <- parameters_list[[aac_source]][[drug]][["non_zero_features"]]
      beta_coefficients <- parameters_list[[aac_source]][[drug]][["beta_coefficients"]]
      lambda_value <- parameters_list[[aac_source]][[drug]][["lambda"]]
      alpha_value <- parameters_list[[aac_source]][[drug]][["alpha"]]
      dev_ratio <- parameters_list[[aac_source]][[drug]][["dev_ratio"]]
      
      # skip the drug if there are no non-zero features
      if (length(non_zero_features) == 0) {
        cat("skipping:", drug, "in:", toupper(aac_source), "due to no non-zero features.\n")
        next
      }
      
      coefficient_data <- data.frame(
        Feature = non_zero_features,
        Coefficient = beta_coefficients
      )
      
      # plot beta coefficients for non-zero features
      plot_title <- paste("Coefficient Plot for", drug, "in", toupper(aac_source))
      lambda_label <- paste("λ:", round(lambda_value, 5))
      alpha_label <- paste("α:", round(alpha_value, 5))
      dev_ratio_label <- paste("Dev. Ratio:", round(dev_ratio, 5))
      
      coefficient_plot <- ggplot(coefficient_data, aes(x = reorder(Feature, Coefficient), y = Coefficient)) +
        geom_point(aes(color = Coefficient), size = 4) +  # use dots instead of bars
        scale_color_viridis(option = "D", direction = -1) +  # apply the viridis palette
        coord_flip() +  # flip the axes for horizontal layout
        geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
        labs(title = plot_title, x = "Features", y = "Coefficients") +
        theme_grey() +
        geom_text(aes(label = round(Coefficient, 5)), vjust = -1.5, size = 2.5) +
        annotate("text", x = Inf, y = Inf, label = alpha_label, hjust = 1.05, vjust = 57, color = "black", size = 4) +
        annotate("text", x = Inf, y = Inf, label = lambda_label, hjust = 1.05, vjust = 58.5, color = "black", size = 4) +
        annotate("text", x = Inf, y = Inf, label = dev_ratio_label, hjust = 1.05, vjust = 60, color = "black", size = 4) 
      
      print(coefficient_plot)
    }
  }
}

coefficient_plots(ccle_rnaseq_model_parameters)
coefficient_plots(ccle_rppa_model_parameters)
coefficient_plots(ccle_ms_model_parameters)
coefficient_plots(mclp_rppa_model_parameters)
