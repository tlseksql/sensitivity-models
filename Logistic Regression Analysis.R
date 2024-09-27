# Load aac data
ctrpv2_aac <- readRDS(file = 'model_data/aac_data/CTRPv2.rds')
prism_aac <- readRDS('model_data/aac_data/PRISM.rds')
gdsc1_aac <- readRDS(file = 'model_data/aac_data/GDSC1.rds')
gdsc2_aac <- readRDS(file = 'model_data/aac_data/GDSC2.rds')

# Load omics data
ccle_rnaseq_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RNAseq.rds')
ccle_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RPPA.rds')
ccle_ms_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE MS.rds')
mclp_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/MCLP RPPA.rds')

library(dplyr)
library(tidyr)
library(tibble)

aac_binarization <- function(aac) {
  aac_long <- aac %>%
    pivot_longer(cols = -c(cell_line, lineage, sublineage), 
                 names_to = "treatment", values_to = "sensitivity")
  
  quartile <- aac_long %>%
    group_by(treatment) %>%
    summarize(third_quartile = quantile(sensitivity, 0.75, na.rm = TRUE)) %>%
    ungroup()
  
  aac_binarized <- aac_long %>%
    left_join(quartile, by = "treatment") %>%
    mutate(sensitivity_bin = ifelse(sensitivity > third_quartile, 1, 0))
  
  return(aac_binarized)
}

ctrpv2_binarized <- aac_binarization(ctrpv2_aac)
prism_binarized <- aac_binarization(prism_aac)
gdsc1_binarized <- aac_binarization(gdsc1_aac)
gdsc2_binarized <- aac_binarization(gdsc2_aac)

library(data.table)
aac_omics_combiner <- function(aac_binarized, omics) {
  aac_dt <- as.data.table(aac_binarized)
  aac_wide <- dcast(aac_dt, cell_line + lineage + sublineage ~ treatment, value.var = "sensitivity_bin")
  aac_wide <- aac_wide[cell_line %in% colnames(omics)]
  
  omics_sub <- omics[, colnames(omics) %in% aac_wide$cell_line, drop = FALSE]
  omics_dt <- as.data.table(t(omics_sub), keep.rownames = "cell_line")
  combined <- merge(aac_wide, omics_dt, by = "cell_line", all.x = TRUE)
  
  return(combined)
}

ctrpv2_ccle_rnaseq_combination <- aac_omics_combiner(ctrpv2_binarized, ccle_rnaseq_variables)
prism_ccle_rnaseq_combination <- aac_omics_combiner(prism_binarized, ccle_rnaseq_variables)
gdsc1_ccle_rnaseq_combination <- aac_omics_combiner(gdsc1_binarized, ccle_rnaseq_variables)
gdsc2_ccle_rnaseq_combination <- aac_omics_combiner(gdsc2_binarized, ccle_rnaseq_variables)

ctrpv2_ccle_rppa_combination <- aac_omics_combiner(ctrpv2_binarized, ccle_rppa_variables)
prism_ccle_rppa_combination <- aac_omics_combiner(prism_binarized, ccle_rppa_variables)
gdsc1_ccle_rppa_combination <- aac_omics_combiner(gdsc1_binarized, ccle_rppa_variables)
gdsc2_ccle_rppa_combination <- aac_omics_combiner(gdsc2_binarized, ccle_rppa_variables)

ctrpv2_ccle_ms_combination <- aac_omics_combiner(ctrpv2_binarized, ccle_ms_variables)
prism_ccle_ms_combination <- aac_omics_combiner(prism_binarized, ccle_ms_variables)
gdsc1_ccle_ms_combination <- aac_omics_combiner(gdsc1_binarized, ccle_ms_variables)
gdsc2_ccle_ms_combination <- aac_omics_combiner(gdsc2_binarized, ccle_ms_variables)

ctrpv2_mclp_rppa_combination <- aac_omics_combiner(ctrpv2_binarized, mclp_rppa_variables)
prism_mclp_rppa_combination <- aac_omics_combiner(prism_binarized, mclp_rppa_variables)
gdsc1_mclp_rppa_combination <- aac_omics_combiner(gdsc1_binarized, mclp_rppa_variables)
gdsc2_mclp_rppa_combination <- aac_omics_combiner(gdsc2_binarized, mclp_rppa_variables)

  # [OPTIONAL] Function to subset omics
    # ! NOTE: Insufficient memory usage to handle parallel processes for ccle_rnaseq (all)
    # Utility: 'max_assays_per_part' capped at 10000 for ccle_rnaseq
omics_combination_splitter <- function(combined, aac_binarized, max_assays_per_part = 10000) {
  treatments <- unique(aac_binarized$treatment)
  treatment_columns <- colnames(combined)[colnames(combined) %in% treatments]
  
  assay_columns <- setdiff(colnames(combined), c("cell_line", "lineage", "sublineage", treatment_columns))
  
  num_parts <- ceiling(length(assay_columns) / max_assays_per_part)
  parts_list <- list()
  
  for (i in seq_len(num_parts)) {
    start_idx <- (i - 1) * max_assays_per_part + 1
    end_idx <- min(i * max_assays_per_part, length(assay_columns))
    
    part_assay_columns <- assay_columns[start_idx:end_idx]
    
    part_combined <- combined[, c("cell_line", "lineage", "sublineage", treatment_columns, part_assay_columns), with = FALSE]
    
    # Name the part and add it to the list
    part_name <- paste0("combination_part", i)
    parts_list[[part_name]] <- part_combined
  }
  
  return(parts_list)
}

library(parallel)
library(e1071)
logistic_regression <- function(combined, aac_binarized) {
  num_cores <- max(1, detectCores() - 1)
  
  treatments <- unique(aac_binarized$treatment)
  treatment_columns <- colnames(combined)[colnames(combined) %in% treatments]
  non_na_treatment_columns <- treatment_columns[colSums(!is.na(combined[, ..treatment_columns])) > 0]
  na_treatments <- setdiff(treatment_columns, non_na_treatment_columns)
  if (length(na_treatments) > 0) {
    cat("No AAC values for:\n", na_treatments, "\n")
  }
  
  assay_columns <- setdiff(colnames(combined), c("cell_line", "lineage", "sublineage", treatment_columns))
  
    # Utility: Small to mid-size data.tables
  # assay_chunks <- split(assay_columns, rep(1:num_cores, length.out = length(assay_columns)))
  
    # Utility: Large data.tables e.g. RNA-seq
      # Chunking manages memory and computational load
  chunk_size <- 25
  assay_chunks <- split(assay_columns, ceiling(seq_along(assay_columns) / chunk_size))
  
  process_chunk <- function(chunk) {
    results <- data.table()
    
    for (assay_id in chunk) {
      id_expression <- combined[[assay_id]]
      valid_indice <- !is.na(id_expression)
      
      if (length(unique(id_expression[valid_indice])) < 2) {
        cat("Insufficient unique id_expression values for", assay_id, ". Skipping assay_id.\n")
        next
      }
      
      binary_subset <- combined[valid_indice, ]
      id_expression <- id_expression[valid_indice]
      
      for (treatment in non_na_treatment_columns) {
        model <- glm(binary_subset[[treatment]] ~ id_expression, family = binomial)
        summary_model <- summary(model)$coefficients
        
        if (nrow(summary_model) < 2) {
          cat("Logistic regression failed for assay", assay_id, "and treatment", treatment, ". Skipping.\n")
          next
        }
        
        beta0 <- summary_model[1, "Estimate"] 
        beta1 <- summary_model[2, "Estimate"]  
        
        odds_ratio <- exp(beta1)
        
        sig_codes <- ifelse(summary_model[2, "Pr(>|z|)"] < 0.001, "***", 
                            ifelse(summary_model[2, "Pr(>|z|)"] < 0.01, "**", 
                                   ifelse(summary_model[2, "Pr(>|z|)"] < 0.05, "*", 
                                          ifelse(summary_model[2, "Pr(>|z|)"] < 0.1, ".", " "))))
        
        mean_probability <- 1 / (1 + exp(-(beta0 + beta1 * mean(id_expression))))
        
        id_range <- seq(min(id_expression), max(id_expression), length.out = 100)
        predicted_probabilities <- 1 / (1 + exp(-(beta0 + beta1 * id_range)))
        
        skewness_probabilities <- skewness(predicted_probabilities)
        
        median_probability <- median(predicted_probabilities)
        max_probability <- max(predicted_probabilities)
        min_probability <- min(predicted_probabilities)
        
        results <- rbindlist(list(results, data.table(
          assay_id = assay_id,
          treatment = treatment,
          intercept = beta0,
          estimate = beta1,
          odds_ratio = odds_ratio,
          std_error = summary_model[2, "Std. Error"],
          z_value = summary_model[2, "z value"],
          p_value = summary_model[2, "Pr(>|z|)"], 
          sig_codes = sig_codes,
          skewness_probabilities = skewness_probabilities,
          mean_probability = mean_probability,
          median_probability = median_probability,
          min_probability = min_probability,
          max_probability = max_probability
        )), use.names = TRUE)
      }
      gc()
    }
    return(results)
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(num_cores)
    on.exit(stopCluster(cl))
    result_list <- parLapply(cl, assay_chunks, process_chunk)
  } else {
    result_list <- mclapply(assay_chunks, process_chunk, mc.cores = num_cores)
  }
  
  results <- rbindlist(result_list, use.names = TRUE)
  
  return(results)
}

logreg_ctrpv2_ccle_rnaseq <- logistic_regression(ctrpv2_ccle_rnaseq_combination, ctrpv2_binarized)
logreg_prism_ccle_rnaseq <- logistic_regression(prism_ccle_rnaseq_combination, prism_binarized)
logreg_gdsc1_ccle_rnaseq <- logistic_regression(gdsc1_ccle_rnaseq_combination, gdsc1_binarized)
logreg_gdsc2_ccle_rnaseq <- logistic_regression(gdsc2_ccle_rnaseq_combination, gdsc2_binarized)

logreg_ctrpv2_ccle_rppa <- logistic_regression(ctrpv2_ccle_rppa_combination, ctrpv2_binarized)
logreg_prism_ccle_rppa <- logistic_regression(prism_ccle_rppa_combination, prism_binarized)
logreg_gdsc1_ccle_rppa <- logistic_regression(gdsc1_ccle_rppa_combination, gdsc1_binarized)
logreg_gdsc2_ccle_rppa <- logistic_regression(gdsc2_ccle_rppa_combination, gdsc2_binarized)

logreg_ctrpv2_ccle_ms <- logistic_regression(ctrpv2_ccle_ms_combination, ctrpv2_binarized)
logreg_prism_ccle_ms <- logistic_regression(prism_ccle_ms_combination, prism_binarized)
logreg_gdsc1_ccle_ms <- logistic_regression(gdsc1_ccle_ms_combination, gdsc1_binarized)
logreg_gdsc2_ccle_ms <- logistic_regression(gdsc2_ccle_ms_combination, gdsc2_binarized)

logreg_ctrpv2_mclp_rppa <- logistic_regression(ctrpv2_mclp_rppa_combination, ctrpv2_binarized)
logreg_prism_mclp_rppa <- logistic_regression(prism_mclp_rppa_combination, prism_binarized)
logreg_gdsc1_mclp_rppa <- logistic_regression(gdsc1_mclp_rppa_combination, gdsc1_binarized)
logreg_gdsc2_mclp_rppa <- logistic_regression(gdsc2_mclp_rppa_combination, gdsc2_binarized)

  # Function to calculate adjusted P value
adj_p_value_calculator <- function(data_table) {
  if (!is.data.table(data_table) || !"p_value" %in% colnames(data_table)) {
    stop("Input must be a data.table containing a 'p_value' column.")
  }
  
  data_table[, adj_p_value := p.adjust(p_value, method = "BH")]
  
  return(data_table)
}

logreg_ctrpv2_ccle_rnaseq <- adj_p_value_calculator(logreg_ctrpv2_ccle_rnaseq)
logreg_prism_ccle_rnaseq <- adj_p_value_calculator(logreg_prism_ccle_rnaseq)
logreg_gdsc1_ccle_rnaseq <- adj_p_value_calculator(logreg_gdsc1_ccle_rnaseq)
logreg_gdsc2_ccle_rnaseq <- adj_p_value_calculator(logreg_gdsc2_ccle_rnaseq)

logreg_ctrpv2_ccle_rppa  <- adj_p_value_calculator(logreg_ctrpv2_ccle_rppa)
logreg_prism_ccle_rppa <- adj_p_value_calculator(logreg_prism_ccle_rppa)
logreg_gdsc1_ccle_rppa <- adj_p_value_calculator(logreg_gdsc1_ccle_rppa) 
logreg_gdsc2_ccle_rppa <- adj_p_value_calculator(logreg_gdsc2_ccle_rppa) 

logreg_ctrpv2_ccle_ms <- adj_p_value_calculator(logreg_ctrpv2_ccle_ms) 
logreg_prism_ccle_ms <- adj_p_value_calculator(logreg_prism_ccle_ms) 
logreg_gdsc1_ccle_ms <- adj_p_value_calculator(logreg_gdsc1_ccle_ms)
logreg_gdsc2_ccle_ms <- adj_p_value_calculator(logreg_gdsc2_ccle_ms)

logreg_ctrpv2_mclp_rppa <- adj_p_value_calculator(logreg_ctrpv2_mclp_rppa) 
logreg_prism_mclp_rppa <- adj_p_value_calculator(logreg_prism_mclp_rppa) 
logreg_gdsc1_mclp_rppa <- adj_p_value_calculator(logreg_gdsc1_mclp_rppa) 
logreg_gdsc2_mclp_rppa <- adj_p_value_calculator(logreg_gdsc2_mclp_rppa) 

  # Logistic regression results files backup
    # !! NOTE: Function takes +/- 48 hours to complete
omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'logreg_analysis', dir), recursive = TRUE))
saveRDS(logreg_ctrpv2_ccle_rnaseq, file = 'model_data/logreg_analysis/ccle_rnaseq/CTRPv2.rds')
saveRDS(logreg_prism_ccle_rnaseq, file = 'model_data/logreg_analysis/ccle_rnaseq/PRISM.rds')
saveRDS(logreg_gdsc1_ccle_rnaseq, file = 'model_data/logreg_analysis/ccle_rnaseq/GDSC1.rds')
saveRDS(logreg_gdsc2_ccle_rnaseq, file = 'model_data/logreg_analysis/ccle_rnaseq/GDSC2.rds')

saveRDS(logreg_ctrpv2_ccle_rppa, file = 'model_data/logreg_analysis/ccle_rppa/CTRPv2.rds')
saveRDS(logreg_prism_ccle_rppa, file = 'model_data/logreg_analysis/ccle_rppa/PRISM.rds')
saveRDS(logreg_gdsc1_ccle_rppa, file = 'model_data/logreg_analysis/ccle_rppa/GDSC1.rds')
saveRDS(logreg_gdsc2_ccle_rppa, file = 'model_data/logreg_analysis/ccle_rppa/GDSC2.rds')

saveRDS(logreg_ctrpv2_ccle_ms, file = 'model_data/logreg_analysis/ccle_ms/CTRPv2.rds')
saveRDS(logreg_prism_ccle_ms, file = 'model_data/logreg_analysis/ccle_ms/PRISM.rds')
saveRDS(logreg_gdsc1_ccle_ms, file = 'model_data/logreg_analysis/ccle_ms/GDSC1.rds')
saveRDS(logreg_gdsc2_ccle_ms, file = 'model_data/logreg_analysis/ccle_ms/GDSC2.rds')

saveRDS(logreg_ctrpv2_mclp_rppa, file = 'model_data/logreg_analysis/mclp_rppa/CTRPv2.rds')
saveRDS(logreg_prism_mclp_rppa, file = 'model_data/logreg_analysis/mclp_rppa/PRISM.rds')
saveRDS(logreg_gdsc1_mclp_rppa, file = 'model_data/logreg_analysis/mclp_rppa/GDSC1.rds')
saveRDS(logreg_gdsc2_mclp_rppa, file = 'model_data/logreg_analysis/mclp_rppa/GDSC2.rds')

logreg_all_ccle_rnaseq <- list(logreg_ctrpv2_ccle_rnaseq = logreg_ctrpv2_ccle_rnaseq, 
                               logreg_prism_ccle_rnaseq = logreg_prism_ccle_rnaseq, 
                               logreg_gdsc1_ccle_rnaseq = logreg_gdsc1_ccle_rnaseq, 
                               logreg_gdsc2_ccle_rnaseq = logreg_gdsc2_ccle_rnaseq)
logreg_all_ccle_rppa <- list(logreg_ctrpv2_ccle_rppa = logreg_ctrpv2_ccle_rppa, 
                               logreg_prism_ccle_rppa = logreg_prism_ccle_rppa, 
                               logreg_gdsc1_ccle_rppa = logreg_gdsc1_ccle_rppa, 
                               logreg_gdsc2_ccle_rppa = logreg_gdsc2_ccle_rppa)
logreg_all_ccle_ms <- list(logreg_ctrpv2_ccle_ms = logreg_ctrpv2_ccle_ms, 
                               logreg_prism_ccle_ms = logreg_prism_ccle_ms, 
                               logreg_gdsc1_ccle_ms = logreg_gdsc1_ccle_ms, 
                               logreg_gdsc2_ccle_ms = logreg_gdsc2_ccle_ms)
logreg_all_mclp_rppa <- list(logreg_ctrpv2_mclp_rppa = logreg_ctrpv2_mclp_rppa, 
                               logreg_prism_mclp_rppa = logreg_prism_mclp_rppa, 
                               logreg_gdsc1_mclp_rppa = logreg_gdsc1_mclp_rppa, 
                               logreg_gdsc2_mclp_rppa = logreg_gdsc2_mclp_rppa)

  # Function to identify assay_ids with significant effect on treatment response
SAI_identifier <- function(list_data, omics_source, omics) {
  criteria <- list(
    ctrpv2 = c("AT406", "BIRINAPANT"),
    prism = c("EMBELIN", "BIRINAPANT", "LCL161", "GDC0152"),
    gdsc1 = c("AZD5582", "EMBELIN", "IAP5620", "IAP7638", "LCL161"),
    gdsc2 = c("AZD5582", "IAP5620", "LCL161")
  )
  
  SAI_list <- list()
  
  for (dataset_name in names(list_data)) {
    name_parts <- strsplit(dataset_name, "_")[[1]]
    aac_source <- name_parts[2]
    omics_source <- name_parts[3]
    omics <- name_parts[4]
    
    if (aac_source %in% names(criteria)) {
      dataset <- list_data[[dataset_name]]
      dataset <- subset(dataset, treatment %in% criteria[[aac_source]])
      
      # Significance definition: P-value < 0.05
      significant_dataset <- subset(dataset, p_value < 0.05)
      
      if (nrow(significant_dataset) > 0) {
        significant_dataset$aac_source <- aac_source
        significant_dataset$omics_source <- omics_source
        significant_dataset$omics <- omics
        SAI_list[[dataset_name]] <- significant_dataset
      }
    }
  }
  
  if (length(SAI_list) == 0) {
    return(data.frame(assay_id = character(), omics_source = character(), omics = character(), stringsAsFactors = FALSE))
  }
  
  SAI_df <- rbindlist(SAI_list, use.names = TRUE)
  return(SAI_df)
}

logreg_IAP_ccle_rnaseq <- SAI_identifier(logreg_all_ccle_rnaseq, "ccle", "rnaseq")
logreg_IAP_ccle_rppa <- SAI_identifier(logreg_all_ccle_rppa, "ccle", "rppa")
logreg_IAP_ccle_ms <- SAI_identifier(logreg_all_ccle_ms, "ccle", "ms")
logreg_IAP_mclp_rppa <- SAI_identifier(logreg_all_mclp_rppa, "mclp", "rppa")

logreg_SAI_ccle_rnaseq_ids <- unique(logreg_IAP_ccle_rnaseq$assay_id) # n = 21004
logreg_SAI_ccle_rppa_ids <- unique(logreg_IAP_ccle_rppa$assay_id) # n = 69
logreg_SAI_ccle_ms_ids <- unique(logreg_IAP_ccle_ms$assay_id) # n = 1425
logreg_SAI_mclp_rppa_ids <- unique(logreg_IAP_mclp_rppa$assay_id) # n = 2

  # Function to analyze unique assay_ids based on odds ratio and median predicted probability
logreg_assay_id_distribution <- function(dt) {
  required_cols <- c("assay_id", "odds_ratio", "median_probability")
  if (!all(required_cols %in% colnames(dt))) {
    stop(paste("The data.table must contain the following columns:", paste(required_cols, collapse = ", ")))
  }
  
  odds_ratio_less_than_one_ids <- dt[odds_ratio < 1, unique(assay_id)]
  odds_ratio_less_than_one_count <- length(odds_ratio_less_than_one_ids)
  
  odds_ratio_greater_than_one_ids <- dt[odds_ratio > 1, unique(assay_id)]
  odds_ratio_greater_than_one_count <- length(odds_ratio_greater_than_one_ids)

  median_prob_greater_than_half_ids <- dt[median_probability > 0.5, unique(assay_id)]
  median_prob_greater_than_half_count <- length(median_prob_greater_than_half_ids)
  
  median_prob_less_than_half_ids <- dt[median_probability < 0.5, unique(assay_id)]
  median_prob_less_than_half_count <- length(median_prob_less_than_half_ids)
  
  assay_id_counts <- dt[, .N, by = assay_id]
  duplicate_assay_ids <- assay_id_counts[N > 1, assay_id]
  
  duplicate_odds_ratio_less_than_one <- dt[odds_ratio < 1 & assay_id %in% duplicate_assay_ids, unique(assay_id)]
  duplicate_odds_ratio_greater_than_one <- dt[odds_ratio > 1 & assay_id %in% duplicate_assay_ids, unique(assay_id)]

  duplicate_median_prob_greater_than_half <- dt[median_probability > 0.5 & assay_id %in% duplicate_assay_ids, unique(assay_id)]
  duplicate_median_prob_less_than_half <- dt[median_probability < 0.5 & assay_id %in% duplicate_assay_ids, unique(assay_id)]
  
  return(list(
    "odds_ratio < 1" = list(
      "count" = odds_ratio_less_than_one_count,
      "assay_ids" = odds_ratio_less_than_one_ids
    ),
    "odds_ratio > 1" = list(
      "count" = odds_ratio_greater_than_one_count,
      "assay_ids" = odds_ratio_greater_than_one_ids
    ),
    "median_probability > 0.5" = list(
      "count" = median_prob_greater_than_half_count,
      "assay_ids" = median_prob_greater_than_half_ids
    ),
    "median_probability < 0.5" = list(
      "count" = median_prob_less_than_half_count,
      "assay_ids" = median_prob_less_than_half_ids
    ),
    "duplicate_assay_ids" = list(
      "overall_count" = length(duplicate_assay_ids),
      "all" = duplicate_assay_ids,
      "odds_ratio < 1" = list(
        "count" = length(duplicate_odds_ratio_less_than_one),
        "assay_ids" = duplicate_odds_ratio_less_than_one
      ),
      "odds_ratio > 1" = list(
        "count" = length(duplicate_odds_ratio_greater_than_one),
        "assay_ids" = duplicate_odds_ratio_greater_than_one
      ),
      "median_probability > 0.5" = list(
        "count" = length(duplicate_median_prob_greater_than_half),
        "assay_ids" = duplicate_median_prob_greater_than_half
      ),
      "median_probability < 0.5" = list(
        "count" = length(duplicate_median_prob_less_than_half),
        "assay_ids" = duplicate_median_prob_less_than_half
      )
    )
  ))
}

logreg_SAI_ccle_rnaseq_ids_distribution <- logreg_assay_id_distribution(logreg_IAP_ccle_rnaseq)
logreg_SAI_ccle_rppa_ids_distribution <- logreg_assay_id_distribution(logreg_IAP_ccle_rppa)
logreg_SAI_ccle_ms_ids_distribution <- logreg_assay_id_distribution(logreg_IAP_ccle_ms)
logreg_SAI_mclp_rppa_ids_distribution <- logreg_assay_id_distribution(logreg_IAP_mclp_rppa)

saveRDS(logreg_SAI_ccle_rnaseq_ids_distribution, file = 'model_data/logreg_analysis/ccle_rnaseq/MPP OR Distribution.rds')
saveRDS(logreg_SAI_ccle_rppa_ids_distribution, file = 'model_data/logreg_analysis/ccle_rppa/MPP OR Distribution.rds')
saveRDS(logreg_SAI_ccle_ms_ids_distribution, file = 'model_data/logreg_analysis/ccle_ms/MPP OR Distribution.rds')
saveRDS(logreg_SAI_mclp_rppa_ids_distribution, file = 'model_data/logreg_analysis/mclp_rppa/MPP OR Distribution.rds')

  # Function to differentiate features into quadrants
quad_assay_id_distribution <- function(dt) {
  if (!("odds_ratio" %in% colnames(dt)) || !("median_probability" %in% colnames(dt))) {
    stop("The data.table must contain 'odds_ratio' and 'median_probability' columns.")
  }
  
  group_1_ids <- dt[odds_ratio < 1 & median_probability > 0.5, unique(assay_id)]
  group_2_ids <- dt[odds_ratio > 1 & median_probability > 0.5, unique(assay_id)]
  group_3_ids <- dt[odds_ratio > 1 & median_probability < 0.5, unique(assay_id)]
  group_4_ids <- dt[odds_ratio < 1 & median_probability < 0.5, unique(assay_id)]
  
  # Identify assay_ids present in more than one group
  all_assay_ids <- c(group_1_ids, group_2_ids, group_3_ids, group_4_ids)
  assay_ids_in_multiple_groups <- all_assay_ids[duplicated(all_assay_ids)]
  
  # Remove assay_ids that are in multiple groups from each group
  group_1_ids_unique <- setdiff(group_1_ids, assay_ids_in_multiple_groups)
  group_2_ids_unique <- setdiff(group_2_ids, assay_ids_in_multiple_groups)
  group_3_ids_unique <- setdiff(group_3_ids, assay_ids_in_multiple_groups)
  group_4_ids_unique <- setdiff(group_4_ids, assay_ids_in_multiple_groups)
  
  group_1_duplicates <- dt[odds_ratio < 1 & median_probability > 0.5, .N, by = assay_id][N > 1, assay_id]
  group_2_duplicates <- dt[odds_ratio > 1 & median_probability > 0.5, .N, by = assay_id][N > 1, assay_id]
  group_3_duplicates <- dt[odds_ratio > 1 & median_probability < 0.5, .N, by = assay_id][N > 1, assay_id]
  group_4_duplicates <- dt[odds_ratio < 1 & median_probability < 0.5, .N, by = assay_id][N > 1, assay_id]
  
  return(list(
    "odds_ratio < 1 & median_probability > 0.5" = list(
      "count" = length(group_1_ids_unique),
      "assay_ids" = group_1_ids_unique,
      "duplicate_assay_ids" = group_1_duplicates
    ),
    "odds_ratio > 1 & median_probability > 0.5" = list(
      "count" = length(group_2_ids_unique),
      "assay_ids" = group_2_ids_unique,
      "duplicate_assay_ids" = group_2_duplicates
    ),
    "odds_ratio > 1 & median_probability < 0.5" = list(
      "count" = length(group_3_ids_unique),
      "assay_ids" = group_3_ids_unique,
      "duplicate_assay_ids" = group_3_duplicates
    ),
    "odds_ratio < 1 & median_probability < 0.5" = list(
      "count" = length(group_4_ids_unique),
      "assay_ids" = group_4_ids_unique,
      "duplicate_assay_ids" = group_4_duplicates
    ),
    "multiple_conditions" = list(
      "assay" = unique(assay_ids_in_multiple_groups)
    )
  ))
}
  
logreg_SAI_ccle_rnaseq_ids_quad <- quad_assay_id_distribution(logreg_IAP_ccle_rnaseq)
logreg_SAI_ccle_rppa_ids_quad <- quad_assay_id_distribution(logreg_IAP_ccle_rppa)
logreg_SAI_ccle_ms_ids_quad <- quad_assay_id_distribution(logreg_IAP_ccle_ms)
logreg_SAI_mclp_rppa_ids_quad <- quad_assay_id_distribution(logreg_IAP_mclp_rppa)

saveRDS(logreg_SAI_ccle_rnaseq_ids_quad, file = 'model_data/logreg_analysis/ccle_rnaseq/MPP OR Quadrant.rds')
saveRDS(logreg_SAI_ccle_rppa_ids_quad, file = 'model_data/logreg_analysis/ccle_rppa/MPP OR Quadrant.rds')
saveRDS(logreg_SAI_ccle_ms_ids_quad, file = 'model_data/logreg_analysis/ccle_ms/MPP OR Quadrant.rds')
saveRDS(logreg_SAI_mclp_rppa_ids_quad, file = 'model_data/logreg_analysis/mclp_rppa/MPP OR Quadrant.rds')

  # Function to subset significant assay_ids across list
SAI_subsetter <- function(results_list, SAI_ids) {
  SAI_list <- list()
  
  for (dataset_name in names(results_list)) {
    name_parts <- strsplit(dataset_name, "_")[[1]]
    aac_source <- name_parts[2]
    omics_source <- name_parts[3]
    omics <- name_parts[4]
    
    dataset <- as.data.table(results_list[[dataset_name]])
    
    SAI_dataset <- dataset[assay_id %in% SAI_ids]
    
    SAI_dataset[, aac_source := aac_source]
    SAI_dataset[, omics_source := omics_source]
    SAI_dataset[, omics := omics]
    
    SAI_list[[dataset_name]] <- SAI_dataset
  }
  
  return(SAI_list)
}

logreg_SAI_ccle_rnaseq <- SAI_subsetter(logreg_all_ccle_rnaseq, logreg_SAI_ccle_rnaseq_ids)
logreg_SAI_ccle_rppa <- SAI_subsetter(logreg_all_ccle_rppa, logreg_SAI_ccle_rppa_ids)
logreg_SAI_ccle_ms <- SAI_subsetter(logreg_all_ccle_ms, logreg_SAI_ccle_ms_ids)
logreg_SAI_mclp_rppa <- SAI_subsetter(logreg_all_mclp_rppa, logreg_SAI_mclp_rppa_ids)

  # Function to convert list to dataset
list_to_dataset <- function(logreg_list) {
  rbindlist(logreg_list, use.names = TRUE)
}

logreg_SAI_ccle_rnaseq <- list_to_dataset(logreg_SAI_ccle_rnaseq)
logreg_SAI_ccle_rppa <- list_to_dataset(logreg_SAI_ccle_rppa)
logreg_SAI_ccle_ms <- list_to_dataset(logreg_SAI_ccle_ms)
logreg_SAI_mclp_rppa <- list_to_dataset(logreg_SAI_mclp_rppa)

# Load aac metadata
aac_metadata <- readRDS(file = 'pset_data/treatment_response/metadata/manual_annotations/COMPLETE LIST.rds')

# Function to identify unique pharmacological classes
pharmacological_class <- aac_metadata %>%
  mutate(
    is_synergy = grepl("SYNERGY: ", manual_annotation, ignore.case = TRUE)
  ) %>%
  group_by(is_synergy) %>%
  do(
    if (.$is_synergy[1]) {
      .
    } else {
      separate_rows(., manual_annotation, sep = ", ")
    }
  ) %>%
  ungroup() %>%
  group_by(manual_annotation) %>%
  summarise(treatment = paste(unique(treatment), collapse = ", ")) %>%
  rename(pharmacological_class = manual_annotation)  %>%
  separate_rows(treatment, sep = ", ")

# Function to clean up pharmacological class metadata
whitespacer_cleaner <- function(metadata) {
  cleaned_metadata <- metadata %>%
    filter(
      !is.na(pharmacological_class) &   # Remove rows where pharmacological_class is NA
        trimws(pharmacological_class) != ""  # Remove rows where pharmacological_class is whitespace
    )
  return(cleaned_metadata)
}

pharmacological_class <- whitespacer_cleaner(pharmacological_class)

reference_treatments <- c("AT406", "BIRINAPANT", "EMBELIN", "LCL161", "GDC0152", "AZD5582", "IAP5620", "IAP7638", "LCL161")

# Function to compute mean probability within pharmacological class
means_calculator <- function(logreg_results, pharmacological_class, reference_treatments) {
  joint <- logreg_results %>%
    inner_join(pharmacological_class, by = "treatment", relationship = "many-to-many")
  
  class_means <- joint %>%
    group_by(pharmacological_class, assay_id) %>%
    summarise(
      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      mean_mean_probability = mean(mean_probability, na.rm = TRUE),
      mean_median_probability = mean(median_probability, na.rm = TRUE),
      .groups = 'drop'
    )
  
  reference_means <- logreg_results %>%
    filter(treatment %in% reference_treatments) %>%
    group_by(assay_id) %>%
    summarise(
      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      mean_mean_probability = mean(mean_probability, na.rm = TRUE),
      mean_median_probability = mean(median_probability, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
  mutate(pharmacological_class = "Reference")
  
  combined_means <- bind_rows(class_means, reference_means)
  return(combined_means)
}

class_means_ccle_rnaseq <- means_calculator(logreg_SAI_ccle_rnaseq, pharmacological_class, reference_treatments)
class_means_ccle_rppa <- means_calculator(logreg_SAI_ccle_rppa, pharmacological_class, reference_treatments)
class_means_ccle_ms <- means_calculator(logreg_SAI_ccle_ms, pharmacological_class, reference_treatments)
class_means_mclp_rppa <- means_calculator(logreg_SAI_mclp_rppa, pharmacological_class, reference_treatments)

# Function to correlate SAI genes expression between pharmacological classes
correlation_calculator <- function(means) {
  reference_means <- means %>%
    filter(pharmacological_class == "Reference") %>%
    select(-pharmacological_class)
  
  class_means <- means %>%
    filter(pharmacological_class != "Reference")
  
  joined_means <- class_means %>%
    inner_join(reference_means, by = "assay_id", suffix = c("_class", "_ref"))
  
  correlations <- joined_means %>%
    group_by(pharmacological_class) %>%
    summarise(
      cor_mean_odds_ratio = cor(mean_odds_ratio_class, mean_odds_ratio_ref, method = "pearson"),
      cor_mean_mean_probability = cor(mean_mean_probability_class, mean_mean_probability_ref, method = "pearson"),
      cor_mean_median_probability = cor(mean_median_probability_class, mean_median_probability_ref, method = "pearson"),
      .groups = 'drop'
    )
  
  return(correlations)
}
       
  # ! NOTE: Correlation with IAP INHIBITORS WILL NOT = 1
    # Compound "IDRONOXIL" is an IAP inhibitor but not a SMAC mimetic
correlations_ccle_rnaseq <- correlation_calculator(class_means_ccle_rnaseq)
correlations_ccle_rppa <- correlation_calculator(class_means_ccle_rppa)
correlations_ccle_ms <- correlation_calculator(class_means_ccle_ms)
correlations_mclp_rppa <- correlation_calculator(class_means_mclp_rppa) # Warning in cor(): Standard deviation is zero

codirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(codirs, function(dir) dir.create(file.path('model_data', 'logreg_correlations', dir), recursive = TRUE))
saveRDS(correlations_ccle_rnaseq , file = 'model_data/logreg_correlations/ccle_rnaseq/Correlations.rds')
saveRDS(correlations_ccle_rppa, file = 'model_data/logreg_correlations/ccle_rppa/Correlations.rds')
saveRDS(correlations_ccle_ms, file = 'model_data/logreg_correlations/ccle_ms/Correlations.rds')
saveRDS(correlations_mclp_rppa, file = 'model_data/logreg_correlations/mclp_rppa/Correlations.rds')
