# load aac data
ctrpv2_aac <- readRDS(file = 'model_data/aac_data/CTRPv2.rds')
prism_aac <- readRDS('model_data/aac_data/PRISM.rds')
gdsc1_aac <- readRDS(file = 'model_data/aac_data/GDSC1.rds')
gdsc2_aac <- readRDS(file = 'model_data/aac_data/GDSC2.rds')

#lLoad omics data
ccle_rnaseq_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RNAseq.rds')
ccle_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE RPPA.rds')
ccle_ms_variables <- readRDS(file = 'model_data/omics_data/selected_features/CCLE MS.rds')
mclp_rppa_variables <- readRDS(file = 'model_data/omics_data/selected_features/MCLP RPPA.rds')

library(dplyr)
library(tidyr)
library(tibble)

aac_formatter <- function(aac) {
  aac_long <- aac %>%
    pivot_longer(cols = -c(cell_line, lineage, sublineage), 
                 names_to = "treatment", values_to = "sensitivity")
  
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

library(parallel)
library(stats)
library(data.table)

linear_regression <- function(combined, aac_data, aac_source, omics_source, omics) {
  num_cores <- max(1, detectCores() - 1)
  
  treatments <- unique(aac_data$treatment)
  treatment_columns <- colnames(combined)[colnames(combined) %in% treatments]
  
  non_na_treatment_columns <- treatment_columns[colSums(!is.na(combined[, ..treatment_columns])) > 0]
  na_treatments <- setdiff(treatment_columns, non_na_treatment_columns)
  if (length(na_treatments) > 0) {
    cat("no AAC values for:\n", na_treatments, "\n")
  }
  
  assay_columns <- setdiff(colnames(combined), c("cell_line", "lineage", "sublineage", treatment_columns))
  
  chunk_size <- 25
  assay_chunks <- split(assay_columns, ceiling(seq_along(assay_columns) / chunk_size))
  
  process_chunk <- function(chunk) {
    results <- data.table()  
    
    for (assay_id in chunk) {
      id_expression <- combined[[assay_id]]
      valid_indice <- !is.na(id_expression)
      
      if (length(unique(id_expression[valid_indice])) < 2) {
        cat("insufficient unique id_expression values for", assay_id, ". skipping assay_id.\n")
        next
      }
      
      subset <- combined[valid_indice, ]
      id_expression <- id_expression[valid_indice]
      
      for (treatment in non_na_treatment_columns) {
        aac_values <- subset[[treatment]]
        
        # fit linear regression model: AAC ~ id_expression
        model <- lm(aac_values ~ id_expression)
        summary_model <- summary(model)
        coefficients <- summary_model$coefficients
        
        if (nrow(summary_model$coefficients) < 2) {
          cat("linear regression failed for assay", assay_id, "and treatment", treatment, ". skipping.\n")
          next
        }
        
        beta0 <- coefficients[1, "Estimate"]  # intercept
        beta1 <- coefficients[2, "Estimate"]  # slope
        
        r_squared <- summary_model$r.squared
        
        # calculate Pearson correlation coefficient (r)
        r_value <- sqrt(r_squared)
        if (coef(model)[2] < 0) {  # adjust sign if slope is negative
          r_value <- -r_value
        }
        
        # Calculate adjusted p-values using Benjamini-Hochberg method
        p_values <- summary_model$coefficients[, "Pr(>|t|)"]
        adj_p_values <- p.adjust(p_values, method = "BH")
        
        results <- rbindlist(list(results, data.table(
          assay_id = assay_id,
          treatment = treatment,
          aac_source = aac_source,
          omics_source = omics_source,
          omics = omics,
          intercept = beta0,
          estimate = beta1,
          r = r_value,  
          r_squared = r_squared,
          adj_r_squared = summary_model$adj.r.squared,
          std_error = summary_model$coefficients[2, "Std. Error"],
          p_value = summary_model$coefficients[2, "Pr(>|t|)"], 
          adj_p_value = adj_p_values[2],
          f_statistic = summary_model$fstatistic[1],
          p_f_statistic = pf(summary_model$fstatistic[1], 
                             summary_model$fstatistic[2], 
                             summary_model$fstatistic[3], 
                             lower.tail = FALSE)
        )), use.names = TRUE)
      }
      gc()  
    }
    
    if (nrow(results) > 0) {
      return(results)
    } else {
      return(NULL) 
    }
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(num_cores)
    on.exit(stopCluster(cl))
    result_list <- parLapply(cl, assay_chunks, process_chunk)
  } else {
    result_list <- mclapply(assay_chunks, process_chunk, mc.cores = num_cores)
  }
  
  result_list <- result_list[!sapply(result_list, is.null)]
  
  if (length(result_list) > 0) {
    results <- rbindlist(result_list, use.names = TRUE)
  } else {
    results <- data.table()  
  }
  
  return(results)
}

ctrpv2_ccle_rnaseq_lr <- linear_regression(combined = ctrpv2_ccle_rnaseq_combination, aac_data = ctrpv2_aac, 
                                           aac_source = "ctrpv2", omics_source = "ccle", omics = "rnaseq")
prism_ccle_rnaseq_lr <- linear_regression(combined = prism_ccle_rnaseq_combination, aac_data = prism_aac, 
                                           aac_source = "prism", omics_source = "ccle", omics = "rnaseq")
gdsc1_ccle_rnaseq_lr <- linear_regression(combined = gdsc1_ccle_rnaseq_combination, aac_data = gdsc1_aac, 
                                          aac_source = "gdsc1", omics_source = "ccle", omics = "rnaseq")
gdsc2_ccle_rnaseq_lr <- linear_regression(combined = gdsc2_ccle_rnaseq_combination, aac_data = gdsc2_aac, 
                                          aac_source = "gdsc2", omics_source = "ccle", omics = "rnaseq")

ctrpv2_ccle_rppa_lr <- linear_regression(combined = ctrpv2_ccle_rppa_combination, aac_data = ctrpv2_aac, 
                                       aac_source = "ctrpv2", omics_source = "ccle", omics = "rppa")
prism_ccle_rppa_lr <- linear_regression(combined = prism_ccle_rppa_combination, aac_data = prism_aac, 
                                      aac_source = "prism", omics_source = "ccle", omics = "rppa")
gdsc1_ccle_rppa_lr <- linear_regression(combined = gdsc1_ccle_rppa_combination, aac_data = gdsc1_aac, 
                                      aac_source = "gdsc1", omics_source = "ccle", omics = "rppa")
gdsc2_ccle_rppa_lr <- linear_regression(combined = gdsc2_ccle_rppa_combination, aac_data = gdsc2_aac, 
                                      aac_source = "gdsc2", omics_source = "ccle", omics = "rppa")

ctrpv2_ccle_ms_lr <- linear_regression(combined = ctrpv2_ccle_ms_combination, aac_data = ctrpv2_aac, 
                                           aac_source = "ctrpv2", omics_source = "ccle", omics = "ms")
prism_ccle_ms_lr <- linear_regression(combined = prism_ccle_ms_combination, aac_data = prism_aac, 
                                          aac_source = "prism", omics_source = "ccle", omics = "ms")
gdsc1_ccle_ms_lr <- linear_regression(combined = gdsc1_ccle_ms_combination, aac_data = gdsc1_aac, 
                                          aac_source = "gdsc1", omics_source = "ccle", omics = "ms")
gdsc2_ccle_ms_lr <- linear_regression(combined = gdsc2_ccle_ms_combination, aac_data = gdsc2_aac, 
                                          aac_source = "gdsc2", omics_source = "ccle", omics = "ms")

ctrpv2_mclp_rppa_lr <- linear_regression(combined = ctrpv2_mclp_rppa_combination, aac_data = ctrpv2_aac, 
                                         aac_source = "ctrpv2", omics_source = "mclp", omics = "rppa")
prism_mclp_rppa_lr <- linear_regression(combined = prism_mclp_rppa_combination, aac_data = prism_aac, 
                                        aac_source = "prism", omics_source = "mclp", omics = "rppa")
gdsc1_mclp_rppa_lr <- linear_regression(combined = gdsc1_mclp_rppa_combination, aac_data = gdsc1_aac, 
                                        aac_source = "gdsc1", omics_source = "mclp", omics = "rppa")
gdsc2_mclp_rppa_lr <- linear_regression(combined = gdsc2_mclp_rppa_combination, aac_data = gdsc2_aac, 
                                        aac_source = "gdsc2", omics_source = "mclp", omics = "rppa")

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'linreg_analysis', dir), recursive = TRUE))
saveRDS(ctrpv2_ccle_rnaseq_lr, file = 'model_data/linreg_analysis/ccle_rnaseq/CTRPv2.rds')
saveRDS(prism_ccle_rnaseq_lr, file = 'model_data/linreg_analysis/ccle_rnaseq/PRISM.rds')
saveRDS(gdsc1_ccle_rnaseq_lr, file = 'model_data/linreg_analysis/ccle_rnaseq/GDSC1.rds')
saveRDS(gdsc2_ccle_rnaseq_lr, file = 'model_data/linreg_analysis/ccle_rnaseq/GDSC2.rds')

saveRDS(ctrpv2_ccle_rppa_lr, file = 'model_data/linreg_analysis/ccle_rppa/CTRPv2.rds')
saveRDS(prism_ccle_rppa_lr, file = 'model_data/linreg_analysis/ccle_rppa/PRISM.rds')
saveRDS(gdsc1_ccle_rppa_lr, file = 'model_data/linreg_analysis/ccle_rppa/GDSC1.rds')
saveRDS(gdsc2_ccle_rppa_lr, file = 'model_data/linreg_analysis/ccle_rppa/GDSC2.rds')

saveRDS(ctrpv2_ccle_ms_lr, file = 'model_data/linreg_analysis/ccle_ms/CTRPv2.rds')
saveRDS(prism_ccle_ms_lr, file = 'model_data/linreg_analysis/ccle_ms/PRISM.rds')
saveRDS(gdsc1_ccle_ms_lr, file = 'model_data/linreg_analysis/ccle_ms/GDSC1.rds')
saveRDS(gdsc2_ccle_ms_lr, file = 'model_data/linreg_analysis/ccle_ms/GDSC2.rds')

saveRDS(ctrpv2_mclp_rppa_lr, file = 'model_data/linreg_analysis/mclp_rppa/CTRPv2.rds')
saveRDS(prism_mclp_rppa_lr, file = 'model_data/linreg_analysis/mclp_rppa/PRISM.rds')
saveRDS(gdsc1_mclp_rppa_lr, file = 'model_data/linreg_analysis/mclp_rppa/GDSC1.rds')
saveRDS(gdsc2_mclp_rppa_lr, file = 'model_data/linreg_analysis/mclp_rppa/GDSC2.rds')

ctrpv2_ccle_rnaseq_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rnaseq/CTRPv2.rds')
prism_ccle_rnaseq_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rnaseq/PRISM.rds')
gdsc1_ccle_rnaseq_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rnaseq/GDSC1.rds')
gdsc2_ccle_rnaseq_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rnaseq/GDSC2.rds')

ctrpv2_ccle_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rppa/CTRPv2.rds')
prism_ccle_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rppa/PRISM.rds')
gdsc1_ccle_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rppa/GDSC1.rds')
gdsc2_ccle_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_rppa/GDSC2.rds')

ctrpv2_ccle_ms_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_ms/CTRPv2.rds')
prism_ccle_ms_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_ms/PRISM.rds')
gdsc1_ccle_ms_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_ms/GDSC1.rds')
gdsc2_ccle_ms_lr <- readRDS(file = 'model_data/linreg_analysis/ccle_ms/GDSC2.rds')

ctrpv2_mclp_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/mclp_rppa/CTRPv2.rds')
prism_mclp_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/mclp_rppa/PRISM.rds')
gdsc1_mclp_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/mclp_rppa/GDSC1.rds')
gdsc2_mclp_rppa_lr <- readRDS(file = 'model_data/linreg_analysis/mclp_rppa/GDSC2.rds')


ccle_rnaseq_lr_list <- list(ctrpv2_ccle_rnaseq_lr = ctrpv2_ccle_rnaseq_lr,
                            prism_ccle_rnaseq_lr = prism_ccle_rnaseq_lr,
                            gdsc1_ccle_rnaseq_lr = gdsc1_ccle_rnaseq_lr,
                            gdsc2_ccle_rnaseq_lr = gdsc2_ccle_rnaseq_lr)

ccle_rppa_lr_list <- list(ctrpv2_ccle_rppa_lr = ctrpv2_ccle_rppa_lr,
                            prism_ccle_rppa_lr = prism_ccle_rppa_lr,
                            gdsc1_ccle_rppa_lr = gdsc1_ccle_rppa_lr,
                            gdsc2_ccle_rppa_lr = gdsc2_ccle_rppa_lr)

ccle_ms_lr_list <- list(ctrpv2_ccle_ms_lr = ctrpv2_ccle_ms_lr,
                            prism_ccle_ms_lr = prism_ccle_ms_lr,
                            gdsc1_ccle_ms_lr = gdsc1_ccle_ms_lr,
                            gdsc2_ccle_ms_lr = gdsc2_ccle_ms_lr)

mclp_rppa_lr_list <- list(ctrpv2_mclp_rppa_lr = ctrpv2_mclp_rppa_lr,
                            prism_mclp_rppa_lr = prism_mclp_rppa_lr,
                            gdsc1_mclp_rppa_lr = gdsc1_mclp_rppa_lr,
                            gdsc2_mclp_rppa_lr = gdsc2_mclp_rppa_lr)

# statistical significance: adj_p_value < 0.05 and p_f_statistic < 0.05
statsig <- function(list_data) {
  treatments <- c("AT406", "BIRINAPANT", "AZD5582", "EMBELIN", "IAP5620", "IAP7638", "LCL161", "GDC0152")
  
  statsig_datasets <- lapply(list_data, function(dataset) {
    filtered_dataset <- dataset[dataset$treatment %in% treatments & dataset$adj_p_value < 0.05 & abs(r) > 0.3, ]
    return(filtered_dataset)
  })
  
  statsig <- rbindlist(statsig_datasets, use.names = TRUE)
  return(statsig)
}

ccle_rnaseq_lr_sigs <- statsig(ccle_rnaseq_lr_list)
ccle_rppa_lr_sigs <- statsig(ccle_rppa_lr_list)
ccle_ms_lr_sigs <- statsig(ccle_ms_lr_list)
mclp_rppa_lr_sigs <- statsig(mclp_rppa_lr_list)

ccle_rnaseq_lr_ids <- unique(ccle_rnaseq_lr_sigs$assay_id) # 10472 (of 23044 obs.)
ccle_rppa_lr_ids <- unique(ccle_rppa_lr_sigs$assay_id) # 90 (of 160 obs.)
ccle_ms_lr_ids <- unique(ccle_ms_lr_sigs$assay_id) # 2093 (of 3534 obs.)
mclp_rppa_lr_ids <- unique(mclp_rppa_lr_sigs$assay_id) # 65 (of 96 obs.)

  # function to subset significant assay_id across all iap inhibitors
statsig_stratifier <- function(dataset) {

  positive_r <- dataset[dataset$r > 0, ]
  negative_r <- dataset[dataset$r < 0, ]
  
  confounding_r <- intersect(positive_r$assay_id, negative_r$assay_id)
  
  # exclude assay_ids with confounding r from each of the sets
  positive_r <- positive_r[!(positive_r$assay_id %in% confounding_r), ]
  negative_r <- negative_r[!(negative_r$assay_id %in% confounding_r), ]

  strata_list <- list(
    positive_r = unique(positive_r$assay_id),  
    negative_r = unique(negative_r$assay_id),  
    confounding_r = confounding_r              
  )
  
  return(strata_list)
}

ccle_rnaseq_lr_sigstrata <- statsig_stratifier(ccle_rnaseq_lr_sigs)
ccle_rppa_lr_sigstrata <- statsig_stratifier(ccle_rppa_lr_sigs)
ccle_ms_lr_sigstrata <- statsig_stratifier(ccle_ms_lr_sigs)
mclp_rppa_lr_sigstrata <- statsig_stratifier(mclp_rppa_lr_sigs)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'linreg_analysis', dir, 'strata'), recursive = TRUE))
saveRDS(ccle_rnaseq_lr_sigstrata, file = 'model_data/linreg_analysis/ccle_rnaseq/strata/significant_strata.rds')
saveRDS(ccle_rppa_lr_sigstrata, file = 'model_data/linreg_analysis/ccle_rppa/strata/significant_strata.rds')
saveRDS(ccle_ms_lr_sigstrata, file = 'model_data/linreg_analysis/ccle_ms/strata/significant_strata.rds')
saveRDS(mclp_rppa_lr_sigstrata, file = 'model_data/linreg_analysis/mclp_rppa/strata/significant_strata.rds')

  # function to subset significant assay_id for iap inhibitors across all treatments
statsig_refer <- function(list_data, lr_ids) {
  statsig_datasets <- lapply(list_data, function(dataset) {
    filtered_dataset <- dataset[dataset$assay_id %in% lr_ids, ]
    return(filtered_dataset)
  })
  statsig_df <- rbindlist(statsig_datasets, use.names = TRUE)
  return(statsig_df)
}

ccle_rnaseq_lr_sigref <- statsig_refer(ccle_rnaseq_lr_list, ccle_rnaseq_lr_ids)
ccle_rppa_lr_sigref <- statsig_refer(ccle_rppa_lr_list, ccle_rppa_lr_ids)
ccle_ms_lr_sigref <- statsig_refer(ccle_ms_lr_list, ccle_ms_lr_ids)
mclp_rppa_lr_sigref <- statsig_refer(mclp_rppa_lr_list, mclp_rppa_lr_ids)

# load aac metadata
aac_metadata <- readRDS(file = 'pset_data/treatment_response/metadata/manual_annotations/COMPLETE LIST.rds')

  # function to identify unique pharmacological (pcl) classes
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

# clean up pharmacological class metadata
whitespacer_cleaner <- function(metadata) {
  cleaned_metadata <- metadata %>%
    filter(
      !is.na(pharmacological_class) &   # remove rows where pharmacological_class is NA
        trimws(pharmacological_class) != ""  # remove rows where pharmacological_class is whitespace
    )
  return(cleaned_metadata)
}

pharmacological_class <- whitespacer_cleaner(pharmacological_class)

# iap reference treatments
iap_refs <- function(data, pharmacological_class) {
  merged_data <- merge(data, pharmacological_class, by = "treatment", all.x = TRUE, allow.cartesian = TRUE)
  iap_inhibitors <-  c("AT406", "BIRINAPANT", "EMBELIN", "LCL161", "GDC0152", "AZD5582", "IAP5620", "IAP7638")
  
  reference_data <- merged_data[merged_data$treatment %in% iap_inhibitors, ]
  
  return(reference_data)
}

ccle_rnaseq_lr_iapref <- iap_refs(ccle_rnaseq_lr_sigref, pharmacological_class)
ccle_rppa_lr_iapref <- iap_refs(ccle_rppa_lr_sigref, pharmacological_class)
ccle_ms_lr_iapref <- iap_refs(ccle_ms_lr_sigref, pharmacological_class)
mclp_rppa_lr_iapref <- iap_refs(mclp_rppa_lr_sigref, pharmacological_class)

library(corrplot)
library(viridis)

# function to corrplot within pharmacological class
intra_corrplotter <- function(data) {
  setDT(data)
  
  # 1: create unique treatment labels by combining treatment and aac_source
  # data[, treatment_label := paste0(treatment, " (", aac_source, ")")]
  
  # 2: pivot data to wide format, using assay_id as rows and treatment_label as columns
    # use mean as the aggregation function for duplicate assay_id/treatment combinations
  
  ## wide_data <- dcast(data, assay_id ~ treatment_label, value.var = "r", fun.aggregate = mean, na.rm = TRUE)
  wide_data <- dcast(data, assay_id ~ treatment, value.var = "r", fun.aggregate = mean, na.rm = TRUE)
  
  # 3: remove rows where all 'r' values are NA
  wide_data <- wide_data[rowSums(is.na(wide_data)) < (ncol(wide_data) - 1)]
  
  # 4: remove columns with zero variance (constant values or all NAs)
  valid_columns <- sapply(wide_data[, -1, with = FALSE], function(col) {
    sd(col, na.rm = TRUE) > 0
  })
  
  if (sum(valid_columns) == 0) {
    stop("no valid columns with variance found for correlation calculation.")
  }
  
  # keep only the columns with variance
  wide_data_filtered <- wide_data[, c(TRUE, valid_columns), with = FALSE]
  
  # 5: calculate the correlation matrix between treatments
  correlation_matrix <- cor(wide_data_filtered[, -1, with = FALSE], use = "pairwise.complete.obs")
  
  # 6: generate a corrplot to visualize the correlations
  corrplot(correlation_matrix, method = "circle",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           col = viridis(256),
           addCoef.col = "black",
           order = "hclust",
           tl.cex = 0.7,
           number.cex = 0.4)
}

ccle_rnaseq_lr_iapplot <- intra_corrplotter(ccle_rnaseq_lr_iapref)
ccle_rppa_lr_iapplot <- intra_corrplotter(ccle_rppa_lr_iapref)
ccle_ms_lr_iapplot <- intra_corrplotter(ccle_ms_lr_iapref)
mclp_rppa_lr_iapplot <- intra_corrplotter(mclp_rppa_lr_iapref)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'linreg_corrplot', 'intra', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_lr_iapplot, file = 'model_data/linreg_corrplot/intra/ccle_rnaseq/merged_correlations.rds')
saveRDS(ccle_rppa_lr_iapplot, file = 'model_data/linreg_corrplot/intra/ccle_rppa/merged_correlations.rds')
saveRDS(ccle_ms_lr_iapplot, file = 'model_data/linreg_corrplot/intra/ccle_ms/merged_correlations.rds')
saveRDS(mclp_rppa_lr_iapplot, file = 'model_data/linreg_corrplot/intra/mclp_rppa/merged_correlations.rds')

# merge withh pharmacological class
all_refs <- function(data, pharmacological_class) {
  merged_data <- merge(data, pharmacological_class, by = "treatment", all.x = TRUE, allow.cartesian = TRUE)
  
  filtered_data <- merged_data[merged_data$treatment != "IDRONOXIL", ]
  
  class_treatment_counts <- filtered_data[, .(unique_treatments = uniqueN(treatment)), by = pharmacological_class]
  valid_classes <- class_treatment_counts[unique_treatments > 1, pharmacological_class]
  final_data <- filtered_data[pharmacological_class %in% valid_classes]
  
  return(final_data)
}

ccle_rnaseq_lr_allref <- all_refs(ccle_rnaseq_lr_sigref, pharmacological_class)
ccle_rppa_lr_allref <- all_refs(ccle_rppa_lr_sigref, pharmacological_class)
ccle_ms_lr_allref <- all_refs(ccle_ms_lr_sigref, pharmacological_class)
mclp_rppa_lr_allref <- all_refs(mclp_rppa_lr_sigref, pharmacological_class)

library(corrplot)
library(viridis)

inter_corrplotter <- function(data) {
  setDT(data)
  
  # 1: Convert to wide format
  wide_data <- dcast(data, assay_id ~ pharmacological_class, value.var = "r", fun.aggregate = mean, na.rm = TRUE)
  
  # 2: Remove rows where all 'r' values are NA
  wide_data <- wide_data[rowSums(is.na(wide_data)) < (ncol(wide_data) - 1)]
  
  # 3: Remove columns with zero variance (constant values or all NAs)
  valid_columns <- sapply(wide_data[, -1, with = FALSE], function(col) {
    sd(col, na.rm = TRUE) > 0
  })
  
  # Ensure valid_columns is only TRUE/FALSE and handle NAs
  valid_columns[is.na(valid_columns)] <- FALSE
  
  if (sum(valid_columns) == 0) {
    stop("No valid columns with variance found for correlation calculation.")
  }
  
  # Keep only the columns with variance
  wide_data_filtered <- wide_data[, c(TRUE, valid_columns), with = FALSE]
  
  # 4: Calculate the correlation matrix between pharmacological classes
  correlation_matrix <- cor(wide_data_filtered[, -1, with = FALSE], use = "pairwise.complete.obs")
  
  # 5: Extract correlations with "IAP INHIBITOR"
  if (!"IAP INHIBITOR" %in% colnames(correlation_matrix)) {
    stop("Pharmacological class 'IAP INHIBITOR' not found in the dataset.")
  }
  
  iap_correlations <- correlation_matrix[,"IAP INHIBITOR", drop = FALSE]
  
  # 6: Filter out classes that contain "SYNERGY:"
  iap_correlations_filtered <- iap_correlations[!grepl("^SYNERGY:", rownames(iap_correlations)), , drop = FALSE]
  
  # 7: Sort the correlations with "IAP INHIBITOR"
  iap_correlations_sorted <- sort(iap_correlations_filtered[,1], decreasing = TRUE, na.last = NA)
  
  # 8: Select the top 10 and bottom 10 correlated classes (excluding IAP INHIBITOR itself)
  top_10 <- head(iap_correlations_sorted, 10)
  bottom_10 <- tail(iap_correlations_sorted, 10)
  
  # 9: Find the 10 classes with correlation closest to 0
  closest_to_zero <- head(iap_correlations_filtered[order(abs(iap_correlations_filtered[,1])), , drop = FALSE], 10)
  
  # 10: Combine the top 10, bottom 10, and closest to 0 classes, ensuring no duplicates
  selected_classes <- unique(c(names(top_10), names(bottom_10), rownames(closest_to_zero)))
  
  # 11: Subset the correlation matrix to only include the selected classes
  selected_matrix <- correlation_matrix[selected_classes, selected_classes]
  
  # 12: Generate a corrplot to visualize the correlations
  corrplot(selected_matrix, method = "circle",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           col = viridis(256),
           addCoef.col = "black",
           order = "hclust",
           tl.cex = 0.7,
           number.cex = 0.4)
}

ccle_rnaseq_lr_pclplot <- inter_corrplotter(ccle_rnaseq_lr_allref)
ccle_rppa_lr_pclplot <- inter_corrplotter(ccle_rppa_lr_allref)
ccle_ms_lr_pclplot <- inter_corrplotter(ccle_ms_lr_allref)
mclp_rppa_lr_pclplot <- inter_corrplotter(mclp_rppa_lr_allref)

omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'linreg_corrplot', 'inter', dir), recursive = TRUE))
saveRDS(ccle_rnaseq_lr_pclplot, file = 'model_data/linreg_corrplot/inter/ccle_rnaseq/filtered_correlations.rds')
saveRDS(ccle_rppa_lr_pclplot, file = 'model_data/linreg_corrplot/inter/ccle_rppa/filtered_correlations.rds')
saveRDS(ccle_ms_lr_pclplot, file = 'model_data/linreg_corrplot/inter/ccle_ms/filtered_correlations.rds')
saveRDS(mclp_rppa_lr_pclplot, file = 'model_data/linreg_corrplot/inter/mclp_rppa/filtered_correlations.rds')
