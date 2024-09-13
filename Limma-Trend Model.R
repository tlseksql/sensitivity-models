# Load aac data
ctrpv2_aac <- readRDS(file = 'model_data/aac_data/CTRPv2.rds')
prism_aac <- readRDS('model_data/aac_data/PRISM.rds')
gdsc1_aac <- readRDS(file = 'model_data/aac_data/GDSC1.rds')
gdsc2_aac <- readRDS(file = 'model_data/aac_data/GDSC2.rds')

library(dplyr)

  # Function to format aac matrices for limma-trend
aac_limma_formatter <- function(aac_list) {
  matrix_list <- list()
  for (dataset_name in names(aac_list)) {
    data <- as.data.frame(aac_list[[dataset_name]])
    data <- data %>%
      dplyr::select(-lineage, -sublineage)
    rownames(data) <- data$cell_line
    data <- data[, !(names(data) %in% c("cell_line"))]
    data <- as.matrix(data)
    data <- t(data)
    
    for (i in 1:nrow(data)) {
      row_name <- rownames(data)[i]
      matrice <- data[i, , drop = FALSE]
      matrice <- matrice[, !is.na(matrice)]
      matrice <- matrix(matrice, nrow = 1, dimnames = list(row_name, colnames(data)[!is.na(data[i, ])]))
      renamed <- paste0(dataset_name, "_", row_name)
      matrix_list[[renamed]] <- matrice
    }
  }
  return(matrix_list)
}

ctrpv2_aac_list <- list(ctrpv2 = ctrpv2_aac)
prism_aac_list <- list(prism = prism_aac)
gdsc1_aac_list <- list(gdsc1 = gdsc1_aac)
gdsc2_aac_list <- list(gdsc2 = gdsc2_aac)

ctrpv2_aac_matrices <- aac_limma_formatter(ctrpv2_aac_list)
prism_aac_matrices <- aac_limma_formatter(prism_aac_list)
gdsc1_aac_matrices <- aac_limma_formatter(gdsc1_aac_list)
gdsc2_aac_matrices <- aac_limma_formatter(gdsc2_aac_list)

# Load omics data
ccle_rnaseq <- readRDS(file = 'model_data/omics_data/CCLE RNAseq.rds')
ccle_rppa <- readRDS(file = 'model_data/omics_data/CCLE RPPA.rds')
ccle_ms <- readRDS(file = 'model_data/omics_data/CCLE MS.rds')
mclp_rppa <- readRDS(file = 'model_data/omics_data/MCLP RPPA.rds')

# Function to remove zero-variance genes in omics
variance_filter <- function(omics) {
  omics <- na.omit(omics)
  gene_variance <- apply(omics, 1, var) # Calculate (row-wise) variance
  variable_genes <- omics[gene_variance != 0, ]
  return(variable_genes)
}

ccle_rnaseq_variables <- variance_filter(ccle_rnaseq)
ccle_rppa_variables <- variance_filter(ccle_rppa)
ccle_ms_variables <- variance_filter(ccle_ms)
mclp_rppa_variables <- variance_filter(mclp_rppa)

# Linear models for microarray data (limma)-trend
library(edgeR)
library(limma)
library(tibble)

limma_trend <- function(aac, omics, aac_source, treatment) {
  available_samples <- intersect(colnames(aac), colnames(omics))
  if (length(available_samples) == 0) {
    cat("No available samples for drug", treatment, "\n")
    return(NULL)
  }
  
  aac <- aac[, available_samples, drop = FALSE]
  omics <- omics[, available_samples, drop = FALSE]
  
  aac <- t(aac)
  design <- model.matrix(~ aac)
  
  if (nrow(design) <= ncol(design)) {
    cat("Skipping", treatment, "due to insufficient samples for design matrix\n")
    return(NULL)
  }

  fit <- lmFit(omics, design)
  fit <- eBayes(fit, trend = TRUE)
  genes <- topTable(fit, coef = ncol(design), number = nrow(omics), adjust.method = "BH")
  
  genes$aac_source <- aac_source
  genes$treatment <- treatment
  
  return(genes)
}

  # Function to run limma trend model
limma_trend_processor <- function(omics, aac) {
  limma_all_results_list <- list()
  for (dataset_name in names(aac)) {
    dataset <- aac[[dataset_name]]
    parts <- strsplit(dataset_name, "_")[[1]]
    aac_source <- parts[1]
    treatment <- parts[2]
    result <- limma_trend(dataset, omics, aac_source, treatment)
    
    cat("Processing:", treatment, "\n")
    
    if (!is.null(result)) {
      limma_all_results_list[[dataset_name]] <- result
    }
    gc()  # Run garbage collection after each treatment
  }
  return(limma_all_results_list)
}

  # !! NOTE: Do not run all aac matrices together
    # Error: vector memory limit of 18.0 Gb reached, see mem.maxVSize()
limma_ctrpv2_ccle_rnaseq <- limma_trend_processor(ccle_rnaseq_variables, ctrpv2_aac_matrices)
limma_prism_ccle_rnaseq <- limma_trend_processor(ccle_rnaseq_variables, prism_aac_matrices)
limma_gdsc1_ccle_rnaseq <- limma_trend_processor(ccle_rnaseq_variables, gdsc1_aac_matrices)
limma_gdsc2_ccle_rnaseq <- limma_trend_processor(ccle_rnaseq_variables, gdsc2_aac_matrices)

limma_ctrpv2_ccle_rppa <- limma_trend_processor(ccle_rppa_variables, ctrpv2_aac_matrices)
limma_prism_ccle_rppa <- limma_trend_processor(ccle_rppa_variables, prism_aac_matrices)
limma_gdsc1_ccle_rppa <- limma_trend_processor(ccle_rppa_variables, gdsc1_aac_matrices)
limma_gdsc2_ccle_rppa <- limma_trend_processor(ccle_rppa_variables, gdsc2_aac_matrices)

limma_ctrpv2_ccle_ms <- limma_trend_processor(ccle_ms_variables, ctrpv2_aac_matrices)
limma_prism_ccle_ms <- limma_trend_processor(ccle_ms_variables, prism_aac_matrices)
limma_gdsc1_ccle_ms <- limma_trend_processor(ccle_ms_variables, gdsc1_aac_matrices)
limma_gdsc2_ccle_ms <- limma_trend_processor(ccle_ms_variables, gdsc2_aac_matrices)

limma_ctrpv2_mclp_rppa <- limma_trend_processor(mclp_rppa_variables, ctrpv2_aac_matrices)
limma_prism_mclp_rppa <- limma_trend_processor(mclp_rppa_variables, prism_aac_matrices)
limma_gdsc1_mclp_rppa <- limma_trend_processor(mclp_rppa_variables, gdsc1_aac_matrices)
limma_gdsc2_mclp_rppa <- limma_trend_processor(mclp_rppa_variables, gdsc2_aac_matrices)

  # Limma results files backup
# !! NOTE: Function takes +/-6 hours to complete
omdirs <- c('ccle_rnaseq', 'ccle_rppa', 'ccle_ms', 'mclp_rppa')
sapply(omdirs, function(dir) dir.create(file.path('model_data', 'limma_results', dir), recursive = TRUE))
saveRDS(limma_ctrpv2_ccle_rnaseq, file = 'model_data/limma_results/ccle_rnaseq/CTRPv2.rds')
saveRDS(limma_prism_ccle_rnaseq, file = 'model_data/limma_results/ccle_rnaseq/PRISM.rds')
saveRDS(limma_gdsc1_ccle_rnaseq, file = 'model_data/limma_results/ccle_rnaseq/GDSC1.rds')
saveRDS(limma_gdsc2_ccle_rnaseq, file = 'model_data/limma_results/ccle_rnaseq/GDSC2.rds')

saveRDS(limma_ctrpv2_ccle_rppa, file = 'model_data/limma_results/ccle_rppa/CTRPv2.rds')
saveRDS(limma_prism_ccle_rppa, file = 'model_data/limma_results/ccle_rppa/PRISM.rds')
saveRDS(limma_gdsc1_ccle_rppa, file = 'model_data/limma_results/ccle_rppa/GDSC1.rds')
saveRDS(limma_gdsc2_ccle_rppa, file = 'model_data/limma_results/ccle_rppa/GDSC2.rds')

saveRDS(limma_ctrpv2_ccle_ms, file = 'model_data/limma_results/ccle_ms/CTRPv2.rds')
saveRDS(limma_prism_ccle_ms, file = 'model_data/limma_results/ccle_ms/PRISM.rds')
saveRDS(limma_gdsc1_ccle_ms, file = 'model_data/limma_results/ccle_ms/GDSC1.rds')
saveRDS(limma_gdsc2_ccle_ms, file = 'model_data/limma_results/ccle_ms/GDSC2.rds')

saveRDS(limma_ctrpv2_mclp_rppa, file = 'model_data/limma_results/mclp_rppa/CTRPv2.rds')
saveRDS(limma_prism_mclp_rppa, file = 'model_data/limma_results/mclp_rppa/PRISM.rds')
saveRDS(limma_gdsc1_mclp_rppa, file = 'model_data/limma_results/mclp_rppa/GDSC1.rds')
saveRDS(limma_gdsc2_mclp_rppa, file = 'model_data/limma_results/mclp_rppa/GDSC2.rds')

limma_all_ccle_rnaseq <- c(limma_ctrpv2_ccle_rnaseq, limma_prism_ccle_rnaseq, limma_gdsc1_ccle_rnaseq, limma_gdsc2_ccle_rnaseq)
limma_all_ccle_rppa <- c(limma_ctrpv2_ccle_rppa, limma_prism_ccle_rppa, limma_gdsc1_ccle_rppa, limma_gdsc2_ccle_rppa)
limma_all_ccle_ms <- c(limma_ctrpv2_ccle_ms, limma_prism_ccle_ms, limma_gdsc1_ccle_ms, limma_gdsc2_ccle_ms)
limma_all_mclp_rppa <- c(limma_ctrpv2_mclp_rppa, limma_prism_mclp_rppa, limma_gdsc1_mclp_rppa, limma_gdsc2_mclp_rppa)

# Function to identify differentially expressed (DE) genes for IAP inhibitors
  # IAP inhibitor
    # ctrpv2: AT406, BIRINAPANT
    # prism: EMBELIN, BIRINAPANT, LCL161, GDC0152
    # gdsc1: AZD5582, EMBELIN, IAP5620, IAP7638, LCL161
    # gdsc2: AZD5582, IAP5620, LCL161

DE_identifier <- function(list_data, omics_source, omics) {
  criteria <- list(
    ctrpv2 = c("AT406", "BIRINAPANT"),
    prism = c("EMBELIN", "BIRINAPANT", "LCL161", "GDC0152"),
    gdsc1 = c("AZD5582", "EMBELIN", "IAP5620", "IAP7638", "LCL161"),
    gdsc2 = c("AZD5582", "IAP5620", "LCL161")
  )
  
  DE_list <- list()
  
  for (dataset_name in names(list_data)) {
    aac_source <- strsplit(dataset_name, "_")[[1]][1]
    treatment <- strsplit(dataset_name, "_")[[1]][2]
    
    if (aac_source %in% names(criteria) && treatment %in% criteria[[aac_source]]) {
      dataset <- list_data[[dataset_name]]
      # Significance definition: Adjusted P-value < 0.05, |logFC| >= 1
      significant_dataset <- subset(dataset, adj.P.Val < 0.05 & abs(logFC) >= 1)
      if (nrow(significant_dataset) > 0) {
        significant_dataset$assay_id <- rownames(significant_dataset)
        significant_dataset$omics_source <- omics_source
        significant_dataset$omics <- omics
        DE_list[[dataset_name]] <- significant_dataset
      }
    }
  }
  
  if (length(DE_list) > 0) {
    combined_DE_list <- do.call(rbind, DE_list)
    rownames(combined_DE_list) <- NULL
    return(combined_DE_list)
  } else {
    return(data.frame(assay_id = character(), omics_source = character(), omics = character(), stringsAsFactors = FALSE))
  }
}

limma_IAP_ccle_rnaseq <- DE_identifier(limma_all_ccle_rnaseq, "ccle", "rnaseq")
limma_IAP_ccle_rppa <- DE_identifier(limma_all_ccle_rppa, "ccle", "rppa")
limma_IAP_ccle_ms <- DE_identifier(limma_all_ccle_ms, "ccle", "ms") # ! n = 0
limma_IAP_mclp_rppa <- DE_identifier(limma_all_mclp_rppa, "mclp", "rppa")

limma_DE_ccle_rnaseq_ids <- unique(limma_IAP_ccle_rnaseq$assay_id) # n = 4281
limma_DE_ccle_rppa_ids <- unique(limma_IAP_ccle_rppa$assay_id) # n = 4
limma_DE_mclp_rppa_ids <- unique(limma_IAP_mclp_rppa$assay_id) # n = 3

  # Function to differentiate differential expression
DE_id_distribution <- function(df) {
  if (!("logFC" %in% colnames(df)) || !("assay_id" %in% colnames(df))) {
    stop("The data must contain 'logFC' and 'assay_id' columns.")
  }

  negative_assay_ids <- unique(df$assay_id[df$logFC < 0])
  negative_count <- length(negative_assay_ids)

  positive_assay_ids <- unique(df$assay_id[df$logFC > 0])
  positive_count <- length(positive_assay_ids)
  
  assay_id_counts <- table(df$assay_id)  # Count occurrences of each assay_id
  duplicate_assay_ids <- names(assay_id_counts[assay_id_counts > 1])  # Filter those with more than 1 row
  
  duplicate_assay_ids_logFC_less_than_zero <- unique(df[df$logFC < 0 & df$assay_id %in% duplicate_assay_ids, "assay_id"])
  duplicate_assay_ids_logFC_greater_than_zero <- unique(df[df$logFC > 0 & df$assay_id %in% duplicate_assay_ids, "assay_id"])
  
  # Return the results as a list containing both assay_ids and their counts
  return(list(
    "logFC < 0" = list(
      "count" = negative_count,
      "assay_ids" = negative_assay_ids
    ),
    "logFC > 0" = list(
      "count" = positive_count,
      "assay_ids" = positive_assay_ids
    ),
    "duplicate_assay_ids" = list(
      "overall_count" = length(duplicate_assay_ids),
      "all" = duplicate_assay_ids,
      "logFC < 0" = list(
        "count" = length(duplicate_assay_ids_logFC_less_than_zero),
        "assay_ids" = duplicate_assay_ids_logFC_less_than_zero
      ),
      "logFC > 0" = list(
        "count" = length(duplicate_assay_ids_logFC_greater_than_zero),
        "assay_ids" = duplicate_assay_ids_logFC_greater_than_zero
      )
    )
  ))
}

limma_DE_ccle_rnaseq_ids_distribution <- DE_id_distribution(limma_IAP_ccle_rnaseq)
limma_DE_ccle_rppa_ids_distribution <- DE_id_distribution(limma_IAP_ccle_rppa)
limma_DE_mclp_rppa_ids_distribution <- DE_id_distribution(limma_IAP_mclp_rppa)

saveRDS(limma_DE_ccle_rnaseq_ids_distribution, file = 'model_data/limma_results/ccle_rnaseq/logFC Distribution.rds')
saveRDS(limma_DE_ccle_rppa_ids_distribution, file = 'model_data/limma_results/ccle_rppa/logFC Distribution.rds')
saveRDS(limma_DE_mclp_rppa_ids_distribution, file = 'model_data/limma_results/mclp_rppa/logFC Distribution.rds')


  # Function to subset DE genes across list
DE_subsetter <- function(results_list, DE_ids) {
  DE_list <- list()
  for (dataset_name in names(results_list)) {
    dataset <- results_list[[dataset_name]]
    
    DE_dataset <- dataset[rownames(dataset) %in% DE_ids, ]
    DE_list[[dataset_name]] <- DE_dataset
  }
  return(DE_list)
}

limma_DE_ccle_rnaseq <- DE_subsetter(limma_all_ccle_rnaseq, limma_DE_ccle_rnaseq_ids)
limma_DE_ccle_rppa <- DE_subsetter(limma_all_ccle_rppa, limma_DE_ccle_rppa_ids)
limma_DE_mclp_rppa <- DE_subsetter(limma_all_mclp_rppa, limma_DE_mclp_rppa_ids)

# Function to convert list to dataset
list_to_dataset <- function(limma_list) {
  row_names <- lapply(limma_list, rownames)
  dataset <- bind_rows(limma_list)
  dataset$assay_id <- unlist(row_names)
  return(dataset)
}

limma_DE_ccle_rnaseq <- list_to_dataset(limma_DE_ccle_rnaseq)
limma_DE_ccle_rppa <- list_to_dataset(limma_DE_ccle_rppa)
limma_DE_mclp_rppa <- list_to_dataset(limma_DE_mclp_rppa)

# Load aac metadata
aac_metadata <- readRDS(file = 'pset_data/treatment_response/metadata/manual_annotations/COMPLETE LIST.rds')

library(tidyr)
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

  # Function to compute mean AveExpr and logFC within pharmacological class
means_calculator <- function(limma_results, pharmacological_class, reference_treatments) {
  joint <- limma_results %>%
    inner_join(pharmacological_class, by = "treatment", relationship = "many-to-many")
  
  class_means <- joint %>%
    group_by(pharmacological_class, assay_id) %>%
    summarise(
      mean_logFC = mean(logFC, na.rm = TRUE),
      mean_AveExpr = mean(AveExpr, na.rm = TRUE),
      .groups = 'drop'
    )
  
  reference_means <- limma_results %>%
    filter(treatment %in% reference_treatments) %>%
    group_by(assay_id) %>%
    summarise(
      mean_logFC = mean(logFC, na.rm = TRUE),
      mean_AveExpr = mean(AveExpr, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(pharmacological_class = "Reference")
  
  combined_means <- bind_rows(class_means, reference_means)
  return(combined_means)
}

class_means_ccle_rnaseq <- means_calculator(limma_DE_ccle_rnaseq, pharmacological_class, reference_treatments)
class_means_ccle_rppa <- means_calculator(limma_DE_ccle_rppa, pharmacological_class, reference_treatments)
class_means_mclp_rppa <- means_calculator(limma_DE_mclp_rppa, pharmacological_class, reference_treatments)

  # Function to correlate DE genes expression between pharmacological classes
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
      cor_logFC = cor(mean_logFC_class, mean_logFC_ref, method = "pearson"),
      cor_AveExpr = cor(mean_AveExpr_class, mean_AveExpr_ref, method = "pearson"),
      .groups = 'drop'
    )
  return(correlations)
}
  
  # ! NOTE: Correlation with IAP INHIBITORS WILL NOT = 1
    # Compound "IDRONOXIL" is an IAP inhibitor but not a SMAC mimetic
correlations_ccle_rnaseq <- correlation_calculator(class_means_ccle_rnaseq)
correlations_ccle_rppa <- correlation_calculator(class_means_ccle_rppa)
correlations_mclp_rppa <- correlation_calculator(class_means_mclp_rppa)

codirs <- c('ccle_rnaseq', 'ccle_rppa', 'mclp_rppa')
sapply(codirs, function(dir) dir.create(file.path('model_data', 'limma_correlations', dir), recursive = TRUE))
saveRDS(correlations_ccle_rnaseq , file = 'model_data/limma_correlations/ccle_rnaseq/Correlations.rds')
saveRDS(correlations_ccle_rppa, file = 'model_data/limma_correlations/ccle_rppa/Correlations.rds')
saveRDS(correlations_mclp_rppa, file = 'model_data/limma_correlations/mclp_rppa/Correlations.rds')

