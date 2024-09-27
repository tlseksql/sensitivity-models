# Model context [n=152]
model_context <- read.csv("Model Context.csv")

# Drug sensitivity data
ctrpv2_complete <- read.csv("pset_data/treatment_response/data/CTRPv2.csv")
prism_complete <- read.csv("pset_data/treatment_response/data/PRISM.csv")
gdsc1_complete <- read.csv("pset_data/treatment_response/data/GDSC1.csv")
gdsc2_complete <- read.csv("pset_data/treatment_response/data/GDSC2.csv")

library(dplyr)
library(stringr)
library(tidyr)

  # Function to extract aac in model context
aac_extractor <- function(data) {
  row_counter <- 0
  data %>%
    dplyr::select(-X) %>%
    dplyr::filter(
      sapply(seq_len(nrow(data)), function(i) {
        row_counter <<- row_counter + 1
        if (row_counter %% 10000 == 0) {
          message(paste("Processing row:", row_counter))
        }
        cl <- data$cell_line[i]
        any(sapply(model_context$cell_line, function(mc) cl %in% str_split(mc, " ")[[1]]))
      })
    ) %>%
    group_by(cell_line, treatment) %>%
    summarise(aac = mean(aac, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = treatment, values_from = aac) %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    left_join(model_context %>% dplyr::select(cell_line, lineage, sublineage),
              by = "cell_line") %>%
    dplyr::select(cell_line, lineage, sublineage, everything())
}

ctrpv2_aac <- aac_extractor(ctrpv2_complete)
prism_aac <- aac_extractor(prism_complete)
gdsc1_aac <- aac_extractor(gdsc1_complete)
gdsc2_aac <- aac_extractor(gdsc2_complete)

# AAC files backup
  # !! NOTE: Function to extract aac takes +/-0.75 hours to complete
dir.create(file.path('model_data', 'aac_data'), recursive = TRUE)
saveRDS(ctrpv2_aac, file = 'model_data/aac_data/CTRPv2.rds')
saveRDS(prism_aac, file = 'model_data/aac_data/PRISM.rds')
saveRDS(gdsc1_aac, file = 'model_data/aac_data/GDSC1.rds')
saveRDS(gdsc2_aac, file = 'model_data/aac_data/GDSC2.rds')


# Omics data
ccle_rnaseq_complete <- read.csv("pset_data/omics/data/CCLE RNAseq.csv")
ccle_rppa_complete <- read.csv("pset_data/omics/data/CCLE RPPA.csv")
ccle_ms_complete <- read.csv("pset_data/omics/data/CCLE MS.csv")
mclp_rppa_complete <- read.csv("pset_data/omics/data/MCLP RPPA.csv")

  # Function to extract omics data in model context
omics_contextor <- function(data, model_context) {
  rownames(data) <- data$X
  data <- data %>%
    dplyr::select(-X)
  data_matrix <- as.matrix(data)
  unique_cell_lines <- unique(unlist(str_split(model_context$cell_line, " ")))
  relevant_columns <- colnames(data_matrix)[colnames(data_matrix) %in% unique_cell_lines]
  contextualized_data <- data_matrix[, relevant_columns, drop = FALSE]
  return(contextualized_data)
}

ccle_rnaseq <- omics_contextor(ccle_rnaseq_complete, model_context)
ccle_rppa <- omics_contextor(ccle_rppa_complete, model_context)
ccle_ms <- omics_contextor(ccle_ms_complete, model_context)
mclp_rppa <- omics_contextor(mclp_rppa_complete, model_context)

  # OMICS files backup
dir.create(file.path('model_data', 'omics_data'), recursive = TRUE)
saveRDS(ccle_rnaseq, file = 'model_data/omics_data/CCLE RNAseq.rds')
saveRDS(ccle_rppa, file = 'model_data/omics_data/CCLE RPPA.rds')
saveRDS(ccle_ms, file = 'model_data/omics_data/CCLE MS.rds')
saveRDS(mclp_rppa, file = 'model_data/omics_data/MCLP RPPA.rds')

  # Function to select for "protein-coding" genes in RNA-seq
protein_coding_filter <- function(omics, metadata){
  protein_coding_metadata <- metadata[metadata$gene_type == "protein_coding", ]
  protein_coding_genes <- protein_coding_metadata$assay_id
  
  omics_filtered <- omics[rownames(omics) %in% protein_coding_genes, ]
  return(omics_filtered)
}
  
ccle_rnaseq_pc <- protein_coding_filter(ccle_rnaseq, ccle_rnaseq_metadata)

  # Function to filter zero variance genes
variance_filter <- function(omics) {
  omics <- na.omit(omics)
  gene_variance <- apply(omics, 1, var) # Calculate (row-wise) variance
  variable_genes <- omics[gene_variance != 0, ]
  return(variable_genes)
}

ccle_rnaseq_pc_variables <- variance_filter(ccle_rnaseq_pc)
ccle_rppa_variables <- variance_filter(ccle_rppa)
ccle_ms_variables <- variance_filter(ccle_ms)
mclp_rppa_variables <- variance_filter(mclp_rppa)

dir.create(file.path('model_data', 'omics_data', 'selected_features'), recursive = TRUE)
saveRDS(ccle_rnaseq_pc_variables, file = 'model_data/omics_data/selected_features/CCLE RNAseq.rds')
saveRDS(ccle_rppa_variables, file = 'model_data/omics_data/selected_features/CCLE RPPA.rds')
saveRDS(ccle_ms_variables, file = 'model_data/omics_data/selected_features/CCLE MS.rds')
saveRDS(mclp_rppa_variables, file = 'model_data/omics_data/selected_features/MCLP RPPA.rds')

  # Function for batch effects detection in omics data
library(tidyverse)
library(caret)
library(ggfortify)
library(viridis)

pca_processor <- function(matrice, model_context, matrice_name) {
  # Rows: Samples | Columns: Genes
  matrix <- t(matrice)
  data <- as.data.frame(matrix)
  
  missing <- colSums(is.na(data)) / nrow(data)
  # Handle NA values
  columns_to_omit <- names(missing[missing > 0.2])  # Omit columns with >20% missing data
  columns_to_impute <- names(missing[missing <= 0.2])  # Impute columns with <=20% missing data
  
  cleaned_data <- data %>%
    dplyr::select(-dplyr::all_of(columns_to_omit))
  
  if (length(columns_to_impute) > 0) {
    imputation <- preProcess(cleaned_data %>%
                               dplyr::select(dplyr::all_of(columns_to_impute)), method = 'medianImpute')
    imputed_data <- predict(imputation, newdata = cleaned_data %>%
                              dplyr::select(dplyr::all_of(columns_to_impute)))
    cleaned_data <- bind_cols(cleaned_data %>% 
                                dplyr::select(-dplyr::all_of(columns_to_impute)), imputed_data)
  }
  
  # Remove zero-variance genes
  non_zero_var_data <- cleaned_data %>%
    dplyr::select(where(~ is.finite(var(.x, na.rm = TRUE)) && var(.x, na.rm = TRUE) != 0))
  
  if (ncol(non_zero_var_data) == 0) {
    stop("No columns left after removing zero variance columns")
  }
  
  pca_result <- prcomp(non_zero_var_data, scale. = TRUE)
  print(summary(pca_result))
  
  pca_scores <- as.data.frame(pca_result$x)  # Principal component scores
  pca_variance <- summary(pca_result)$importance[2, ]  # Proportion of variance explained by PCs
  
  model_context_expanded <- model_context %>%
    separate_rows(cell_line, sep = " ")
  
  pca_scores <- pca_scores %>%
    rownames_to_column(var = "cell_line") %>%
    left_join(model_context_expanded, by = "cell_line")
  
  formatted_title <- gsub("_", " ", matrice_name)
  pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = lineage)) +
    geom_point(size = 3) +
    labs(
      title = formatted_title,
      x = paste("PC1 (", round(pca_variance[1] * 100, 2), "% variance)", sep = ""),
      y = paste("PC2 (", round(pca_variance[2] * 100, 2), "% variance)", sep = ""),
      color = "Lineage"
    ) +
    theme_grey() +
    scale_color_viridis_d(option = "viridis") +
    theme(legend.position = "right")
  
  print(pca_plot)
  
  output_dir <- file.path('model_data', 'pca_plot', 'selected_features')
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, paste0(formatted_title, ".png"))

  ggsave(output_file, plot = pca_plot, width = 8, height = 6)
}

ccle_rnaseq_pc_variables_pca <- pca_processor(ccle_rnaseq_pc_variables, model_context, "CCLE_RNAseq")
ccle_rppa_variables_pca <- pca_processor(ccle_rppa_variables, model_context, "CCLE_RPPA")
ccle_ms_variables_pca <- pca_processor(ccle_ms_variables, model_context, "CCLE_MS")
mclp_rppa_variables_pca <- pca_processor(mclp_rppa_variables, model_context, "MCLP_RPPA")
