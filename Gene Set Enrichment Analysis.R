  # Load omics features metadata
ccle_rnaseq_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RNAseq.csv')
ccle_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE RPPA.csv')
ccle_ms_metadata <- read.csv(file = 'pset_data/omics/metadata/CCLE MS.csv')
mclp_rppa_metadata <- read.csv(file = 'pset_data/omics/metadata/MCLP RPPA.csv')

   # Load limma differentially expressed features
limma_DE_ccle_rnaseq_ids_distribution <- readRDS(file = 'model_data/limma_results/ccle_rnaseq/logFC Distribution.rds')
limma_DE_ccle_rppa_ids_distribution <- readRDS(file = 'model_data/limma_results/ccle_rppa/logFC Distribution.rds')
limma_DE_mclp_rppa_ids_distribution <- readRDS(file = 'model_data/limma_results/mclp_rppa/logFC Distribution.rds')

  # Load logistic regression significant features
logreg_SAI_ccle_rnaseq_ids_distribution <- readRDS(file = 'model_data/logreg_analysis/ccle_rnaseq/MPP OR Distribution.rds')
logreg_SAI_ccle_rppa_ids_distribution <- readRDS(file = 'model_data/logreg_analysis/ccle_rppa/MPP OR Distribution.rds')
logreg_SAI_ccle_ms_ids_distribution <- readRDS(file = 'model_data/logreg_analysis/ccle_ms/MPP OR Distribution.rds')
logreg_SAI_mclp_rppa_ids_distribution <- readRDS(file = 'model_data/logreg_analysis/mclp_rppa/MPP OR Distribution.rds')

  # Load logistic regression significant features by quadrant
logreg_SAI_ccle_rnaseq_ids_quad <- readRDS(file = 'model_data/logreg_analysis/ccle_rnaseq/MPP OR Quadrant.rds')
logreg_SAI_ccle_rppa_ids_quad <- readRDS(file = 'model_data/logreg_analysis/ccle_rppa/MPP OR Quadrant.rds')
logreg_SAI_ccle_ms_ids_quad <- readRDS(file = 'model_data/logreg_analysis/ccle_ms/MPP OR Quadrant.rds')
logreg_SAI_mclp_rppa_ids_quad <- readRDS(file = 'model_data/logreg_analysis/mclp_rppa/MPP OR Quadrant.rds')

  # Function to retrieve ENTREZ ID based on ENSEMBL Gene ID and HGNC Symbols
library(biomaRt)
library(tidyr)

# Function to map Ensembl gene IDs with fallback to HGNC symbols to Entrez IDs
entrez_id_retriever <- function(metadata) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Step 1: Handle multiple Ensembl Gene IDs
  metadata$ensembl_gene <- strsplit(as.character(metadata$ensembl_gene), " ")
  
  metadata_long <- metadata[rep(seq_len(nrow(metadata)), sapply(metadata$ensembl_gene, length)), ]
  metadata_long$ensembl_gene <- unlist(metadata$ensembl_gene)
  
  # Step 2: Map through Ensembl Gene ID
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = unique(metadata_long$ensembl_gene),  # Map only unique Ensembl IDs
    mart = ensembl
  )
  
  metadata_long <- merge(metadata_long, gene_mapping, 
                         by.x = "ensembl_gene", by.y = "ensembl_gene_id", 
                         all.x = TRUE)
  
  names(metadata_long)[names(metadata_long) == "entrezgene_id"] <- "entrez_id"
  
  # Identify unmapped Ensembl genes
  unmapped_ensembl <- metadata_long$ensembl_gene[is.na(metadata_long$entrez_id)]
  
  if (length(unmapped_ensembl) > 0) {
    # Step 3: Handle multiple HGNC symbols
    metadata_long$gene <- strsplit(as.character(metadata_long$gene), " ")
    metadata_long <- metadata_long[rep(seq_len(nrow(metadata_long)), sapply(metadata_long$gene, length)), ]
    metadata_long$gene <- unlist(metadata_long$gene)
    
    # Step 4: Map using HGNC symbols for unmapped entries
    gene_mapping_hgnc <- getBM(
      attributes = c("hgnc_symbol", "entrezgene_id"),
      filters = "hgnc_symbol",
      values = unique(metadata_long$gene[is.na(metadata_long$entrez_id)]),  # Only for unmapped entries
      mart = ensembl
    )
    
    if (nrow(gene_mapping_hgnc) > 0) {
      metadata_long <- merge(metadata_long, gene_mapping_hgnc, 
                             by.x = "gene", by.y = "hgnc_symbol", 
                             all.x = TRUE, suffixes = c("", "_hgnc"))
      
      if ("entrezgene_id_hgnc" %in% colnames(metadata_long)) {
        metadata_long$entrez_id[is.na(metadata_long$entrez_id)] <- metadata_long$entrezgene_id_hgnc
        metadata_long$entrezgene_id_hgnc <- NULL
      }
    } else {
      cat("No HGNC symbols mapped to Entrez IDs.\n")
    }
  }
  
  # Step 5: Concatenate multiple Entrez IDs
  metadata_final <- aggregate(entrez_id ~ assay_id, data = metadata_long, 
                              FUN = function(x) paste(unique(na.omit(x)), collapse = " "))
  
  metadata <- merge(metadata, metadata_final, by = "assay_id", all.x = TRUE)
  
  # Calculate % that failed to map via both Ensembl Gene ID and HGNC symbol
  total_genes <- nrow(metadata)
  unmapped_final <- sum(metadata$entrez_id == "" | is.na(metadata$entrez_id))
  percent_unmapped <- (unmapped_final / total_genes) * 100
  
  cat(sprintf("Failed to map to Entrez IDs: %.2f%%\n", percent_unmapped))
  
  return(metadata)
}

ccle_rnaseq_metadata <- entrez_id_retriever(ccle_rnaseq_metadata) # Failed to map to Entrez IDs: 56.59%
ccle_rppa_metadata <- entrez_id_retriever(ccle_rppa_metadata) # Failed to map to Entrez IDs: 0.00%
ccle_ms_metadata <- entrez_id_retriever(ccle_ms_metadata) # Failed to map to Entrez IDs: 0.49%
mclp_rppa_metadata <- entrez_id_retriever(mclp_rppa_metadata) # Failed to map to Entrez IDs: 0.00%

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

enrichment_analyzer <- function(model_data_list, metadata, model, omics_source, omics, analysis_type) {
  protein_coding_metadata <- metadata[metadata$gene_type == "protein_coding", ]
  
  dir.create("model_data/gsea_plot", showWarnings = FALSE)

  if (model == "limma") {
    lists_of_interest <- list("logFC < 0" = model_data_list$`logFC < 0`, 
                              "logFC > 0" = model_data_list$`logFC > 0`)
  } else if (model == "logreg") {
    if (analysis_type == "distribution") {
      lists_of_interest <- list("odds_ratio > 1" = model_data_list$`odds_ratio > 1`, 
                                "odds_ratio < 1" = model_data_list$`odds_ratio < 1`)
    } else if (analysis_type == "quad") {
      lists_of_interest <- list(
        "odds_ratio > 1 & median_probability > 0.5" = model_data_list$`odds_ratio > 1 & median_probability > 0.5`,
        "odds_ratio > 1 & median_probability < 0.5" = model_data_list$`odds_ratio > 1 & median_probability < 0.5`,
        "odds_ratio < 1 & median_probability > 0.5" = model_data_list$`odds_ratio < 1 & median_probability > 0.5`,
        "odds_ratio < 1 & median_probability < 0.5" = model_data_list$`odds_ratio < 1 & median_probability < 0.5`
      )
    } else {
      stop("Unsupported analysis type for logreg model")
    }
  } else {
    stop("Unsupported model type")
  }
  
  gsea_results <- list()
  
  for (list_name in names(lists_of_interest)) {
    assay_ids <- lists_of_interest[[list_name]]$assay_ids
    
    matched_metadata <- protein_coding_metadata[protein_coding_metadata$assay_id %in% assay_ids, ]
    
    entrez_ids <- unlist(strsplit(as.character(matched_metadata$entrez_id), " "))
    entrez_ids <- na.omit(entrez_ids)  
    
    if (length(entrez_ids) == 0) {
      cat(paste("No matching protein-coding entrez IDs found for", list_name, "in", model, "\n"))
      next
    }
    
    kegg_result <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05,
                              pAdjustMethod = "BH", qvalueCutoff = 0.05)

    reactome_result <- enrichPathway(gene = entrez_ids, organism = 'human', pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH", qvalueCutoff = 0.05)
    
    if (!is.null(kegg_result) && nrow(as.data.frame(kegg_result)) > 0) {
      plot_title <- paste("KEGG:", model, "of", omics_source, omics, "(condition:", list_name, ")")
      kegg_plot <- dotplot(kegg_result, showCategory = 10) + 
        ggtitle(plot_title)
      
      kegg_plot_filename <- paste0("model_data/gsea_plot/", gsub("[:()]", "", plot_title), ".png")
      ggsave(kegg_plot_filename, plot = kegg_plot, width = 10, height = 6)
      
      print(kegg_plot)
    } else {
      cat(paste("No significant KEGG enrichment found for", list_name, "\n"))
    }
    
    if (!is.null(reactome_result) && nrow(as.data.frame(reactome_result)) > 0) {
      plot_title <- paste("Reactome:", model, "of", omics_source, omics, "(condition:", list_name, ")")
      reactome_plot <- dotplot(reactome_result, showCategory = 10) + 
        ggtitle(plot_title)
      
      reactome_plot_filename <- paste0("model_data/gsea_plot/", gsub("[:()]", "", plot_title), ".png")
      ggsave(reactome_plot_filename, plot = reactome_plot, width = 10, height = 6)
      
      print(reactome_plot)
    } else {
      cat(paste("No significant Reactome enrichment found for", list_name, "\n"))
    }
    
    gsea_results[[list_name]] <- list(
      kegg_result = kegg_result,
      reactome_result = reactome_result
    )
  }
  
  return(gsea_results)
}

limma_DE_ccle_rnaseq_gsea <- enrichment_analyzer(limma_DE_ccle_rnaseq_ids_distribution, ccle_rnaseq_metadata, model = "limma", omics_source = "ccle", omics = "rnaseq", analysis_type = "distribution")
limma_DE_ccle_rppa_gsea <- enrichment_analyzer(limma_DE_ccle_rppa_ids_distribution, ccle_rppa_metadata, model = "limma", omics_source = "ccle", omics = "rppa", analysis_type = "distribution")
limma_DE_mclp_rppa_gsea <- enrichment_analyzer(limma_DE_mclp_rppa_ids_distribution, mclp_rppa_metadata, model = "limma", omics_source = "mclp", omics = "rppa", analysis_type = "distribution")

logreg_SAI_ccle_rnaseq_gsea <- enrichment_analyzer(logreg_SAI_ccle_rnaseq_ids_distribution, ccle_rnaseq_metadata, model = "logreg", omics_source = "ccle", omics = "rnaseq", analysis_type = "distribution")
logreg_SAI_ccle_rppa_gsea <- enrichment_analyzer(logreg_SAI_ccle_rppa_ids_distribution, ccle_rppa_metadata, model = "logreg", omics_source = "ccle", omics = "rppa", analysis_type = "distribution")
logreg_SAI_ccle_ms_gsea <- enrichment_analyzer(logreg_SAI_ccle_ms_ids_distribution, ccle_ms_metadata, model = "logreg", omics_source = "ccle", omics = "ms", analysis_type = "distribution")
logreg_SAI_mclp_rppa_gsea <- enrichment_analyzer(logreg_SAI_mclp_rppa_ids_distribution, mclp_rppa_metadata, model = "logreg", omics_source = "mclp", omics = "rppa", analysis_type = "distribution")

logreg_quad_ccle_rnaseq_gsea <- enrichment_analyzer(logreg_SAI_ccle_rnaseq_ids_quad, ccle_rnaseq_metadata, model = "logreg", omics_source = "ccle", omics = "rnaseq", analysis_type = "quad")
logreg_quad_ccle_rppa_gsea <- enrichment_analyzer(logreg_SAI_ccle_rppa_ids_quad, ccle_rppa_metadata, model = "logreg", omics_source = "ccle", omics = "rppa", analysis_type = "quad")
logreg_quad_ccle_ms_gsea <- enrichment_analyzer(logreg_SAI_ccle_ms_ids_quad, ccle_ms_metadata, model = "logreg", omics_source = "ccle", omics = "ms", analysis_type = "quad")
logreg_quad_mclp_rppa_gsea <- enrichment_analyzer(logreg_SAI_mclp_rppa_ids_quad, mclp_rppa_metadata, model = "logreg", omics_source = "mclp", omics = "rppa", analysis_type = "quad")
