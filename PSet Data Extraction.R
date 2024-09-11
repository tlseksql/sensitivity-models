library(SparseArray)
library(SummarizedExperiment)
library(CoreGx)
library(PharmacoGx)
library(pander)

# Treatment response
CTRPv2 <- downloadPSet("CTRPv2_2015")
PRISM <- downloadPSet("PRISM_2020")
GDSC1 <- downloadPSet("GDSC_2020(v1-8.2)")
GDSC2 <- downloadPSet("GDSC_2020(v2-8.2)")

# Functions to uniform-ize data formats for treatment response
library(dplyr)
  # Primary 
drug_response_formatter <- function(dataset) {
  dataset <- dataset %>%
    mutate(cell_line = toupper(gsub("\\[.*?\\]", "", cell_line)),
           cell_line = toupper(gsub("[^A-Za-z0-9]", "", cell_line)),
           treatment = toupper(gsub("[^A-Za-z0-9]", "", treatment))) 
  return(dataset)
}
  # Metadata 
drug_metadata_formatter <- function(dataset) {
  dataset <- dataset %>%
    mutate(treatment = toupper(gsub("[^A-Za-z0-9]", "", treatment)),
           targets = toupper(gsub(";", " , ", targets)))
}

  # Directory
trdirs <- c('data', 'metadata')
sapply(trdirs, function(dir) dir.create(file.path('pset_data', 'treatment_response', dir), recursive = TRUE))

# CTRPv2 data extraction
  # Primary
CTRPv2_cell_line <- CTRPv2@treatmentResponse[["info"]][["sampleid"]]
CTRPv2_treatment <- CTRPv2@treatmentResponse[["info"]][["treatmentid"]]
CTRPv2_auc <- CTRPv2@treatmentResponse[["profiles"]][["aac_published"]]
CTRPv2_aac <- CTRPv2@treatmentResponse[["profiles"]][["aac_recomputed"]]
CTRPv2_ic50 <- CTRPv2@treatmentResponse[["profiles"]][["ic50_recomputed"]]

CTRPv2_drug_response <- data.frame(cell_line = CTRPv2_cell_line,
                                   treatment = CTRPv2_treatment,
                                   auc = CTRPv2_auc,
                                   aac = CTRPv2_aac,
                                   ic50 = CTRPv2_ic50)
CTRPv2_drug_response <- drug_response_formatter(CTRPv2_drug_response)
write.csv(CTRPv2_drug_response, file = "pset_data/treatment_response/data/CTRPv2.csv")

  # Metadata
CTRPv2_compound <- CTRPv2@treatment[["cpd_name"]]
CTRPv2_moa <- CTRPv2@treatment[["target_or_activity_of_compound"]]
CTRPv2_target <- CTRPv2@treatment[["gene_symbol_of_protein_target"]]
  
CTRPv2_drug_metadata <- data.frame(treatment = CTRPv2_compound,
                                   moa = CTRPv2_moa,
                                   targets = CTRPv2_target)
CTRPv2_drug_metadata <- drug_metadata_formatter(CTRPv2_drug_metadata)
write.csv(CTRPv2_drug_metadata, file = "pset_data/treatment_response/metadata/CTRPv2.csv")

# PRISM data extraction
  # Primary
PRISM_cell_line <- PRISM@treatmentResponse[["info"]][["sampleid"]]
PRISM_treatment <- PRISM@treatmentResponse[["info"]][["treatmentid"]]
PRISM_auc <- PRISM@treatmentResponse[["profiles"]][["auc_published"]]
PRISM_aac <- PRISM@treatmentResponse[["profiles"]][["aac_recomputed"]]
PRISM_ic50 <- PRISM@treatmentResponse[["profiles"]][["ic50_recomputed"]]

PRISM_drug_response <- data.frame(cell_line = PRISM_cell_line,
                                   treatment = PRISM_treatment,
                                   auc = PRISM_auc,
                                   aac = PRISM_aac,
                                   ic50 = PRISM_ic50)
PRISM_drug_response <- drug_response_formatter(PRISM_drug_response)
write.csv(PRISM_drug_response, file = "pset_data/treatment_response/data/PRISM.csv")

  # Metadata
PRISM_compound <- PRISM@treatment[["treatmentid"]]
PRISM_moa <- PRISM@treatment[["moa"]]
PRISM_target <- PRISM@treatment[["target"]]

PRISM_drug_metadata <- data.frame(treatment = PRISM_compound,
                                   moa = PRISM_moa,
                                   targets = PRISM_target)
PRISM_drug_metadata <- drug_metadata_formatter(PRISM_drug_metadata)
write.csv(PRISM_drug_metadata, file = "pset_data/treatment_response/metadata/PRISM.csv")

# GDSC1 data extraction
  # Primary
GDSC1_cell_line <- GDSC1@treatmentResponse[["info"]][["sampleid"]]
GDSC1_treatment <- GDSC1@treatmentResponse[["info"]][["treatmentid"]]
GDSC1_aac <- GDSC1@treatmentResponse[["profiles"]][["aac_recomputed"]]
GDSC1_ic50 <- GDSC1@treatmentResponse[["profiles"]][["ic50_recomputed"]]
  
GDSC1_drug_response <- data.frame(cell_line = GDSC1_cell_line,
                                 treatment = GDSC1_treatment,
                                 aac = GDSC1_aac,
                                 ic50 = GDSC1_ic50)
GDSC1_drug_response <- drug_response_formatter(GDSC1_drug_response)
write.csv(GDSC1_drug_response, file = "pset_data/treatment_response/data/GDSC1.csv")
  
  # Metadata
GDSC1_compound <- GDSC1@treatment[["treatmentid"]]
GDSC1_target <- GDSC1@treatment[["TARGET"]]
GDSC1_target_pathway <- GDSC1@treatment[["TARGET_PATHWAY"]]
  
GDSC1_drug_metadata <- data.frame(treatment = GDSC1_compound,
                                 targets = GDSC1_target,
                                 target_pathway = GDSC1_target_pathway)
GDSC1_drug_metadata <- drug_metadata_formatter(GDSC1_drug_metadata)
write.csv(GDSC1_drug_metadata, file = "pset_data/treatment_response/metadata/GDSC1.csv")

# GDSC2 data extraction
  # Primary
GDSC2_cell_line <- GDSC2@treatmentResponse[["info"]][["sampleid"]]
GDSC2_treatment <- GDSC2@treatmentResponse[["info"]][["treatmentid"]]
GDSC2_aac <- GDSC2@treatmentResponse[["profiles"]][["aac_recomputed"]]
GDSC2_ic50 <- GDSC2@treatmentResponse[["profiles"]][["ic50_recomputed"]]

GDSC2_drug_response <- data.frame(cell_line = GDSC2_cell_line,
                                  treatment = GDSC2_treatment,
                                  aac = GDSC2_aac,
                                  ic50 = GDSC2_ic50)
GDSC2_drug_response <- drug_response_formatter(GDSC2_drug_response)
write.csv(GDSC2_drug_response, file = "pset_data/treatment_response/data/GDSC2.csv")

  # Metadata
GDSC2_compound <- GDSC2@treatment[["treatmentid"]]
GDSC2_target <- GDSC2@treatment[["TARGET"]]
GDSC2_target_pathway <- GDSC2@treatment[["TARGET_PATHWAY"]]

GDSC2_drug_metadata <- data.frame(treatment = GDSC2_compound,
                                  targets = GDSC2_target,
                                  target_pathway = GDSC2_target_pathway)
GDSC2_drug_metadata <- drug_metadata_formatter(GDSC2_drug_metadata)
write.csv(GDSC2_drug_metadata, file = "pset_data/treatment_response/metadata/GDSC2.csv")

# Omics
CCLE_RNA <- downloadPSet("CCLE_2015")
CCLE_RPPA <- readRDS(file = 'pset_data/bhk_resource/CCLE_RPPA_SE.rds')
MCLP_RPPA <- readRDS('pset_data/bhk_resource/MCLP_RPPA_SE.rds')
CCLE_MS <- readRDS('pset_data/bhk_resource/CCLE_MS_SE.rds')

# Functions to uniform-ize data formats for omics
  # Primary
  # Colname: sample name | Rowname: assay id
omics_formatter <- function(matrice) {
  colnames(matrice) <- toupper(gsub("[[:punct:] ]", "", colnames(matrice)))
  return(matrice)
}
  # Metadata
  # ID format: ENS(species)(object type)(identifier).(version)
omics_metadata_formatter <- function(dataset) {
  dataset <- dataset %>%
    mutate(across(contains("ensembl"), ~ gsub("\\.[0-9]+", "", .))) %>%  # Remove version number
    mutate(across(everything(), ~ gsub("\\||///", " ", .)))  # Replace "|" or "///" with a space
  return(dataset)
}

  # Directory
omdirs <- c('data', 'metadata')
sapply(omdirs, function(dir) dir.create(file.path('pset_data', 'omics', dir), recursive = TRUE))

# Function to populate missing IDs
library(biomaRt)
ensembl_retriever <- function(dataset) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  if (!("ensembl_gene" %in% names(dataset))) dataset$ensembl_gene <- NA
  if (!("ensembl_transcript" %in% names(dataset))) dataset$ensembl_transcript <- NA
  if (!("ensembl_protein" %in% names(dataset))) dataset$ensembl_protein <- NA
  if (!("gene_type" %in% names(dataset))) dataset$gene_type <- NA
  if (!("transcript" %in% names(dataset))) dataset$transcript <- NA
  
  get_ensembl_data <- function(genes) {
    results <- getBM(
      attributes = c(
        'hgnc_symbol', 
        'ensembl_gene_id', 
        'gene_biotype', 
        'ensembl_transcript_id', 
        'ensembl_peptide_id',
        'external_transcript_name'
      ),
      filters = 'hgnc_symbol',
      values = genes,
      mart = ensembl
    )
    return(results)
  }
  
  for (i in 1:nrow(dataset)) {
    if (any(trimws(dataset[i, c("ensembl_gene", "ensembl_transcript", "ensembl_protein")]) == "")) {
      genes <- unlist(strsplit(as.character(dataset$gene[i]), " "))
      
      if (length(genes) > 0 && all(genes != "")) {
        ensembl_data <- get_ensembl_data(genes)
        
        if (nrow(ensembl_data) > 0) {
          print(paste("Row", i, "HGNC symbols:", paste(genes, collapse = ", ")))  # Debugging output

          ensembl_genes <- unique(ensembl_data$ensembl_gene_id)
          if (!all(is.na(ensembl_genes))) {
            dataset$ensembl_gene[i] <- paste(ensembl_genes, collapse = " ")
          }
          ensembl_transcripts <- unique(ensembl_data$ensembl_transcript_id)
          if (!all(is.na(ensembl_transcripts))) {
            dataset$ensembl_transcript[i] <- paste(ensembl_transcripts, collapse = " ")
          }
          ensembl_proteins <- unique(ensembl_data$ensembl_peptide_id)
          if (!all(is.na(ensembl_proteins))) {
            dataset$ensembl_protein[i] <- paste(ensembl_proteins, collapse = " ")
          }
          gene_types <- unique(ensembl_data$gene_biotype)
          if (!all(is.na(gene_types))) {
            dataset$gene_type[i] <- paste(gene_types, collapse = " ")
          }
          transcripts <- unique(ensembl_data$external_transcript_name)
          if (!all(is.na(transcripts))) {
            dataset$transcript[i] <- paste(transcripts, collapse = " ")
          }
        }
      } else {
        print(paste("Row", i, "has no valid HGNC symbol. Skipping."))
      }
    }
  }
  return(dataset)
}

# CCLE RNAseq data extraction
  # CCLE RNAseq experiment information (run)
  # Source: Expression Atlas (https://www.ebi.ac.uk/gxa/experiments/E-MTAB-2770/Experiment%20Design)
library(readr)
CCLE_RNA_run <- read_tsv("CCLE-RNAseq-experiment-design.tsv")

  # Primary
CCLE_RNA_data <- CCLE_RNA@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]
CCLE_rnaseq <- assay(CCLE_RNA_data)
    # Colname: Replace run serials with sample names
colnames(CCLE_rnaseq) <- CCLE_RNA_run$`Factor Value[cell line]`[match(colnames(CCLE_rnaseq), CCLE_RNA_run$Run)]
CCLE_rnaseq <- omics_formatter(CCLE_rnaseq)
write.csv(CCLE_rnaseq, file = "pset_data/omics/data/CCLE RNAseq.csv")

  # Metadata
CCLE_rnaseq_assay_id <- CCLE_RNA@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@NAMES
CCLE_rnaseq_gene <- CCLE_RNA@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@elementMetadata@listData[["gene_name"]]
CCLE_rnaseq_ensembl_gene <- CCLE_RNA@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@elementMetadata@listData[["gene_id"]]
CCLE_rnaseq_gene_type <- CCLE_RNA@molecularProfiles[["Kallisto_0.46.1.rnaseq"]]@elementMetadata@listData[["gene_type"]]

CCLE_rnaseq_metadata <- data.frame(assay_id = CCLE_rnaseq_assay_id,
                                 gene = CCLE_rnaseq_gene,
                                 ensembl_gene = CCLE_rnaseq_ensembl_gene,
                                 gene_type = CCLE_rnaseq_gene_type)
CCLE_rnaseq_metadata <- omics_metadata_formatter(CCLE_rnaseq_metadata)

write.csv(CCLE_rnaseq_metadata, file = "pset_data/omics/metadata/CCLE RNAseq.csv")

# CCLE RPPA data extraction
  # Primary
CCLE_rppa <- assay(CCLE_RPPA)
CCLE_rppa <- omics_formatter(CCLE_rppa)
write.csv(CCLE_rppa, file = "pset_data/omics/data/CCLE RPPA.csv")

  # Metadata
CCLE_rppa_assay_id <- CCLE_RPPA@metadata[["Antibody_Name_Assay"]]
CCLE_rppa_gene <- CCLE_RPPA@metadata[["Target_Genes"]]
CCLE_rppa_ensembl_gene <- CCLE_RPPA@metadata[["gene_id"]]
CCLE_rppa_gene_type <- CCLE_RPPA@metadata[["gene_type"]]
CCLE_rppa_transcript <- CCLE_RPPA@metadata[["transcript_name"]]
CCLE_rppa_ensembl_transcript <- CCLE_RPPA@metadata[["transcript_id"]]
CCLE_rppa_ensembl_protein <- CCLE_RPPA@metadata[["protein_id"]]
CCLE_rppa_uniprot <- CCLE_RPPA@metadata[["uniprot_id"]]

CCLE_rppa_metadata <- data.frame(assay_id = CCLE_rppa_assay_id,
                                 gene = CCLE_rppa_gene,
                                 ensembl_gene = CCLE_rppa_ensembl_gene,
                                 gene_type = CCLE_rppa_gene_type,
                                 transcript = CCLE_rppa_transcript,
                                 ensembl_transcript = CCLE_rppa_ensembl_transcript,
                                 ensembl_protein = CCLE_rppa_ensembl_protein,
                                 uniprot = CCLE_rppa_uniprot)
CCLE_rppa_metadata <- ensembl_retriever(CCLE_rppa_metadata)
CCLE_rppa_metadata <- omics_metadata_formatter(CCLE_rppa_metadata)
    # Fail to map assay IDs: 
      # Mre11_Caution (previous HGNC symbol: MRE11A | approved HGNC symbol: MRE11) 
      # TIGAR (previous HGNC symbol: C12orf5 | approved HGNC symbol: TIGAR)
CCLE_rppa_metadata$gene[CCLE_rppa_metadata$assay_id == "Mre11_Caution"] <- "MRE11"
CCLE_rppa_metadata$gene[CCLE_rppa_metadata$assay_id == "TIGAR"] <- "TIGAR"

CCLE_rppa_metadata <- ensembl_retriever(CCLE_rppa_metadata)
write.csv(CCLE_rppa_metadata, file = "pset_data/omics/metadata/CCLE RPPA.csv")

# CCLE MS data extraction
  # Primary
CCLE_ms <- assay(CCLE_MS)
    # Colname: Replace CCLE code with sample names
colnames(CCLE_ms) <- CCLE_MS@colData@listData[["Cell.Line"]][match(colnames(CCLE_ms), CCLE_MS@colData@rownames)]
CCLE_ms <- omics_formatter(CCLE_ms)
write.csv(CCLE_ms, file = "pset_data/omics/data/CCLE MS.csv")

  # Metadata
CCLE_ms_assay_id <- CCLE_MS@metadata[["Uniprot_Acc"]]
CCLE_ms_gene <- CCLE_MS@metadata[["Gene_Symbol"]]
CCLE_ms_ensembl_gene <- CCLE_MS@metadata[["gene_id"]]
CCLE_ms_gene_type <- CCLE_MS@metadata[["gene_type"]]
CCLE_ms_transcript <- CCLE_MS@metadata[["transcript_name"]]
CCLE_ms_ensembl_transcript <- CCLE_MS@metadata[["transcript_id"]]
CCLE_ms_ensembl_protein <- CCLE_MS@metadata[["protein_id"]]
CCLE_ms_uniprot <- CCLE_MS@metadata[["uniprot_id"]]

CCLE_ms_metadata <- data.frame(assay_id = CCLE_ms_assay_id,
                               gene = CCLE_ms_gene,
                               ensembl_gene = CCLE_ms_ensembl_gene,
                               gene_type = CCLE_ms_gene_type,
                               transcript = CCLE_ms_transcript,
                               ensembl_transcript = CCLE_ms_ensembl_transcript,
                               ensembl_protein = CCLE_ms_ensembl_protein,
                               uniprot = CCLE_ms_uniprot)
CCLE_ms_metadata <- ensembl_retriever(CCLE_ms_metadata)
CCLE_ms_metadata <- omics_metadata_formatter(CCLE_ms_metadata)

    # Missing gene names 
      # Q96FF7 (approved HGNC symbol: MISP3)
      # C9J7I0 (approved HGNC symbol: UMAD1)
      # H7BZ55 (approved HGNC symbol: CROCC2)
      # P01614 (approved HGNC symbol: IGKV2D-40)
      # Q69YL0 (approved HGNC symbol: NCBP2AS2)
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q96FF7"] <- "MISP3"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "C9J7I0"] <- "UMAD1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H7BZ55"] <- "CROCC2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P01614"] <- "IGKV2D-40"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q69YL0"] <- "NCBP2AS2"

    # Missing all information
      # P63135 (approved HGNC symbol: ERVK-7)
      # H3BRJ5 (approved HGNC symbol: TRIM59-IFT80)
      # C9JNC2 (erged into Q2M2H8, approved HGNC symbol: MGAM2)
      # P01620 (merged into P01619, approved HGNC symbol: IGKV3-20)
      # P18135 (merged into P01619, approved HGNC symbol: IGKV3-20)
      # I3L1I5 (merged into Q6ZSZ5, approved HGNC symbol: ARHGEF18)
      # P01781 (merged into P01780, approved HGNC symbol: IGHV3-7)
      # Un-annotatable [n=19] + no ensembl ids [n=1]
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P63135"] <- "ERVK-7"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H3BRJ5"] <- "TRIM59-IFT80"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "C9JNC2"] <- "MGAM2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P01620"] <- "IGKV3-20"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P18135"] <- "IGKV3-20"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "I3L1I5"] <- "ARHGEF18"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P01781"] <- "IGHV3-7"

    # Fail to map assay IDs: 
      # Q6ZTU2-6 (previous HGNC symbol: EP400NL | approved HGNc symbol: EP400P1)
      # F5HFY4 (merged into Q99733, previous HGNC symbol: NAP1L4b | approved HGNC symbol: NAP1L4)
      # E9PJD7 (merged into P0DTL6, previous HGNC symbol: CYHR1 | approved HGNC symbol: ZFTRAF1)
      # Q9Y4K1 (previous HGNC symbol: AIM1 | approved HGNC symbol: CRYBG1)
      # G3V325 (incorrect HGNC symbol: ATP5J2-PTCD1 | approved HGNC symbol: ATP5MF-PTCD1)
      # Q5SNT6 (merged into Q641Q2, previous HGNC symbol: FAM21B | approved HGNC symbol: WASHC2A)
      # F8W7U8 (previous HGNC symbol: MRE11A | approved HGNC symbol: MRE11)
      # G3V599 (previous HGNC symbol: CTAGE5| approved HGNC symbol: MIA2)
      # A6PVN8 (previous HGNC symbol: PPP2R4 | approved HGNC symbol: PTPA)
      # Q8N655 (merged into Q96JN0, previous HGNC symbol: C10orf12 | approved HGNC symbol: LCOR)
      # E9PEI6 (previous HGNC symbol: DPCR1 | approved HGNC symbol: MUCL3)
      # B5MDU6 (previous HGNC symbol: C2orf43 | approved HGNC symbol: LDAH)
      # F8W0P7 (previous HGNC symbol: ATP5B | approved HGNC symbol: ATP5F1B)
      # O94854 (previous HGNC symbol: KIAA0754 | approved HGNC symbol: KIAA0754)
      # Q5LJA0 (merged into P10155, previous HGNC symbol: TROVE2 | approved HGNC symbol: RO60)
      # G5E9R9 (previous HGNC symbol: TROVE2 | approved HGNC symbol: RO60)
      # H7C0W7 (previous HGNC symbol: FAM126A | approved HGNC symbol: HYCC1)
      # H3BMR9 (incorrect HGNC symbol: KARS | approved HGNC symbol: KARS1)
      # Q5XLA6 (previous HGNC symbol: CARD17 | approved HGNC symbol: CARD17P)
      # A4D1W8 (alias symbol: UCC1 | approved HGNC symbol: EPDR1)
      # Q6ZW33 (merged into O94851, previous HGNC symbol: MICALCL| approved HGNC symbol: MICAL2)
      # E7ET48 (previous HGNC symbol: AIM1L | approved HGNC symbol: CRYBG2)
      # J3KQD0 (merged into Q8TBF2, previous HGNC symbol: FAM213B | approved HGNC symbol: PRXL2B)
      # P20848 (previous HGNC symbol: SERPINA2P | approved HGNC symbol: SERPINA2)
      # Q9BZK3 (previous HGNC symbol: NACAP1 | approved HGNC symbol: NACA4P)
      # B4DEC1 (merged into Q96B96, previous HGNC symbol: TMEM159 | approved HGNC symbol: LDAF1)
      # I3L4Q0 (previous HGNC symbol: FAM195B | approved HGNC symbol: MCRIP1)
      # P0CW22 (merged into P08708, previous HGNC symbol: RPS17L | approved HGNC symbol: RPS17)
      # Q5T5C7 (previous HGNC symbol: SARS | approved HGNC symbol: SARS1)
      # B7ZLF7 (merged into Q6QNK2, previous HGNC symbol: GPR133 | approved HGNC symbol: ADGRD1)
      # H3BPK4 (previous HGNC symbol: MYLPF | approved HGNC symbol: MYL11)
      # A8MSY1 (previous HGNC symbol: TMEM110-MUSTN1 | approved HGNC symbol: STIMATE-MUSTN1)
      # H7C525 (previous HGNC symbol: CCDC58 | approved HGNC symbol: MIX23)
      # G3XAE9 (previous HGNC symbol: FAM179B | approved HGNC symbol: TOGARAM1)
      # J3KN01 (previous HGNC symbol: MLLT4 | approved HGNC symbol: AFDN)
      # Unsalvagable [n=10] <- Withdrawn records [n=9], no ensembl ids [n=1]
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q6ZTU2-6"] <- "EP400P1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "F5HFY4"] <- "NAP1L4"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "E9PJD7"] <- "ZFTRAF1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q9Y4K1"] <- "CRYBG1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "G3V325"] <- "ATP5MF-PTCD1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q5SNT6"] <- "WASHC2A"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "F8W7U8"] <- "MRE11"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "G3V599"] <- "MIA2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "A6PVN8"] <- "PTPA"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q8N655"] <- "LCOR"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "E9PEI6"] <- "MUCL3"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "B5MDU6"] <- "LDAH"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "F8W0P7"] <- "ATP5F1B"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "O94854"] <- "KIAA0754"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q5LJA0"] <- "RO60"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "G5E9R9"] <- "RO60"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H7C0W7"] <- "HYCC1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H3BMR9"] <- "KARS1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q5XLA6"] <- "CARD17P"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "A4D1W8"] <- "EPDR1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q6ZW33"] <- "MICAL2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "E7ET48"] <- "CRYBG2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "J3KQD0"] <- "PRXL2B"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P20848"] <- "SERPINA2"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q9BZK3"] <- "NACA4P"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "B4DEC1"] <- "LDAF1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "I3L4Q0"] <- "MCRIP1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "P0CW22"] <- "RPS17"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "Q5T5C7"] <- "SARS1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "B7ZLF7"] <- "ADGRD1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H3BPK4"] <- "MYL11"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "A8MSY1"] <- "STIMATE-MUSTN1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "H7C525"] <- "MIX23"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "G3XAE9"] <- "TOGARAM1"
CCLE_ms_metadata$gene[CCLE_ms_metadata$assay_id == "J3KN01"] <- "AFDN"

CCLE_ms_metadata <- ensembl_retriever(CCLE_ms_metadata)
write.csv(CCLE_ms_metadata, file = "pset_data/omics/metadata/CCLE MS.csv")

# MCLP RPPA data extraction
  # Primary
MCLP_rppa <- assay(MCLP_RPPA)
MCLP_rppa <- omics_formatter(MCLP_rppa)
write.csv(MCLP_rppa, file = "pset_data/omics/data/MCLP RPPA.csv")

  # Metadata
MCLP_rppa_assay_id <- MCLP_RPPA@rowRanges@ranges@NAMES
MCLP_rppa_gene <- MCLP_RPPA@rowRanges@elementMetadata@listData[["gene_name"]]
MCLP_rppa_ensembl_gene <- MCLP_RPPA@rowRanges@elementMetadata@listData[["gene_id"]]
MCLP_rppa_gene_type <- MCLP_RPPA@rowRanges@elementMetadata@listData[["gene_type"]]
MCLP_rppa_transcript <- MCLP_RPPA@rowRanges@elementMetadata@listData[["transcript_name"]]
MCLP_rppa_ensembl_transcript <- MCLP_RPPA@rowRanges@elementMetadata@listData[["transcript_id"]]
MCLP_rppa_ensembl_protein <- MCLP_RPPA@rowRanges@elementMetadata@listData[["protein_id"]]
MCLP_rppa_uniprot <- MCLP_RPPA@rowRanges@elementMetadata@listData[["uniprot_id"]]

MCLP_rppa_metadata <- data.frame(assay_id = MCLP_rppa_assay_id,
                                 gene = MCLP_rppa_gene,
                                 ensembl_gene = MCLP_rppa_ensembl_gene,
                                 gene_type = MCLP_rppa_gene_type,
                                 transcript = MCLP_rppa_transcript,
                                 ensembl_transcript = MCLP_rppa_ensembl_transcript,
                                 ensembl_protein = MCLP_rppa_ensembl_protein,
                                 uniprot = MCLP_rppa_uniprot)
MCLP_rppa_metadata <- ensembl_retriever(MCLP_rppa_metadata)
MCLP_rppa_metadata <- omics_metadata_formatter(MCLP_rppa_metadata)

write.csv(MCLP_rppa_metadata, file = "pset_data/omics/metadata/MCLP RPPA.csv")
