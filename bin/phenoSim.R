#!/usr/bin/env Rscript
rm(list = ls())
library(ontologyIndex)
library(ontologySimilarity)
args <- commandArgs(trailingOnly = TRUE)

PATIENT_HPO <- args[1]
OMIM_HGMD <- args[2]
OMIM_OBO <- args[3]
OMIM_GENEMAP2 <- args[4]
OMIM_PHENO <- args[5]
OUTFILE_CZ_NAME <- args[6]
OUTFILE_DX_NAME <- args[7]

dat <- read.csv(OMIM_HGMD, sep = "\t")

library(dplyr)
get_HPO_list <- function(df1) {
  df2 <- df1 %>%
    dplyr::group_by(acc_num, phen_id, gene_sym) %>%
    dplyr::summarise(HPO = paste0(hpo_id, collapse = " ")) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  df2$HPO_list <- lapply(df2$HPO, function(i) unlist(strsplit(i, " ")))

  return(df2)
}


# Load HPO_obo
HPO_obo <- get_OBO(OMIM_OBO, propagate_relationships = c("is_a", "part_of"), extract_tags = "minimal")

# set simi_thresh
simi_thresh <- 0


# In public release, there might be empty HGMD phenotype file
if (dim(dat)[1] == 0) {
  col_names <- c("acc_num", "phen_id", "gene_sym", "HPO", "HPO_list", "Similarity_Score")
  dat2 <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
  colnames(dat2) <- col_names
} else {
  dat <- dat[!is.na(dat$hpo_id), ]
  dat2_ori <- get_HPO_list(dat)
  dat2 <- dat2_ori

  # Load patient HPO
  HPO.orig <- read.table(PATIENT_HPO, sep = "\t", fill = T, header = F, stringsAsFactors = F)
  HPO <- HPO.orig$V1

  # remove terms without a HPO ID
  HPO <- HPO[grepl("HP:", HPO)]
  HPO <- list(HPO)

  sim_mat <- get_asym_sim_grid(HPO, dat2$HPO_list, ontology = HPO_obo)
  dat2$Similarity_Score <- as.vector(sim_mat)
  dat2$HPO_list <- unlist(lapply(dat2$HPO_list, function(x) paste0(unlist(x), collapse = "|")))
  dat2 <- dat2[order(dat2$Similarity_Score, decreasing = T), ]
}


write.table(dat2, OUTFILE_CZ_NAME, sep = "\t", quote = F, row.names = F)




# Load OMIM gene-disease data
genemap2 <- readRDS(OMIM_GENEMAP2)

## ---- Load OMIM Phenotype  ----
HPO_orig <- read.table(OMIM_PHENO, sep = "\t", header = T, stringsAsFactors = FALSE, comment.char = "", fill = TRUE, quote = "\"")
OMIM_HPO <- HPO_orig[, c("OMIM_ID", "DiseaseName", "HPO_ID")]
# rename colnames
colnames(OMIM_HPO) <- c("OMIM_ID", "Disease_Name", "HPO_ID")
OMIM_HPO_cl <- unique(OMIM_HPO)


# Function to get_HPO_list
get_HPO_list <- function(df1) {
  df2 <- df1 %>%
    dplyr::group_by(OMIM_ID, Disease_Name) %>%
    dplyr::summarise(HPO = paste0(HPO_ID, collapse = " ")) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  df2$HPO_list <- lapply(df2$HPO, function(i) unlist(strsplit(i, " ")))
  # rownames(df2) <- ifelse(is.na(df2$OMIM_ID), df2$Disease_Name, df2$OMIM_ID)
  return(df2)
}
# Prepare OMIM HPO list for ontology similarity comparision
OMIM_HPO_all <- get_HPO_list(OMIM_HPO_cl[, c("OMIM_ID", "Disease_Name", "HPO_ID")])

HPO.orig <- read.table(PATIENT_HPO, sep = "\t", fill = T, header = F, stringsAsFactors = F)
HPO <- HPO.orig$V1
# remove terms without a HPO ID
HPO <- HPO[grepl("HP:", HPO)]
HPO <- list(HPO)


sim_mat <- get_asym_sim_grid(HPO, OMIM_HPO_all$HPO_list, ontology = HPO_obo)


OMIM_HPO_all$Similarity_Score <- as.vector(sim_mat)

# convert HPO ID to HPO term
OMIM_HPO_all$HPO_term <- unlist(lapply(OMIM_HPO_all$HPO_list, function(x) paste0(HPO_obo$name[unlist(x)], collapse = "|")))

## Add gene - disease relationship
OMIM_HPO_all_wGene <- merge(unique(genemap2[, c("Pheno_ID", "Approved_Gene_Symbol", "Ensembl_Gene_ID", "Entrez_Gene_ID")]),
  OMIM_HPO_all[, c("OMIM_ID", "Disease_Name", "Similarity_Score", "HPO_term")],
  by.y = "OMIM_ID", by.x = "Pheno_ID"
)

colnames(OMIM_HPO_all_wGene)[2] <- "Gene_Symbol"
OMIM_HPO_all_order <- OMIM_HPO_all_wGene[order(OMIM_HPO_all_wGene$Similarity_Score, decreasing = T), ]
# write.table(OMIM_HPO_all_order, paste0(unknown_disease_path,"HPOsimi_all_",output_file_name), sep = "\t", quote = F, row.names = F)

# OMIM_HPO_all_filt <- head(OMIM_HPO_all_order, n = No_candidate)
OMIM_HPO_all_filt <- OMIM_HPO_all_order[OMIM_HPO_all_order$Similarity_Score >= simi_thresh, ]

write.table(OMIM_HPO_all_filt, OUTFILE_DX_NAME, sep = "\t", quote = F, row.names = F)
