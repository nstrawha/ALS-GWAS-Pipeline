#!/usr/bin/env Rscript

# transcriptomic_analysis.R
# Code by Noah Strawhacker
# Sep. 2025
# Script to find differences in gene transcription between case and control


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))


# Main --------------------------------------------------------------------

main <- function() {
  
  num_chromosomes <- 22
  
  # load in transcriptomic data
  all_transmic_data <- fread(here(
    "original_data", 
    "4_matrix", 
    "AnswerALS-651-T-v1-release6_raw-counts.csv"
    ))
  
  # filter out low counts
  row_sums <- rowSums(all_transmic_data[, -c(1, 2)])
  to_filter <- row_sums < 5
  all_transmic_data <- all_transmic_data[-to_filter, ]
  
  # remove unnecessary cols
  all_counts <- all_transmic_data[, -c(1, 2)]
  
  # get als and control cols
  case_cols <- grep("CASE", colnames(all_counts))
  ctrl_cols <- grep("CTRL", colnames(all_counts))
  
  # subset into als and control
  als_data <- all_counts[, -ctrl_cols, with = FALSE]
  control_data <- all_counts[, -case_cols, with = FALSE]
  
  # reform all counts, sorting case and ctrl
  all_counts <- bind_cols(als_data, control_data)
  
  # sample metadata
  sample_info <- data.frame(
    condition = c(rep("ALS", ncol(als_data)), rep("CTRL", ncol(control_data)))
  )
  rownames(sample_info) <- colnames(all_counts)
  
  ncol(all_counts)
  ncol(t(sample_info))
  
  # get deseq data
  dds <- DESeqDataSetFromMatrix(
    countData = all_counts,
    colData = sample_info,
    design = ~ condition
  )
  
  # run deseq
  dds <- DESeq(dds)
  res <- results(dds)
  
  # reformat output
  resdf <- as.data.frame(res)
  resdf$Geneid <- all_transmic_data$Geneid
  resdf <- resdf %>% 
    arrange(padj)
  
  info_to_merge <- data.frame(
    ensembl_gene_id = resdf$Geneid, 
    tmic_l2fc = resdf$log2FoldChange, 
    tmic_p = resdf$padj
  )
  
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Appending transcriptomic data to chr ", chrom_idx)
  
    # read snp data
    snp_data <- fread(here(
      "summary_statistics_categorized", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".categorized.csv")
    ))
    
    # merge pval and log2 fold change into snp data by gene id
    snp_data <- merge(
      snp_data, 
      info_to_merge, 
      by = "ensembl_gene_id", 
      all.x = TRUE, 
      all.y = FALSE
    )
    
    # fill NAs with "none"
    snp_data <- snp_data %>% 
      mutate(tmic_p = replace_na(as.character(tmic_p), "none")) %>% 
      mutate(tmic_l2fc = replace_na(as.character(tmic_l2fc), "none"))
    
    # sort by p
    snp_data <- snp_data %>% 
      arrange(p)
    
    out_floc <- here("summary_statistics_tmic_ps")
    out_faddress <- here(
      out_floc, 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".tmic.ps.csv")
      )
    
    # write output
    dir.create(
      out_floc, 
      recursive = TRUE, 
      showWarnings = FALSE
    )
    
    fwrite(
      snp_data, 
      file = out_faddress,
      col.names = TRUE, 
      quote = TRUE
    )
  }

  invisible(NULL)
}

main()
