#!/usr/bin/env Rscript

# call_high_conf_snps.R
# Code by Noah Strawhacker
# Aug. 2025
# Script to call out relevant ALS SNPs within significant peaks for regular 1e5 bins
# Also generates genome-wide summaries for phenotype and clinical significance


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))


# Main --------------------------------------------------------------------

main <- function() {
  
  # set parameters
  num_chromosomes <- 22
  
  class_types <- c(
    "m6a_lof", 
    "m6a_gof", 
    "m5c_lof", 
    "other"
  )
  

  # Read and format data --------------------------------------------------
  
  # read snp data
  all_snp_data <- list()
  
  for (chrom_idx in 22:num_chromosomes) {
    faddress <- here(
      "summary_statistics_tmic_ps", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".tmic.ps.csv")
    )
    fdata <- fread(file = faddress)
    all_snp_data[[chrom_idx]] <- fdata
  }
  
  all_snp_data <- rbindlist(all_snp_data)
  
  # deduplicate SNPs
  all_snp_data <- all_snp_data[!duplicated(all_snp_data$snp), ]
  
  message("SNP data read")
  
  # read bin data
  all_bin_data <- list()
  
  for (class in class_types) {
    class_sublist <- list()
    
    for (chrom_idx in 22:num_chromosomes) {
      fdata <- fread(
        here("bin_ps", class, paste0("bin_ps_chr", chrom_idx, ".csv"))
        )
      class_sublist[[chrom_idx]] <- fdata
    }
    
    data_to_record <- rbindlist(class_sublist)
    data_to_record <- data_to_record[grepl("gene_id", data_to_record$gene_id), ]
    data_to_record$gene_id <- sub('.*gene_id \\"([^\\"]+)\\".*', '\\1', data_to_record$gene_id)
    
    all_bin_data[[class]] <- data_to_record
  }
  
  all_bin_data <- rbindlist(all_bin_data)
  gene_data_to_merge <- all_bin_data[!duplicated(all_bin_data$gene_id) & all_bin_data$gene_id != "none", ]
  gene_data_to_merge <- data.frame(
    ensembl_gene_id = gene_data_to_merge$gene_id, 
    gene_p = gene_data_to_merge$gene_p
  )
  
  message("Bin data read")
  
  
  # Call out SNPs -------------------------------------------------------
  
  merged_data <- merge(all_snp_data, gene_data_to_merge, all.x = TRUE, by = "ensembl_gene_id")
  merged_data$gene_p[is.na(merged_data$gene_p)] <- "none"
  
  # write output
  for (class in class_types) {
    data_to_write <- merged_data[merged_data$mut_cat == class, ]
    fwrite(
      data_to_write, 
      file = here("outputs", paste0("result_snps_", class, ".csv")), 
      col.names = TRUE,
      quote = TRUE
    )
  }
  
  message("Data written")
  
  invisible(NULL)
}

main()