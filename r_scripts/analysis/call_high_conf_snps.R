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
  
  for (chrom_idx in 1:num_chromosomes) {
    faddress <- here(
      "summary_statistics_categorized", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".categorized.csv")
    )
    fdata <- fread(file = faddress)
    all_snp_data[[chrom_idx]] <- snp_class_data
  }
  
  all_snp_data <- rbindlist(all_snp_data)
  
  # deduplicate SNPs
  all_snp_data <- all_snp_data[!duplicated(all_snp_data$snp), ]
  
  # construct named list of classes
  snp_data_list <- list()
  for (class in class_types) {
    snp_data_list[[class]] <- all_snp_data[all_snp_data$mut_cat == class, ]
  }
  
  message("SNP data read")
  
  # read significant bin data
  bin_types <- c("") # TODO
  
  for (type in bin_types) {
    
    all_bin_data <- list()
    
    for (class in class_types) {
      bin_chromlist <- list()
      for (chrom_idx in 1:num_chromosomes) {
        faddress <- here(
          "significant_bins", 
          class, 
          paste0(type, "significant_peaks_chr", chrom_idx, "_", class, ".csv")
        )
        fdata <- fread(file = faddress)
        bin_chromlist[[chrom_idx]] <- fdata
      }
      chrom_bins <- rbindlist(bin_chromlist)
      all_bin_data[[class]] <- chrom_bins
    }
    
    message("Bin data read")
    
    
    # Call out SNPs -------------------------------------------------------
    
    # find snps within significant bins and generate genome-wide summaries
    for (class in mutation_classes) {
      current_class_snp_df <- all_snp_data[[class]]
      current_class_bin_df <- all_bin_data[[class]]
      
      # initialize storage
      genome_pheno_counts <- list()
      genome_clinsig_counts <- list()
      
      sig_snps_list <- list()
      
      for (chrom_idx in 1:num_chromosomes) {
        current_chr_snps <- current_class_snp_df[current_class_snp_df$chr == chrom_idx, ]
        current_chr_bins <- current_class_bin_df[current_class_bin_df$chr == chrom_idx, ]
        
        current_snps_gr <- if (nrow(current_chr_snps) > 0) GRanges(
          seqnames = paste0("chr", chrom_idx),
          ranges = IRanges(start = current_chr_snps$bp, end = current_chr_snps$bp)
        ) else GRanges()
        
        current_bins_gr <- if (nrow(current_chr_bins) > 0) GRanges(
          seqnames = paste0("chr", chrom_idx),
          ranges = IRanges(start = current_chr_bins$bin_start, end = current_chr_bins$bin_end)
        ) else GRanges()
        
        hits <- findOverlaps(current_snps_gr, current_bins_gr, type = "within")
        sig_snps <- current_chr_snps[queryHits(hits), ]
        sig_snps <- current_chr_snps # TODO just this
  
        # filter for significant p-values
        sig_snps <- sig_snps[sig_snps$p <= 0.05, ]
        
        # store snps for later writing
        sig_snps_list[[chrom_idx]] <- sig_snps
        
        # store genome-wide counts
        genome_pheno_counts[[chrom_idx]] <- sig_snps[, .N, by = phenotype_description]
        genome_clinsig_counts[[chrom_idx]] <- sig_snps[, .N, by = clinical_significance]
      }
      
      # aggregate all snps
      all_sig_snps_df <- rbindlist(sig_snps_list)
      
      # aggregate genome-wide counts
      genome_pheno_df <- rbindlist(genome_pheno_counts)[, .(count = sum(N)), by = phenotype_description][order(-count)]
      genome_clinsig_df <- rbindlist(genome_clinsig_counts)[, .(count = sum(N)), by = clinical_significance][order(-count)]
      
      # filters
      # all_sig_snps_df <- all_sig_snps_df[all_sig_snps_df$freq > 0.90, ]
      
      # write SNPs
      output_faddress <- here(
        "significant_snps", 
        paste0("significant_snps_", class, ".csv")
      )
      
      dir.create(here("significant_snps"), recursive = TRUE, showWarnings = FALSE)
      
      # sort by pval
      all_sig_snps_df <- all_sig_snps_df %>% 
        arrange(p)
      
      fwrite(
        all_sig_snps_df, 
        file = output_faddress, 
        col.names = TRUE, 
        quote = TRUE
      )
      
      # write gene p-values
      dir.create(here("significant_snps", "genes_summary"), recursive = TRUE, showWarnings = FALSE)
      
      gene_pvals <- all_sig_snps_df %>%
        group_by(ensembl_gene_id) %>%
        summarise(
          min_p = min(p, na.rm = TRUE),
          snp_count = n()
        ) %>%
        arrange(min_p)
      
      # write to file
      fwrite(
        gene_pvals,
        file = here("significant_snps", "genes_summary", paste0("gene_pvals_", class, ".csv")),
        col.names = TRUE,
        quote = FALSE
      )
    }
  }
  
  message("Data written and genome-wide summaries generated")
  
  invisible(NULL)
}

main()