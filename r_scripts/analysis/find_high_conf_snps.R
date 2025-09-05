#!/usr/bin/env Rscript

# find_high_conf_snps.R
# Code by Noah Strawhacker
# Aug. 2025
# Script to call out relevant ALS SNPs within significant peaks
# Also generates genome-wide summaries for phenotype and clinical significance

# load packages
library(here)
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

main <- function() {
  
  # set parameters
  num_chromosomes <- 22
  
  functions_list <- c(
    "intron", 
    "exon", 
    "cds",
    "promoter", 
    "terminator", 
    "utr3", 
    "utr5", 
    "unannotated"
  )
  
  mutation_classes <- c(
    "m6a_lof",
    "m5c_lof", 
    "m6a_gof", 
    "other"
  )
  
  # read snp data
  all_snp_data <- list()
  
  for (class in mutation_classes) {
    snp_funclist <- list()
    for (func in functions_list) {
      snp_chromlist <- list()
      for (chrom_idx in 1:num_chromosomes) {
        faddress <- here(
          "summary_statistics_categorized", 
          class, 
          func, 
          paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt")
        )
        fdata <- fread(file = faddress)
        snp_chromlist[[chrom_idx]] <- fdata
      }
      snp_func_data <- rbindlist(snp_chromlist)
      snp_funclist[[func]] <- snp_func_data
    }
    snp_class_data <- rbindlist(snp_funclist)
    all_snp_data[[class]] <- snp_class_data
  }
  
  # deduplicate SNPs
  all_snp_data <- lapply(all_snp_data, function(df) df[!duplicated(df$snp), ])
  
  # create "all" category
  all_snps_df <- rbindlist(all_snp_data)
  all_snps_df <- all_snps_df[!duplicated(all_snps_df$snp)]
  all_snp_data[["all"]] <- all_snps_df
  all_snp_data[["other"]] <- NULL
  
  mutation_classes <- c("m6a_lof", "m5c_lof", "m6a_gof", "all")
  message("SNP data read")
  
  # read significant bin data
  all_bin_data <- list()
  
  for (class in mutation_classes) {
    bin_chromlist <- list()
    for (chrom_idx in 1:num_chromosomes) {
      faddress <- here(
        "significant_bins", 
        class, 
        paste0("significant_peaks_chr", chrom_idx, "_", class, ".txt")
      )
      fdata <- fread(file = faddress)
      bin_chromlist[[chrom_idx]] <- fdata
    }
    chrom_bins <- rbindlist(bin_chromlist)
    all_bin_data[[class]] <- chrom_bins
  }
  
  message("Bin data read")
  
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
    
    # write info/no-info SNPs
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
      quote = FALSE
    )
    
    # write gene p-values
    dir.create(here("significant_snps", "genes_summary"), recursive = TRUE, showWarnings = FALSE)
    
    gene_pvals <- all_sig_snps_df %>%
      group_by(hgnc_symbol) %>%
      summarise(
        avg_p = mean(p, na.rm = TRUE),
        snp_count = n()
      ) %>%
      arrange(avg_p)
    
    # write to file
    fwrite(
      gene_pvals,
      file = here("significant_snps", "genes_summary", paste0("gene_pvals_", class, ".csv")),
      col.names = TRUE,
      quote = FALSE
    )
  }
  
  message("Data written and genome-wide summaries generated")
  invisible(NULL)
}

main()