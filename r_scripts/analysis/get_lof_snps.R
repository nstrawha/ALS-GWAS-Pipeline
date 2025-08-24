#!/usr/bin/env Rscript

# get_mutation_types.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to categorize snps based on m6a/m5c status

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

# function to load peak data as GRanges
load_peaks <- function(path, prefix_to_add) {
  df <- fread(path)
  GRanges(
    seqnames = paste0(prefix_to_add, df$V1), 
    ranges = IRanges(start = df$V2, end = df$V3)
  )
}

# function to categorize snps
categorize_snps <- function(snp_gr, snp_df, peaks_gr, allele) {
  hits <- findOverlaps(snp_gr, peaks_gr, type = "within")
  candidates <- unique(snp_df[queryHits(hits), ])
  filter(candidates, a1 == allele)
}


main <- function() {
  
  # set params
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
  
  class_types <- c(
    "m6a_lof",
    "m5c_lof", 
    "other"
  )
  
  # get m6a and m5c peaks
  m6a_gr <- load_peaks(here("original_data", "m6as_m5cs", "merged_m6a_peaks.bed"), prefix_to_add = "")
  m5c_gr <- load_peaks(here("original_data", "m6as_m5cs", "merged_m5c_peaks.bed"), prefix_to_add = "chr")
  
  # iterate through functional regions
  for (func in functions_list) {
    
    # iterate through chroms
    for (chrom_idx in 3:num_chromosomes) {
      
      message("Categorizing ", func, "s for chr ", chrom_idx, "; LOFs")
      
      # unpack uncategorized data
      in_path <- here("summary_statistics_annotated", func, paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".txt"))
      
      snp_df <- fread(in_path)
      snp_gr <- GRanges(
        seqnames = paste0("chr", snp_df$chr), 
        ranges = IRanges(start = snp_df$bp, end = snp_df$bp)
        )
      
      # helper function for snp handling
      get_snps <- function(peaks_gr, allele) {
        categorize_snps(
          snp_gr = snp_gr,
          snp_df = snp_df, 
          peaks_gr = peaks_gr,
          allele = allele
        )
      }
      
      # categorize snps
      m6a_lof_snps <- get_snps(m6a_gr, "A")
      m5c_lof_snps <- get_snps(m5c_gr, "C")
      
      other_snps <- fsetdiff(
        snp_df, 
        bind_rows(m6a_lof_snps, m5c_lof_snps)
        )
      
      # structure snps
      all_snps <- setNames(
        list(
          m6a_lof_snps, 
          m5c_lof_snps, 
          other_snps
          ), 
        class_types
        )
      
      for (type in class_types) {
        
        # write output files
        out_dir <- here("summary_statistics_categorized", type, func)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        
        fwrite(
          all_snps[[type]],
          file = file.path(
            out_dir, 
            paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", type, ".txt")),
          sep = "\t", 
          quote = FALSE)
      }
      
    }
  }
  
  # get m6a and m5c peak info
  n_peaks_m6a <- length(m6a_gr)
  n_peaks_m5c <- length(m5c_gr)
  avg_peak_len_m6a <- mean(width(m6a_gr))
  avg_peak_len_m5c <- mean(width(m5c_gr))
  total_bps_in_peaks_m6a <- sum(width(m6a_gr))
  total_bps_in_peaks_m5c <- sum(width(m5c_gr))
  
  # format into df and write output
  peak_info <- data.frame(
    n_peaks = c(n_peaks_m6a, n_peaks_m5c),
    avg_peak_len = c(avg_peak_len_m6a, avg_peak_len_m5c), 
    total_bps_in_peaks = c(total_bps_in_peaks_m6a, total_bps_in_peaks_m5c)
    )
  rownames(peak_info) <- c("m6As", "m5Cs")
  
  fwrite(
    peak_info, 
    file = here("outputs", "mutation_peaks_comparison.csv"), 
    col.names = TRUE, 
    row.names = TRUE, 
    quote = FALSE
  )
  
  invisible(NULL)
}

main()