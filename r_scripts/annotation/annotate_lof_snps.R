#!/usr/bin/env Rscript

# annotate_lof_snps.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to categorize snps based on m6a/m5c status


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

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
  message(typeof(hits))
  candidates <- unique(snp_df[queryHits(hits), ])
  filter(candidates, a1 == allele)
}


# Main --------------------------------------------------------------------

main <- function() {
  
  # set params
  num_chromosomes <- 22
  
  class_types <- c(
    "m6a_lof",
    "m5c_lof", 
    "other"
  )
  
  # get m6a and m5c peaks
  m6a_gr <- load_peaks(here("original_data", "m6as_m5cs", "merged_m6a_peaks.bed"), prefix_to_add = "")
  m5c_gr <- load_peaks(here("original_data", "m6as_m5cs", "merged_m5c_peaks.bed"), prefix_to_add = "chr")
  
  
  # iterate through chroms
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Categorizing LOFs for chr ", chrom_idx)
    
    # unpack uncategorized data
    in_path <- here(
      "summary_statistics_annotated", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".annotated.csv")
      )
    
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
    
    # add labels
    for (type in class_types) {
      current_snps <- all_snps[[type]]
      current_snps$mut_cat <- type
      all_snps[[type]] <- current_snps
    }
    
    all_snps <- rbindlist(all_snps)
      
    # write output files
    out_floc <- here("summary_statistics_categorized")
    out_fname <- paste0("als.sumstats.lmm.chr", chrom_idx, ".categorized.csv")
    out_faddress <- here(out_floc, out_fname)
    
    dir.create(
      out_floc, 
      recursive = TRUE, 
      showWarnings = FALSE
      )
    
    fwrite(
      all_snps,
      file = out_faddress,
      col.names = TRUE, 
      quote = TRUE
      )
    
  }
  
  invisible(NULL)
}

main()