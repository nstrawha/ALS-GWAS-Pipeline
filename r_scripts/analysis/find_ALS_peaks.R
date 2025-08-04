#!/usr/bin/env Rscript

# find_ALS_peaks.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to identify ALS peaks using a control distribution
# Note: script HAS to be run in linux - breaks with the fread/cmd line

rm(list = ls())

# load packages
library(data.table)
library(here)
library(tidyverse)
source(here("r_scripts", "analysis", "peak_calling_functions.R"))

main <- function() {
  
  system("which zcat", intern=TRUE)
  system("which awk", intern=TRUE)
  
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
  
  # get col names from small dataset to initialize storage
  dummy <- fread(file = here(
    "summary_statistics_categorized", 
    "m6a_lof", 
    "utr5",
    "als.sumstats.lmm.chr22.utr5.m6a_lof.txt"
  ))
  
  cnames <- colnames(dummy)
  
  als_floc <- here("summary_statistics_categorized")
  
  # gather data (only rows 2 & 3 of vcf, skip meta rows labeled w/ "#")
  message("Reading control data")
  
  control_data_path <- here("original_data", "freq.vcf.gz")
  cmd_to_run <- "zcat /home/sstrawha/ALS_GWAS_Pipeline/original_data/freq.vcf.gz | awk '!/^#/ && NF>=3 { print $2, $3 }'"
  control_data <- fread(
    cmd = cmd_to_run, 
    header = FALSE, 
    fill = TRUE
  )
  
  colnames(control_data) <- c("bp", "snp")
  
  # iterate over chrs, classes, and functions to reconstruct full data
  als_data_list <- list()
  dt_idx <- 1
  
  for (chrom_idx in 1:num_chromosomes) {
    for (class in mutation_classes) {
      for (func in functions_list) {
        
        message("Extracting data from ", class, " ", func, "s")
        
        # get categorized als data
        fdata <- fread(file = here(
          als_floc, 
          class, 
          func, 
          paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt")), 
          verbose = TRUE
        )
        
        # append to record of all als data
        als_data_list[[dt_idx]] <- fdata
        
      }
    }
  }
  
  als_data <- rbindlist(als_data_list, use.names = TRUE, fill = TRUE)
  
  message("All files read")
  
  # deduplicate als data (artifact of double-counting exons and cds, utr, etc)
  als_data <- unique(als_data)
  
  message("Data deduplicated")
  
  # bin data into histograms
  bin_size <- 1e5
  
  binned_control_data <- bin_data(df = control_data, binwidth = bin_size)
  binned_als_data <- bin_data(df = als_data, binwidth = bin_size)
  
  message("Data binned")
  
  # pull significant peaks bin-wise
  num_bootstraps <- 10000
  pval_cutoff <- 0.05
  
  significant_peaks <- call_relative_peaks(
    pop = binned_control_data, 
    obs = binned_als_data, 
    nboots = num_bootstraps, 
    p_cutoff = pval_cutoff
  )
  
  # write results
  fwrite(
    significant_peaks, 
    file = "significant_peaks_als_vs_control.txt", 
    sep = "\t",
    col.names = TRUE, 
    quote = FALSE
  )
  
  invisible(NULL)
}

main()