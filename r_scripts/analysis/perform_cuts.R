#!/usr/bin/env Rscript

# perform_cuts.R
# Code by Noah Strawhacker
# Sep. 2025
# Script perform p-value (and other) cuts for fully annotated SNP data


# Setup and functions -----------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

# short function to clean up readability later on
cut_data <- function(df, col, cutoff, dir) {
  if (dir == "greater") cut_df <- df[df[[col]] > cutoff, ]
  else if (dir == "less") cut_df <- df[df[[col]] < cutoff, ]
  
  return(cut_df)
}


# Main --------------------------------------------------------------------

main <- function() {
  
  num_chromosomes <- 22
  
  # read all data
  all_data <- list()
  
  floc_in <- "summary_statistics_tmic_ps" # TODO change later
  
  for (chrom_idx in 1:num_chromosomes) {
    faddress_in <- here(
      floc_in, 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".tmic.ps.csv")
      )
    
    fdata <- fread(file = faddress_in)
    all_data[[chrom_idx]] <- fdata
  }
  
  all_data_orig <- rbindlist(all_data)
  
  # perform cuts
  snp_data <- cut_data(all_data_orig, col = "p", cutoff = 0.05, dir = "less")
  # snp_data <- cut_data(snp_data, col = "freq", cutoff = 0.90, dir = "greater")
  
  mutation_cats <- c(
    "m6a_lof", 
    "m6a_gof", 
    "m5c_lof", 
    "other"
  )
  
  for (cat in mutation_cats) {
    current_data <- snp_data[mut_cat == cat, ]
    
    if (!(cat == "other")) {
      current_data <- current_data[!(current_data$tmic_p == "none"), ]
      current_data <- cut_data(current_data, col = "tmic_p", cutoff = 0.05, dir = "less")
    }
    
    current_data <- current_data %>% 
      arrange(p)
    
    current_data <- current_data[!duplicated(current_data$snp), ]
    
    fwrite(
      current_data, 
      file = here("outputs", paste0("snps_of_interest_", cat, ".csv")), 
      quote = TRUE, 
      col.names = TRUE
    )
  }
  
  invisible(NULL)
}

main()