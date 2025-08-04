#!/usr/bin/env Rscript

# proj_mine_aggregation.R
# Code by Noah Strawhacker
# Apr. 2025
# Script to deduplicate and aggregate Project MinE ALS data

rm(list = ls())

# load packages
library(here)
library(data.table)

# function to load data
load_mine_data <- function(year, idx) {
  fname <- paste0("als.sumstats.lmm.chr", idx, ".txt")
  floc <- here("original_data", paste0("summary_statistics_", year))
  faddress <- here(floc, fname, fname) # nested inside folder of same name
  fread(faddress, header = TRUE)
}


main <- function() {
  
  # set paths for file handling
  num_chromosomes <- 22
  
  floc_new  <- "summary_statistics_combined"
  
  # iterate over chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Aggregating chr ", chrom_idx, " from Project MinE data")
    
    # load data
    data_2016 <- load_mine_data(year = 2016, idx = chrom_idx)
    data_2018 <- load_mine_data(year = 2018, idx = chrom_idx)
    
    # remove indels
    data_2016 <- data_2016[(nchar(data_2016$a1) == 1 & nchar(data_2016$a2) == 1), ]
    data_2018 <- data_2018[(nchar(data_2018$a1) == 1 & nchar(data_2018$a2) == 1), ]
   
    # merge data by rsIDs
    merged_data <- merge(
      data_2018, data_2016, 
      by = "snp", 
      suffixes = c("_2018", "_2016"), 
      all.x = TRUE
      )
    
    # record snps for which alleles differ
    diff_data_2016 <- merged_data[(a1_2016 != a1_2018 | a2_2016 != a2_2018)]
    diff_snps_2016 <- diff_data_2016$snp
    diff_data_2016 <- data_2016[snp %in% diff_snps_2016]
    
    # record snps which don't appear in 2018 data
    extra_data_2016 <- data_2016[!snp %in% data_2018$snp]
    
    # append all data
    data_to_write <- rbind(data_2018, diff_data_2016, extra_data_2016)
    
    # write to output file
    colnames(data_to_write) <- colnames(data_2018)
    fname_new <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.txt")
    faddress_new <- here(floc_new, fname_new)
    
    dir.create(
      dirname(faddress_new), 
      recursive = TRUE, 
      showWarnings = FALSE
    )
    
    fwrite(
      data_to_write,
      file = faddress_new,
      sep = "\t",
      col.names = TRUE,
      quote = FALSE
      )
  }

  invisible(NULL)
}

main()