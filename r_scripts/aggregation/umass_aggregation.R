#!/usr/bin/env Rscript

# umass_aggregation.R
# Code by Noah Strawhacker
# Apr. 2025
# Script to handle and aggregate UMass ALS data

rm(list = ls())

# load packages
library(here)
library(data.table)


main <- function() {
  
  # set paths for file handling
  num_chromosomes <- 22
  
  floc_mine <- "summary_statistics_combined"
  file_loc_umass <- here("original_data", "umass_data")
  floc_new  <- "summary_statistics_combined"
  
  fname_umass <- "alsMetaSummaryStats_march21st2018.txt"
  faddress_umass <- here(file_loc_umass, fname_umass)
  
  # read and filter umass data
  data_umass_orig <- fread(faddress_umass)
  cols_to_drop <- c(
    "Direction", 
    "HetISq", 
    "HetChiSq", 
    "HetDf", 
    "HetPVal", 
    "SNPforpos", 
    "effectAlleleMinFreq", 
    "effectAlleleMaxFreq", 
    "effectAlleleFreqStdErr"
    )
  
  data_umass <- data_umass_orig[, setdiff(names(data_umass_orig), cols_to_drop), with = FALSE]
  
  # reformat cols
  colnames(data_umass) <- c("snp", "a1", "a2", "b", "se", "p", "chr", "bp", "freq")
  setcolorder(data_umass, c("chr", "snp", "bp", "a1", "a2", "freq", "b", "se", "p"))
  
  # reformat data to MinE standard
  data_umass$a1 <- toupper(data_umass$a1)
  data_umass$a2 <- toupper(data_umass$a2)
  
  # remove indels
  data_umass[(nchar(data_umass$a1) == 1 & nchar(data_umass$a2) == 1), ]
  
  # add tag to identify data origin
  data_umass$orig <- "umass"
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Aggregating chr ", chrom_idx, " from UMass data")
    
    # open the MinE file and read data
    fname_mine <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.txt")
    faddress_mine <- here(floc_mine, fname_mine)
    data_mine <- fread(faddress_mine)
    data_mine$orig <- "mine" # add tag
    
    # append umass data to mine data
    chrom_data_umass <- data_umass[data_umass$chr == chrom_idx]
    data_to_write <- rbind(data_mine, chrom_data_umass)
    
    # write to output file
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