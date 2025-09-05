#!/usr/bin/env Rscript

# lmm_to_bed_convert.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to reformat the combined data for the ALS GWAS into BED files

rm(list = ls())

# load packages
library(here)
library(data.table)


main <- function() {

  # disable scientific notation
  options(scipen = 999)
  
  # set paths for file handling
  num_chromosomes = 22;
  
  old_floc <- "summary_statistics_combined"
  new_floc <- "summary_statistics_combined_bed_format"
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Converting chr ", chrom_idx, " data from LMM to BED")
    
    # access old data
    old_fname <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.txt")
    old_faddress <- here(old_floc, old_fname)
    old_data <- fread(old_faddress)
    
    # reformat to BED
    cols_to_drop <- c("a1", "a2", "freq", "b", "se", "p", "orig")
    bed_data <- old_data[, setdiff(names(old_data), cols_to_drop), with = FALSE] # remove irrelevant cols
    bed_data$bpzb <- (bed_data$bp - 1) # create zero-based bp coord col
    bed_data$chr <- paste0("chr", bed_data$chr)
    setcolorder(bed_data, c("chr", "bpzb", "bp", "snp"))
    
    # write to output file
    fname_new <- paste0("als.sumstats.chr", chrom_idx, ".combined.bed")
    faddress_new <- here(new_floc, fname_new)
    
    dir.create(
      dirname(faddress_new), 
      recursive = TRUE, 
      showWarnings = FALSE
    )
    
    fwrite(
      bed_data, 
      file = faddress_new, 
      sep = "\t",
      col.names = FALSE,
      quote = FALSE
      )
  }
  
  invisible(NULL)
}

main()