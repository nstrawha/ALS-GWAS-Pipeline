#!/usr/bin/env Rscript

# lmm_to_bed_convert.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to reformat the combined data for the ALS GWAS into BED files


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
library(here)
library(data.table)


# Main --------------------------------------------------------------------

main <- function() {

  # disable scientific notation
  options(scipen = 999)
  
  # set paths for file handling
  num_chromosomes = 22;
  
  in_floc <- "summary_statistics_combined"
  out_floc <- here("summary_statistics_combined_bed_format")
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Converting chr ", chrom_idx, " data from LMM to BED")
    
    # access old data
    in_fname <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.csv")
    in_faddress <- here(in_floc, in_fname)
    in_data <- fread(in_faddress)
    
    # reformat to BED
    cols_to_drop <- c("a1", "a2", "freq", "b", "se", "p", "orig")
    bed_data <- in_data[, setdiff(names(in_data), cols_to_drop), with = FALSE] # remove irrelevant cols
    bed_data$bpzb <- (bed_data$bp - 1) # create zero-based bp coord col
    bed_data$chr <- paste0("chr", bed_data$chr)
    setcolorder(bed_data, c("chr", "bpzb", "bp", "snp"))
    
    # write to output file
    fname_out <- paste0("als.sumstats.chr", chrom_idx, ".combined.bed")
    faddress_out <- here(out_floc, fname_out)
    
    dir.create(
      out_floc, 
      recursive = TRUE, 
      showWarnings = FALSE
    )
    
    fwrite(
      bed_data, 
      file = faddress_out, 
      sep = "\t",
      col.names = FALSE,
      quote = FALSE
      )
  }
  
  invisible(NULL)
}

main()