#!/usr/bin/env Rscript

# umass_aggregation.R
# Code by Noah Strawhacker
# Aug. 2025
# Script to handle and aggregate Answer ALS data


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))


# Main --------------------------------------------------------------------

main <- function() {
  
  # set paths for file handling
  num_chromosomes <- 22
  
  floc_orig <- "summary_statistics_combined"
  floc_ans <- here("original_data", "AnswerALS_GWAS")
  floc_out  <- "summary_statistics_combined"
  
  fname_ans <- "VanRheenen_2021_ALL.processed.tsv.gz"
  faddress_ans <- here(floc_ans, fname_ans)
  
  # read and filter answerALS data
  data_ans_orig <- fread(faddress_ans)
  cols_to_drop <- c(
    "MarkerName", 
    "Direction",        
    "HetISq", 
    "HetChiSq",
    "HetDf",
    "HetPVal", 
    "N_effective_EurChiJap", 
    "N_effective"
    )
  
  data_ans <- data_ans_orig[, setdiff(names(data_ans_orig), cols_to_drop), with = FALSE]
  
  # reformat cols
  colnames(data_ans) <- c("a1", "a2", "b", "se", "p", "chr", "bp", "snp")
  data_ans$freq = 0 # dummy freq value
  setcolorder(data_ans, c("chr", "snp", "bp", "a1", "a2", "freq", "b", "se", "p"))
  
  # remove indels
  data_ans <- data_ans[(nchar(data_ans$a1) == 1 & nchar(data_ans$a2) == 1), ]
  
  # remove malformed rows
  data_ans <- data_ans[grepl(":", data_ans$snp), ]
  
  # force alleles to upper case
  data_ans$a1 <- toupper(data_ans$a1)
  data_ans$a2 <- toupper(data_ans$a2)
  
  # add tag to identify data origin
  data_ans$orig <- "answerALS"
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Aggregating chr ", chrom_idx, " from AnswerALS data")
    
    # open the original file and read data
    fname_orig <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.csv")
    faddress_orig <- here(floc_orig, fname_orig)
    data_orig <- fread(faddress_orig)
    
    # append answer als data to mine data
    chrom_data_ans <- data_ans[data_ans$chr == chrom_idx]
    data_to_write <- rbind(data_orig, chrom_data_ans)
    
    # write to output file
    fname_out <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.csv")
    faddress_out <- here(floc_out, fname_out) 
    
    dir.create(
      dirname(faddress_out), 
      recursive = TRUE, 
      showWarnings = FALSE
    )
    
    fwrite(
      data_to_write, 
      file = faddress_out, 
      col.names = TRUE, 
      quote = FALSE
      )
  }
  
  invisible(NULL)
}

main()