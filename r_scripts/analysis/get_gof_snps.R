#!/usr/bin/env Rscript

# get_gof_snps.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to pull out the gain-of-function m6A mutations based on the DRACH consensus motif

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)


main <- function () {
  
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
    "m6a_gof", 
    "other"
  )
    
  genome_ref <- BSgenome.Hsapiens.UCSC.hg38
  
  # iterate through functional regions
  for (func in functions_list) {
    floc <- here("summary_statistics_categorized", "other", func) # only look through uncategorized snps
    
    # iterate through chromosomes
    for (chrom_idx in 1:num_chromosomes) {
      
      message("Categorizing ", func, "s for chr ", chrom_idx, "; GOFs")
      
      # read data
      faddress <- here(floc, paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".other.txt"))
      fdata <- fread(file = faddress)
      
      # filter for only valid snp positions
      chr_name <- paste0("chr", chrom_idx)
      chr_length <- seqlengths(genome_ref)[chr_name]
      valid_idx <- which(fdata$bp >= 1 & fdata$bp <= chr_length)
      fdata <- fdata[valid_idx, ]
      
      # create granges obj
      data_gr <- GRanges(
        seqnames = chr_name, 
        ranges = IRanges(start = fdata$bp - 2, end = fdata$bp + 2)
      )
      
      # find context for all SNPs
      context_vec <- as.character(getSeq(genome_ref, data_gr))
      
      m6a_gof_locs <- grepl("^[AGT][AG][TGC][C][ACT]$", as.character(context_vec)) # match to DR-CH motif
      m6a_gof_locs <- (m6a_gof_locs & (fdata$a2 == "A"))
      m6a_gof_snps <- fdata[m6a_gof_locs, ]
      other_snps <- fdata[!m6a_gof_locs, ]
      all_snps <- setNames(list(m6a_gof_snps, other_snps), class_types)
      
      # write to output files
      for (class in class_types) {
        current_floc <- here("summary_statistics_categorized", class, func)
        dir.create(current_floc, recursive = TRUE, showWarnings = FALSE)
        current_faddress <- here(current_floc, paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt"))
      
        fwrite(
          all_snps[[class]], 
          file = current_faddress, 
          sep = "\t", 
          col.names = TRUE, 
          quote = FALSE
        )
      }
      
    }
  }
  
  invisible(NULL)
}

main()