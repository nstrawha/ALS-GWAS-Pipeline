#!/usr/bin/env Rscript

# annotate_gof_snps.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to pull out the gain-of-function m6A mutations based on the DRACH consensus motif


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))


# Main --------------------------------------------------------------------

main <- function () {
  
  # set params
  num_chromosomes <- 22
  
  class_types <- c(
    "m6a_gof", 
    "other"
  )
    
  genome_ref <- BSgenome.Hsapiens.UCSC.hg38
  
  floc <- here("summary_statistics_categorized")
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Categorizing GOFs for chr ", chrom_idx)
    
    # read data
    faddress <- here(floc, paste0("als.sumstats.lmm.chr", chrom_idx, ".categorized.csv"))
    fdata <- fread(file = faddress)
    fdata_saved <- fdata[!(fdata$mut_cat == "other"), ]
    fdata <- fdata[fdata$mut_cat == "other", ]
    
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
    
    # categorize snps
    m6a_gof_locs <- grepl("^[AGT][AG][TGC][C][ACT]$", as.character(context_vec)) # match to DR-CH motif
    m6a_gof_locs <- (m6a_gof_locs & (fdata$a2 == "A"))
    m6a_gof_snps <- fdata[m6a_gof_locs, ]
    other_snps <- fdata[!m6a_gof_locs, ]
    
    # add mutation category labels
    m6a_gof_snps$mut_cat <- "m6a_gof"
    other_snps$mut_cat <- "other"
    
    # bind and sort by p-value
    all_snps <- rbindlist(list(fdata_saved, m6a_gof_snps, other_snps))
    all_snps <- all_snps %>% 
      arrange(p)
    
    # write to output file
    out_floc <- here("summary_statistics_categorized")
    out_faddress <- here(
      out_floc, 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".categorized.csv")
      )
  
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