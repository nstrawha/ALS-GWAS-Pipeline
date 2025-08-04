#!/usr/bin/env Rscript

# bed_to_lmm_convert.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to reformat the realigned BED files into LMM files

rm(list = ls())

# load packages
library(here)
library(data.table)

# function to load data
load_bed_data <- function(floc, fname, unmap) {
  faddress <- here(floc, fname)
  bed_data <- fread(faddress)
  
  if (unmap) {
    bed_data$V5 = NULL
  }
  
  colnames(bed_data) <- c("chrom", "start", "end", "snp")
  
  return(bed_data)
}

# function for data writing
write_lmm_data <- function(lmm, bed, mapped_indices, hg_ver, idx, floc) {
  
  # separate data based on mapping
  if (hg_ver == 38) {
    data_to_write <- lmm[mapped_indices, ]
    
  } else if (hg_ver == 19) {
    data_to_write <- lmm[-mapped_indices, , ]
    
  }
  
  # assign col names and new bp
  colnames(data_to_write) <- colnames(lmm)
  data_to_write$bp <- bed$end
  
  # write output
  fname <- paste0("als.sumstats.lmm.chr", idx, ".combined.hg", hg_ver, ".txt")
  faddress <- here(floc, fname)
  
  dir.create(
    dirname(faddress), 
    recursive = TRUE, 
    showWarnings = FALSE
  )
  
  fwrite(
    data_to_write, 
    file = faddress, 
    sep = "\t", 
    col.names = TRUE, 
    quote = FALSE
  )
  
  invisible(NULL)
}


main <- function() {
  
  # set paths for file handling
  num_chromosomes <- 22
  
  floc_old_bed <- "summary_statistics_combined_bed_format_hg38"
  floc_old_lmm <- "summary_statistics_combined"
  floc_hg38    <- "summary_statistics_hg38_orig";
  floc_hg19    <- "summary_statistics_hg19_orig";
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Converting chr ", chrom_idx, " data from BED to LMM")
    
    # access bed data
    fname_bed_hg38 <- paste0("als.sumstats.chr", chrom_idx, ".combined.bed")
    bed_data_hg38  <- load_bed_data(floc = floc_old_bed, fname = fname_bed_hg38, unmap = FALSE)
    
    # access unmap data
    fname_bed_hg19 <- paste0("als.sumstats.chr", chrom_idx, ".combined.bed.unmap")
    bed_data_hg19  <- load_bed_data(floc = floc_old_bed, fname = fname_bed_hg19, unmap = TRUE)
    
    # access the old lmm data
    fname_lmm <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.txt")
    faddress_lmm <- here(floc_old_lmm, fname_lmm)
    lmm_data <- fread(faddress_lmm)
    
    # find which snips were successfully mapped
    mapped_snips <- bed_data_hg38$snp
    row_indices_hg38 <- match(mapped_snips, lmm_data$snp)
    
    # organize and write output files
    write_lmm_data(
      lmm = lmm_data, 
      bed = bed_data_hg38, 
      mapped_indices = row_indices_hg38, 
      hg_ver = 38, 
      idx = chrom_idx, 
      floc = floc_hg38
      )
    
    write_lmm_data(
      lmm = lmm_data, 
      bed = bed_data_hg19, 
      mapped_indices = row_indices_hg38, 
      hg_ver = 19, 
      idx = chrom_idx, 
      floc = floc_hg19
    )
  }
  
  invisible(NULL)
}

main()