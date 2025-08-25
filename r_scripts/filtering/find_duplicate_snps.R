#!/usr/bin/env Rscript

# find_duplicate_snps.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to handle duplicate snps in the combined hg38 and hg19 data

rm(list = ls())

# load packages
library(here)
library(data.table)


main <- function() {
  
  # set paths for file handling
  num_chromosomes <- 22
  
  hg_types  <- setNames(list("hg38", "hg19"), c("38", "19"))
  snp_types <- c("duplicate_snps", "unique_snps")
  data_types <- c("hg38.dups", "hg38.uniq", "hg19.dups", "hg19.uniq")
  
  flocs_orig <- paste0("summary_statistics_", hg_types, "_orig")
  
  new_flocs <- list(
    here(paste0("summary_statistics_", hg_types[["38"]]), snp_types), 
    here(paste0("summary_statistics_", hg_types[["19"]]), snp_types)
  )
  
  # create four unique file paths for file categorization
  new_flocs <- setNames(unlist(new_flocs), c("38dups", "38uniq", "19dups", "19uniq"))
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Finding duplicates for ", chrom_idx, " data")
    
    # access original data
    fnames_orig <- paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.", hg_types, ".txt")
    
    faddresses_orig <- setNames(as.list(
      Map(here, flocs_orig, fnames_orig)), 
      hg_types
    )
    
    datas_orig <- setNames(as.list(
      lapply(faddresses_orig, fread)), 
      hg_types
    )
    
    # find duplicate and unique snps
    dup_snps <- lapply(datas_orig, function(df) {
      
      # keep only duplicate snps
      dups <- df[snp %in% df[duplicated(snp) | duplicated(snp, fromLast = TRUE), snp]]
      
      # collapse duplicates by averaging b, se, and p
      dups_avg <- dups[, .(
        bp      = unique(bp)[1],          # bp should be same, but keep unique
        b       = mean(b, na.rm = TRUE),  # average effect sizes
        se      = mean(se, na.rm = TRUE), # same for SE
        p       = mean(p, na.rm = TRUE)   # same for p
      ), by = snp]
      
      return(dups_avg)
    })
    
    dup_bps <- lapply(dup_snps, function(df) df$bp)
    
    uniq_snps <- mapply(
      function(df, dups) df[!(bp %in% dups)],
      datas_orig, dup_bps,
      SIMPLIFY = FALSE
    )
    
    # repackage into list of all snps
    all_snps <- list(
      dup_snps$hg38, 
      uniq_snps$hg38, 
      dup_snps$hg19, 
      uniq_snps$hg19
    )
    
    all_snps <- setNames(all_snps, data_types)
    
    # remove all umass data
    all_snps <- lapply(all_snps, function(snps) snps[orig != "umass"]) # remove umass snps
    all_snps <- lapply(all_snps, function(snps) snps[orig != "answerALS"]) # remove answerALS snps
    all_snps <- lapply(all_snps, function(snps) snps[, orig := NULL])
    all_snps <- setNames(all_snps, data_types )
    
    # write to output files
    new_fnames <- paste0("als.sumstats.lmm.chr", chrom_idx, ".", data_types, ".txt")
    faddresses <- Map(here, new_flocs, new_fnames)
    faddresses <- setNames(faddresses, data_types )
    
    for (type in data_types ) {
      fwrite(
        all_snps[[type]], 
        file = faddresses[[type]], 
        sep = "\t", 
        col.names = TRUE, 
        quote = FALSE
      )
    }
  }
  
  invisible(NULL)
}

main()