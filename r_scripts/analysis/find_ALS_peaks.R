#!/usr/bin/env Rscript

# find_ALS_peaks.R
# Code by Noah Strawhacker
# Aug. 2025
# Script to identify ALS peaks using a control distribution
# Note: script HAS to be run in linux - breaks with the fread/cmd line
# cmd to run: Rscript r_scripts/analysis/find_ALS_peaks.R

rm(list = ls())

# load packages
library(data.table)
library(here)
library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
source(here("r_scripts", "analysis", "peak_calling_functions.R"))


main <- function() {
  
  set.seed(123)
  
  system("which zcat", intern = TRUE)
  system("which awk", intern= TRUE)
  
  # set parameters
  num_chromosomes <- 22
  
  genome_ref <- BSgenome.Hsapiens.UCSC.hg38
  
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
  
  mutation_classes <- c(
    "m6a_lof",
    "m5c_lof", 
    "m6a_gof", 
    "other"
  )
  
  # initialize hg38 chr names in the vcf file
  vcf_chr_names <- c(
    "NC_000001.11", 
    "NC_000002.12",
    "NC_000003.12", 
    "NC_000004.12", 
    "NC_000005.10", 
    "NC_000006.12", 
    "NC_000007.14", 
    "NC_000008.11", 
    "NC_000009.12", 
    "NC_000010.11", 
    "NC_000011.10", 
    "NC_000012.12", 
    "NC_000013.11", 
    "NC_000014.9", 
    "NC_000015.10", 
    "NC_000016.10", 
    "NC_000017.11", 
    "NC_000018.10", 
    "NC_000019.10", 
    "NC_000020.11", 
    "NC_000021.9", 
    "NC_000022.11"
  )
  
  # file locations
  als_floc <- here("summary_statistics_categorized")
  vcf_loc <- here("original_data", "freq.vcf.gz")
  
  # bootstrapping params
  bin_size <- 1e5
  num_bootstraps <- 10000
  pval_cutoff <- 0.05
  
  # helper function to get GRanges from peaks bed file
  get_peaks <- function(type) {
    df <- fread(file = here("original_data", "m6as_m5cs", paste0("merged_", type, "_peaks.bed")))
    
    # handle differences in chr naming conventions
    if (type == "m6a") seqnames_to_use = df$V1
    else if (type == "m5c") seqnames_to_use = paste0("chr", df$V1)
    
    GRanges(
      seqnames = seqnames_to_use, 
      ranges = IRanges(start = df$V2, end = df$V3)
    )
  }
  
  m6a_peaks_gr <- get_peaks("m6a")
  m5c_peaks_gr <- get_peaks("m5c")
  
  # helper function to get snps of different mutation classes
  get_class_snps <- function(data_df, data_gr, peaks_gr, allele) {
    hits <- findOverlaps(data_gr, peaks_gr, type = "within")
    data_df <- data_df[queryHits(hits), ]
    data_df <- data_df[data_df$a1 == allele, ]
    return(data_df)
  }
  
  # initialize storage of binned data for later writing
  all_binned_data <- list()
  
  for (chrom_idx in 1:num_chromosomes) {
    
    # iterate over chrs, classes, and functions to reconstruct full data
    als_data_list <- list()
    dt_idx <- 1
    
    # gather control data for the chr
    message("Reading control data for chr ", chrom_idx)
    
    cmd_to_run <- paste0(
      "zcat \"", 
      vcf_loc, 
      "\" | awk -v name=\"",
      vcf_chr_names[[chrom_idx]], 
      "\" '!/^#/ && $1 == name { print $2, $3, $4, $5 }'"
    )
    
    control_data <- fread(
      cmd = cmd_to_run, 
      header = FALSE, 
      fill = TRUE, 
      sep = " "
    )
    
    colnames(control_data) <- c("bp", "snp", "a1", "a2")
    
    for (class in mutation_classes) {
      for (func in functions_list) {
        
        message("Extracting data from ", class, " ", func, "s")
        
        # get categorized als data
        fdata <- fread(file = here(
          als_floc, 
          class, 
          func, 
          paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt"))
        )
        
        # append to record of all als data
        als_data_list[[dt_idx]] <- fdata
        dt_idx <- dt_idx + 1
        
      }
    }
    
    message("Files for chr ", chrom_idx, " read")
    
    als_data <- rbindlist(als_data_list, use.names = TRUE, fill = TRUE)
    
    rm(als_data_list) # clean to reduce overhead
    
    # deduplicate als data (artifact of double-counting exons and cds, utr, etc)
    als_data <- als_data[!duplicated(als_data$snp), ]
    
    # get m6a and m5c lof snps
    als_data_gr <- GRanges(
      seqnames = paste0("chr", chrom_idx),
      ranges = IRanges(start = als_data$bp, end = als_data$bp)
    )
    
    als_m6a_lof_data <- get_class_snps(
      data_df = als_data, 
      data_gr = als_data_gr, 
      peaks_gr = m6a_peaks_gr, 
      allele = "A"
      )
    
    als_m5c_lof_data <- get_class_snps(
      data_df = als_data, 
      data_gr = als_data_gr, 
      peaks_gr = m5c_peaks_gr, 
      allele = "C"
      )
    
    # get m6a gof snps
    chr_name <- paste0("chr", chrom_idx)
    chr_length <- seqlengths(genome_ref)[chr_name]
    valid_idx <- which(als_data$bp >= 1 & als_data$bp <= chr_length)
    als_data_for_gof <- als_data[valid_idx, ]
    
    # create granges obj
    gof_data_gr <- GRanges(
      seqnames = chr_name, 
      ranges = IRanges(start = als_data_for_gof$bp - 2, end = als_data_for_gof$bp + 2)
    )
    
    # find context for all SNPs
    als_context_vec <- as.character(getSeq(genome_ref, gof_data_gr))
    
    als_m6a_gof_locs <- grepl("^[AGT][AG][TGC][C][ACT]$", as.character(als_context_vec)) # match to DR-CH motif
    als_m6a_gof_locs <- (als_m6a_gof_locs & (als_data_for_gof$a2 == "A"))
    als_m6a_gof_data <- als_data_for_gof[als_m6a_gof_locs, ]
    
    # do the same for all control data
    control_data_gr <- GRanges(
      seqnames = paste0("chr", chrom_idx),
      ranges = IRanges(start = control_data$bp, end = control_data$bp)
    )
    
    control_m6a_lof_data <- get_class_snps(
      data_df = control_data, 
      data_gr = control_data_gr, 
      peaks_gr = m6a_peaks_gr, 
      allele = "A"
    )
    
    control_m5c_lof_data <- get_class_snps(
      data_df = control_data, 
      data_gr = control_data_gr, 
      peaks_gr = m5c_peaks_gr, 
      allele = "C"
    )

    # get m6a gof snps
    valid_idx <- which(control_data$bp >= 1 & control_data$bp <= chr_length)
    control_data_for_gof <- control_data[valid_idx, ]
    
    control_gof_data_gr <- GRanges(
      seqnames = chr_name, 
      ranges = IRanges(start = control_data_for_gof$bp - 2, end = control_data_for_gof$bp + 2)
    )
    
    # find context for all SNPs
    control_context_vec <- as.character(getSeq(genome_ref, control_gof_data_gr))
    
    control_m6a_gof_locs <- grepl("^[AGT][AG][TGC][C][ACT]$", as.character(control_context_vec)) # match to DR-CH motif
    control_m6a_gof_locs <- (control_m6a_gof_locs & (control_data_for_gof$a2 == "A"))
    control_m6a_gof_data <- control_data_for_gof[control_m6a_gof_locs, ]
    
    rm(als_data_gr, gof_data_gr, control_data_gr, control_gof_data_gr); # clean to reduce overhead
    
    # we now have:
    # als_data
    # als_m6a_lof_data
    # als_m5c_lof_data
    # als_m6a_gof_data
    # control_data
    # control_m6a_lof_data
    # control_m5c_lof_data
    # control_m6a_gof_data
  
    # package into named list
    output_classes <- c(
      "all", 
      "m6a_lof", 
      "m5c_lof", 
      "m6a_gof"
    )
    
    data_classes <- c("als", "control")
    
    data_to_analyze <- lapply(list(
      list(als_data, control_data), 
      list(als_m6a_lof_data, control_m6a_lof_data), 
      list(als_m5c_lof_data, control_m5c_lof_data), 
      list(als_m6a_gof_data, control_m6a_gof_data)
      ), function(list) {setNames(list, data_classes)})
    
    names(data_to_analyze) <- output_classes
    
    message("Data restructured for binning")
    
    # bin data into histograms
    bin_types <- c("als", "control")
    binned_data <- list()
    
    for (list_idx in 1:length(data_to_analyze)) {
      current_list <- data_to_analyze[[list_idx]]
      current_class <- output_classes[[list_idx]]
     
      for (type_idx in 1:length(bin_types)) {
        current_df <- current_list[[type_idx]]
        current_type <- bin_types[[type_idx]]
        
        current_bins <- bin_data(
          df = current_df, 
          binwidth = bin_size, 
          chr = chrom_idx, 
          class = output_classes[list_idx], 
          type = current_type
        )
        
        current_list[[type_idx]]<- current_bins
      }
      
      binned_data[[list_idx]] <- current_list
      names(binned_data[[list_idx]]) <- bin_types
    }
    
    message("chr ", chrom_idx, " data binned")
    
    rm(data_to_analyze) # clean to reduce overhead
    
    # bootstrap and write
    names(binned_data) <- output_classes
    
    for (class in output_classes) {
      
      # unpack als and control data
      current_dfs <- binned_data[[class]]
      
      current_dfs <- lapply(current_dfs, function(df) {
        df$data_type <- NULL
        df
      })
      
      als_bins <- current_dfs[["als"]]
      control_bins <- current_dfs[["control"]]
      
      # carry out bootstrap analysis
      significant_peaks <- call_relative_peaks(
        pop = control_bins, 
        obs = als_bins, 
        nboots = num_bootstraps, 
        p_cutoff = pval_cutoff, 
        chr = chrom_idx
      )
      
      significant_peaks <- data.table(significant_peaks)
      
      if (nrow(significant_peaks) > 0) {
        
        # add confidence level tags to bins based on p-value
        significant_peaks[, conf_lvl := fifelse(pval <= 5e-8, "5e-8",
                                                fifelse(pval <= 1e-5, "1e-5",
                                                        fifelse(pval <= 1e-3, "1e-3",
                                                                fifelse(pval <= 0.05, "0.05", ">0.05"))))]
      }
      
      sig_rate <- (nrow(significant_peaks) / nrow(control_bins)) * 100
      message("Significance rate for chr", chrom_idx, ", ", class, ": ", sig_rate, " %")
      
      # write results
      output_floc <- here("significant_bins", class)
      
      fwrite(
        significant_peaks, 
        file = here(
          output_floc, 
          paste0("significant_peaks_chr", chrom_idx, "_", class,".txt")
        ), 
        sep = "\t",
        col.names = TRUE, 
        quote = FALSE
      )
    }
    
    # un-nest binned data for writing
    binned_data <- lapply(binned_data, function(list) {rbindlist(list)})
    binned_data <- rbindlist(binned_data)
    
    fwrite(
      binned_data, 
      file = here("binned_data", paste0("binned_counts_chr", chrom_idx, ".txt")), 
      sep = "\t", 
      quote = FALSE, 
      col.names = TRUE
    )
    
    rm(binned_data) # clean to reduce overhead
  }
  
  invisible(NULL)
}

main()