#!/usr/bin/env Rscript

# snp_dist_histograms.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to create linear histograms of snp distributions

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)
library("BSgenome.Hsapiens.UCSC.hg38")
source(here("r_scripts", "analysis", "peak_calling_functions.R"))

main <- function() {
  
  message("Setting up data for linear histograms")
  
  # set params
  num_chromosomes <- 22
  
  class_types <- c(
    "all",
    "m6a_lof",
    "m5c_lof", 
    "m6a_gof"
  )
  
  bin_size <- 1e5
  
  # build idx list of chromosome start pos
  genome <- BSgenome.Hsapiens.UCSC.hg38
  chrom_names <- paste0("chr", 1:num_chromosomes)
  chrom_lengths <- seqlengths(genome)[chrom_names]
  chrom_starts <- c(0, cumsum(as.numeric(chrom_lengths[-length(chrom_lengths)])))
  names(chrom_starts) <- 1:num_chromosomes
  
  # read files with bins
  data_store <- list()
  
  for (chrom_idx in 1:num_chromosomes) {
    fdata <- fread(file = here(
      "binned_data", 
      paste0("binned_counts_chr", chrom_idx, ".txt")
      ))
    
    data_store[[chrom_idx]] <- fdata
  }
  
  data_store <- rbindlist(data_store)
  
  # correct for chr starts
  data_store[, `:=`(
    bin_start = as.numeric(bin_start),
    bin_end   = as.numeric(bin_end),
    bin_mid   = as.numeric(bin_mid)
  )]
  
  data_df_als <- data_store[data_store$data_type == "als", ]
  data_df_control <- data_store[data_store$data_type == "control", ]
  
  message("Data read and formatted")
  
  # read significant bins results
  significant_bins <- list()
  append_idx <- 1
  
  for (mut in class_types) {
    for (chrom_idx in 1:num_chromosomes) {
      current_df <- fread(file = here(
        "significant_bins", 
        mut, 
        paste0("significant_peaks_chr", chrom_idx, "_", mut, ".txt")
        ))
      
      significant_bins[[append_idx]] <- current_df
      append_idx <- append_idx + 1
    }
  }
  
  significant_bins <- rbindlist(significant_bins)
  
  # mark which rows are significant once
  setkey(data_df_als, chr, bin_start, mut_class)
  setkey(significant_bins, chr, bin_start, mut_class)
  
  # ensure bin_start numeric for join
  data_df_als[, bin_start := as.numeric(bin_start)]
  significant_bins[, bin_start := as.numeric(bin_start)]
  
  # shift bins
  data_df_als[, bin_start_shifted := bin_start + chrom_starts[as.character(chr)]]
  significant_bins[, bin_start_shifted := bin_start + chrom_starts[as.character(chr)]]
  data_df_control[, bin_start_shifted := bin_start + chrom_starts[as.character(chr)]]
  
  # find significant bins
  data_df_als[, significant := !is.na(
    significant_bins[data_df_als, 
                     on = .(bin_start_shifted, mut_class), 
                     which = TRUE]
  )]
  
  x_breaks <- chrom_starts
  x_labels <- paste0("chr", names(chrom_starts))
  
  for (class in class_types) {
    current_df_als <- data_df_als[mut_class == class, ]
    
    als_snp_chromwide_hists <- ggplot(current_df_als, aes(x = bin_start_shifted, y = count, fill = interaction(mut_class, significant))) +
      geom_col(width = bin_size * 0.95, position = position_identity()) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0, 0))) +
      labs(
        x = "Genomic Position",
        y = "SNP Frequency",
        title = paste("ALS Genome-Wide SNP Distribution for", class)
      ) +
      theme(
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),  # vertical labels
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0.6, "lines")
      ) +
      scale_fill_manual(values = setNames(
        c("black", "red"),
        c(paste0(class, ".FALSE"), paste0(class, ".TRUE"))
      ))
    
    ggsave(
      filename = here("outputs", paste0("als_chr_dists_histogram_", class, ".png")), 
      plot = als_snp_chromwide_hists, 
      width = 16, height = 6
    )
  }
  
  # create plots of control snps
  control_snp_chromwide_hists <- ggplot(data_df_control, aes(x = bin_start_shifted, y = count, fill = mut_class)) +
    geom_col(width = bin_size * 0.95) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0, 0))) +
    facet_grid(rows = vars(mut_class), scales = "free_y") +
    labs(
      x = "Genomic Position",
      y = "SNP Frequency",
      title = "Control Genome-Wide SNP Distribution by Mutation Category"
    ) +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.spacing = unit(0.6, "lines"),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 8)
    ) + 
    scale_fill_manual(values = c(
      "m6a_lof" = "red",
      "m6a_gof" = "darkgreen",
      "m5c_lof" = "blue",
      "all" = "black"
    ))
  
  
  ggsave(
    filename = here("outputs", "control_chr_dists_histogram.png"), 
    plot = control_snp_chromwide_hists, 
    width = 16, height = 6
  )
  
  invisible(NULL)
}

main()