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

main <- function() {
  
  message("Setting up data for linear histograms")
  
  # set params
  num_chromosomes <- 22
  
  functions_list <- c(
    "intron", 
    "exon", 
    "cds", # remove; already present in exons
    "promoter", 
    "terminator", 
    "utr3", # remove; already present in exons
    "utr5", # remove; already present in exons
    "unannotated"
  )
  
  class_types <- c(
    "m6a_lof",
    "m5c_lof", 
    "m6a_gof", 
    "other"
  )
  
  # get col names from small dataset
  dummy <- fread(file = here(
    "summary_statistics_categorized", 
    "m6a_lof", 
    "utr5",
    "als.sumstats.lmm.chr22.utr5.m6a_lof.txt"
  ))
  
  cnames <- colnames(dummy)
  
  # initialize named lists and sublists of dfs to store data for rcircos plot
  data_dfs_store <- setNames(
    lapply(class_types, function(class) {
      setNames(
        lapply(functions_list, function(func) {
          df <- data.frame(matrix(ncol = length(cnames), nrow = 0))
          colnames(df) <- cnames
          return(df)
        }),
        functions_list
      )
    }),
    class_types
  )
  
  # iterate through mutation categories
  for (class in class_types) {
    current_class_dfs <- data_dfs_store[[class]]
    
    # iterate through functional regions
    for (func in functions_list) {
      current_func_dfs <- current_class_dfs[[func]]
      
      floc <- here("summary_statistics_categorized", class, func)
      
      # iterate through chromosomes
      for (chrom_idx in 1:num_chromosomes) {
        
        message("Counting chr ", chrom_idx, " for ", class, " ", func, "s")
        
        # get data categorized by mutation class and functional region
        faddress <- here(floc, paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt"))
        fdata <- fread(file = faddress)
        
        # record df
        current_func_dfs <- rbind(current_func_dfs, fdata)
      }
      
      current_class_dfs[[func]] <- current_func_dfs
    }
    
    data_dfs_store[[class]] <- current_class_dfs
  }
  
  new_class_types <- c(
    "m6A LOF",
    "m5C LOF", 
    "m6A GOF", 
    "All SNPs"
  )
  
  # unlist and load data to plot as large dfs for each mutation class
  data_dfs_store <- lapply(data_dfs_store, function(dfs) rbindlist(dfs))
  
  # replace other mutation category with all mutations
  all_snps_df <- rbindlist(data_dfs_store)
  colnames(all_snps_df) <- cnames
  data_dfs_store[["other"]] <- all_snps_df
  data_dfs_store <- setNames(data_dfs_store, new_class_types)

  # remove dups
  data_dfs_store <- lapply(data_dfs_store, function(df) {
    df %>%
      distinct(snp, .keep_all = TRUE)
  })
  
  bin_size <- 1e4
  binned_all_classes <- list()
  
  # bin snps
  for (class in new_class_types) {
    current_mat <- data_dfs_store[[class]]
    
    snp_data <- data.frame(
      chr = current_mat$chr,
      bp = current_mat$bp
    )
    
    snp_data$bin_start <- as.integer(floor(snp_data$bp / bin_size) * bin_size)
    snp_data$bin_end <- as.integer(snp_data$bin_start + bin_size)
    
    total_snps <- nrow(snp_data)
    
    binned_counts <- snp_data %>%
      group_by(chr, bin_start, bin_end) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(
        rel_freq = count / total_snps,
        mut_class = class,
        bin_mid = bin_start + bin_size / 2
      )
    
    binned_all_classes[[class]] <- binned_counts
  }
  
  # combine into single dataframe
  df_binned_plot <- bind_rows(binned_all_classes)
  
  # create plot
  snp_chromwide_hists <- ggplot(df_binned_plot, aes(x = bin_mid, y = count, fill = mut_class)) +
    geom_col(width = bin_size * 0.95) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))) +
    facet_grid(rows = vars(mut_class), scales = "free_y") +
    labs(
      x = "Genomic Position",
      y = "SNP Frequency",
      title = "Genome-Wide SNP Distribution by Mutation Category"
    ) +
    theme(
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.spacing = unit(0.6, "lines"),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8)
    ) + 
    scale_fill_manual(values = c(
      "m6A LOF" = "red",
      "m6A GOF" = "darkgreen",
      "m5C LOF" = "blue",
      "All SNPs" = "black"
    ))
  
  ggsave(
    filename = here("outputs", "chr_dists_histogram.png"), 
    plot = snp_chromwide_hists, 
    width = 16, height = 4
  )
  
  # write data for future analysis
  fwrite(
    df_binned_plot, 
    file = here("binned_data.txt"), 
    sep = "\t", 
    col.names = TRUE, 
    quote = FALSE
  )
  
  invisible(NULL)
}

main()