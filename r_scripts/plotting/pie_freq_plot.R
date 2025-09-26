#!/usr/bin/env Rscript

# pie_freq_plot.R
# Code by Noah Strawhacker
# Sep. 2025
# Script to display mutation category frequencies in a pie plot


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# Main --------------------------------------------------------------------

main <- function() {
  
  num_chromosomes <- 22
  
  in_floc <- "summary_statistics_tmic_ps" # TODO change this
  
  plot_store <- list()
  
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Recording SNP frequency data to plot for chr ", chrom_idx)
    
    # read data
    in_faddress <- here(
      in_floc, 
      paste0(paste0("als.sumstats.lmm.chr", chrom_idx, ".tmic.ps.csv"))
    )
    
    snp_data <- fread(file = in_faddress)
    
    # store data to plot
    plot_store[[chrom_idx]] <- snp_data
    
  }
  
  # bind data
  plot_store <- rbindlist(plot_store)
  
  # deduplicate
  plot_store <- plot_store[!duplicated(plot_store$snp), ]
  
  # snp category frequency plot
  df_to_plot <- plot_store %>% 
    group_by(mut_cat) %>% 
    dplyr::count() # dplyr::count is masked by another function
  
  # rearrange results
  df_to_plot <- df_to_plot[match(c(
    "m6a_lof", 
    "m6a_gof", 
    "m5c_lof", 
    "other"
  ), df_to_plot$mut_cat), ]
  
  # create the plot
  labels_to_plot <- c(
    "m5C LOF SNPs", 
    "m6A GOF SNPs", 
    "m6A LOF SNPs", 
    "Other SNPs"
  )
  
  colors_to_plot <- setNames(list(
    "#FF7F78", 
    "#9CE674", 
    "#F59BF5", 
    "#5EAEE1"
  ), df_to_plot$mut_cat)
  
  ggplot(df_to_plot, aes(x = "", y = n, fill = mut_cat)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", clip = "off", direction = -1) + 
    scale_fill_manual(
      name = "SNP Category", 
      values = colors_to_plot, 
      labels = labels_to_plot
    ) +
    labs(title = "SNP Frequency by Mutation Category") + 
    ggplot2::theme(
      panel.spacing = unit(0, "lines")
    ) +
    theme_void()
  
  ggsave(filename = here("outputs", "snp_frequency_pieplot.png"))
  
  invisible(NULL)
}

main()
