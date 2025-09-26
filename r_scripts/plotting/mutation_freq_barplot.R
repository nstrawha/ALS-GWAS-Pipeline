#!/usr/bin/env Rscript

# mutation_freq_barplot.R
# Code by Noah Strawhacker
# Sep. 2025
# Script to display mutation category frequencies by functional regions in a bar plot


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
  plot_store <- plot_store[!duplicated(plot_store), ]
  
  # create frequency table
  df_to_plot <- as.data.frame(table(plot_store$func, plot_store$mut_cat))
  colnames(df_to_plot) <- c("func", "mut_cat", "freq")
  df_to_plot <- df_to_plot %>% 
    group_by(mut_cat) %>% 
    mutate(relfreq = freq / sum(freq))
  
  # rename functional regions
  df_to_plot$func <- rep(c(
    "CDS", 
    "Exon", 
    "Intron",
    "Promoter", 
    "Terminator", 
    "Unannotated", 
    "3' UTR", 
    "5' UTR"
  ), length(unique(df_to_plot$mut_cat)))
  
  labels_to_plot <- c(
    "m5C LOF SNPs", 
    "m6A GOF SNPs", 
    "m6A LOF SNPs", 
    "Other SNPs"
  )
  
  colors_to_plot <- setNames(list(
    "#F59BF5", 
    "#9CE674", 
    "#FF7F78", 
    "#5EAEE1"
  ), unique(df_to_plot$mut_cat))
  
  # plot
  ggplot(df_to_plot, aes(x = func, y = relfreq, fill = mut_cat)) + 
    geom_bar(
      width = 0.8, 
      stat = "identity", 
      position = position_dodge(width = 0.8), 
      color = "black"
      ) + 
    geom_text(aes(label = freq), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2) + 
    scale_fill_manual(
      name = "SNP Category", 
      values = colors_to_plot, 
      labels = labels_to_plot
    ) + 
    labs(
      title = "Relative SNP Frequency by Mutation Category and Functional Region", 
      x = "Functional Region",
      y = "Relative SNP Frequency Within Mutation Categories",
    ) + 
    theme_minimal()
  
  ggsave(
    filename = here("outputs", "mutation_classes_by_func_barplot.png"), 
    width = 14, height = 5
    )
  
  
  invisible(NULL)
}

main()