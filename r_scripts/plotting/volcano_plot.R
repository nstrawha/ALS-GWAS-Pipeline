#!/usr/bin/env Rscript

# volcano_plot.R
# Code by Noah Strawhacker
# Sep. 2025
# Script to display tmic data in a volcano plot

# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(EnhancedVolcano))


# Main --------------------------------------------------------------------

main <- function() {
  
  num_chromosomes <- 22
  plot_p_cutoff <- 0.05 # cut
  floc <- here("summary_statistics_tmic_ps")
  
  snps_store <- list()
  
  # read data
  for (chrom_idx in 1:num_chromosomes) {
    
    faddress <- here(
      floc, 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".tmic.ps.csv")
    )
    
    snps_store[[chrom_idx]] <- fread(faddress)
  }
  
  snps_store <- rbindlist(snps_store)
  
  # separate by mutation class
  mut_classes <- c(
    "m6a_lof", 
    "m6a_gof", 
    "m5c_lof", 
    "other"
  )
  
  mut_names <- c(
    "m6A LoF", 
    "m6A GoF", 
    "m5C LoF",
    "Other"
  )
  
  names(mut_names) <- mut_classes
  
  classes_store <- list()
  
  for (class in mut_classes) {
    
    # only pull snps with gene symbols and tmic data
    snps_to_store <- filter(snps_store, mut_cat == class)
    
    snps_to_store <- snps_to_store %>%
      filter(hgnc_symbol != "none", tmic_l2fc != "none", tmic_p != "none") %>%
      mutate(tmic_l2fc = as.numeric(tmic_l2fc))
    
    snps_to_store <- snps_to_store %>%
      distinct(hgnc_symbol, .keep_all = TRUE) %>% 
      filter(p < plot_p_cutoff)
    
    data_to_store <- data.frame(
      Gene = snps_to_store$hgnc_symbol, 
      log2FC = as.numeric(snps_to_store$tmic_l2fc), 
      p_value = as.numeric(snps_to_store$tmic_p)
      )
    
    classes_store[[class]] <- data_to_store
  }
  
  # create plots
  for (class in mut_classes) {
    current_data <- classes_store[[class]]
    name <- mut_names[[class]]
      
    plt <- EnhancedVolcano(current_data,
                    lab = current_data$Gene,
                    x = "log2FC",
                    y = "p_value",
                    title = paste("Volcano Plot of Differential Expression for", name),
                    pCutoff = 0.05,
                    FCcutoff = 0.75,
                    pointSize = 1.0,
                    labSize = 3.0, 
                    xlim = c(-2, 2)
                    )
    
    ggsave(
      file = here("outputs", paste0("volcano_plot_", class, ".png")), 
      plot = plt, 
      width = 10, 
      height = 8
      )
  }
  
  
  invisible(NULL)
}

main()