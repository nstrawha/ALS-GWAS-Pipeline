#!/usr/bin/env Rscript

# mutation_classes_by_function.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to plot mutation categories by functional region

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)

main <- function() {
  
  # toggle whether to remove unans from barplot
  remove_unans <- TRUE
  
  # set params
  num_chromosomes <- 22
  
  functions_list <- c(
    "intron", 
    "exon", 
    # "cds", # remove; already present in exons
    "promoter", 
    "terminator", 
    # "utr3", # remove; already present in exons
    # "utr5", # remove; already present in exons
    "unannotated"
  )
  
  class_types <- c(
    "m6a_lof",
    "m5c_lof", 
    "m6a_gof", 
    "other"
  )
  
  # initialize count stores of mutation types by functional region
  mutation_counts_store <- setNames(
    lapply(class_types, function(class) {
      setNames(as.list(numeric(length(functions_list))), functions_list)
    }),
    class_types
  )
  
  func_counts_store <- setNames(
    as.list(numeric(length(functions_list))), 
    functions_list
  )
  
  # iterate through mutation categories
  for (class in class_types) {
    current_class_counts <- mutation_counts_store[[class]]
    
    # iterate through functional regions
    for (func in functions_list) {
      current_func_counts <- current_class_counts[[func]]
      
      floc <- here("summary_statistics_categorized", class, func)
      
      # iterate through chromosomes
      for (chrom_idx in 1:num_chromosomes) {
        
        message("Counting chr ", chrom_idx, " for ", class, " ", func, "s")
        
        # get data categorized by mutation class and functional region
        faddress <- here(floc, paste0("als.sumstats.lmm.chr", chrom_idx, ".", func, ".", class, ".txt"))
        fdata <- fread(file = faddress)
        
        # check if unannotateds are to be counted
        if (func == "unannotated") {
          if (!remove_unans){
            n_snps <- nrow(fdata)
            current_func_counts <- current_func_counts + n_snps
          }
        } else {
          n_snps <- nrow(fdata)
          current_func_counts <- current_func_counts + n_snps
        }

      }
      
      func_counts_store[[func]] <- func_counts_store[[func]] + current_func_counts
      current_class_counts[[func]] <- current_func_counts
    }
    
    mutation_counts_store[[class]] <- current_class_counts
  }
  
  new_class_types <- c(
    "m6A LOF",
    "m5C LOF", 
    "m6A GOF", 
    "All SNPs"
  )
  
  # replace other snps with counts of all snps
  mutation_counts_store[["other"]] <- func_counts_store
  mutation_counts_store <- setNames(
    mutation_counts_store, 
    new_class_types
  )
  
  # get relative frequencies of each function in each mutation class
  class_totals <- lapply(mutation_counts_store, function(cat) sum(unlist(cat)))
  
  rel_freqs_store <- mapply(function(counts, total) {
    lapply(counts, function(x) x / total)
  }, mutation_counts_store, class_totals, SIMPLIFY = FALSE)
  
  # set up data frame for plotting
  df_classes <- rep(new_class_types, each = length(functions_list))
  df_funcs <- rep(functions_list, times = length(new_class_types))
  df_counts <- unlist(mutation_counts_store)
  df_relfreqs <- unlist(rel_freqs_store)
  
  df_to_plot <- data.frame(
    mut_cat = df_classes, 
    func_reg = df_funcs, 
    counts = df_counts, 
    relfreqs = df_relfreqs
  )
  
  snp_cat_barplot <- ggplot(df_to_plot, aes(x = func_reg, y = relfreqs, fill = mut_cat)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
    geom_text(aes(label = counts), 
              position = position_dodge(width = 0.8), 
              vjust = -0.5, size = 2) + 
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Functional Region",
      y = "Relative Frequency",
      fill = "Mutation Category"
    )
  
  # save resulting freq plot
  ggsave(
    filename = here("outputs", "mutation_classes_by_func_barplot.png"), 
    plot = snp_cat_barplot, 
    width = 14, height = 5
  )
  
  invisible(NULL)
}

main()