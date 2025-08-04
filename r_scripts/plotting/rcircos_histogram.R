#!/usr/bin/env Rscript

# rcircos_histogram.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to create rcircos circular histogram of snp distributions

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)
library(RCircos)
source(here("r_scripts", "plotting", "rcircos_modified_functions.R"))

main <- function() {
  
  message("Setting up data for RCircos histogram")
  
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
  
  # set up rcircos
  pdf(here("outputs", "rcircos_plot.pdf"), width = 7, height = 7)
  
  # load data for rcircos plot
  data(UCSC.HG38.Human.CytoBandIdeogram)
  
  # load ideogram
  cyto <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(
    cyto.info = cyto, 
    chr.exclude = c("chrX", "chrY"),
    tracks.inside = 0,
    tracks.outside = 4
  )
  
  # reset plot params
  params <- RCircos.Get.Plot.Parameters()
  params$chr.ideo.pos <- 0.7          # shrink ideogram
  params$track.out.start <- 0.9       # shrink hists around ideogram
  params$track.background <- "white"  # transparent hist bgs
  params$grid.line.color <- "white"   # transparent hist bgs
  params$track.height <- 0.35         # taller histograms
  RCircos.Reset.Plot.Parameters.Custom(new.params = params)
  
  # reset ideogram params
  ideo <- RCircos.Get.Plot.Ideogram()
  
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot.Custom()
  
  text(0, 0.05, "ALS SNP Distribution", cex = 0.8, font = 2)
  text(0, -0.10, "(Bin = 1 MB)", cex = 0.6, font = 1)
  
  bin_size <- 1e6
  track_idxs <- setNames(seq(from = 4, to = 1), new_class_types)
  
  # named list of hist colors
  cat_colors <- c("red", "blue", "darkgreen", "black")
  hist_colors <- setNames(as.list(cat_colors), new_class_types)
  
  for (class in new_class_types) {
    current_mat <- data_dfs_store[[class]]
    
    # bin data
    snp_data <- data.frame(
      chr = current_mat$chr, 
      bp = current_mat$bp
    )
    
    snp_data$bin_start <- as.integer(floor(snp_data$bp / bin_size) * bin_size)
    snp_data$bin_end <- as.integer(snp_data$bin_start + bin_size)
    
    # count snps per bin per chromosome
    binned_counts <- snp_data %>%
      group_by(chr, bin_start, bin_end) %>%
      summarise(count = n(), .groups = "drop") %>%
      dplyr::rename(
        Chromosome = chr, 
        chromStart = bin_start, 
        chromEnd = bin_end, 
        Data = count
      ) # need to specify dplyr
    
    binned_counts$Chromosome <- paste0("chr", binned_counts$Chromosome)
    binned_counts <- as.data.frame(binned_counts)
    
    # trim bins to chr widths
    chr_limits <- cyto %>%
      group_by(Chromosome) %>%
      dplyr::summarise(
        chr_start = min(chromStart),
        chr_end = max(chromEnd)
      )
    
    binned_counts_trimmed <- binned_counts %>%
      left_join(chr_limits, by = "Chromosome") %>%
      mutate(
        chromStart = pmax(chromStart, chr_start),
        chromEnd = pmin(chromEnd, chr_end)
      ) %>%
      dplyr::select(-chr_start, -chr_end) # need to specify dplyr
    
    # set color
    binned_counts_trimmed$PlotColor <- hist_colors[[class]]
    
    # plot data
    RCircos.Histogram.Plot.Custom(
      hist.data = binned_counts_trimmed,
      data.col = 4,
      track.num = track_idxs[[class]],
      side = "out"
    )
  }
  
  # add legend for mutation categories
  legend(
    "bottomleft", 
    legend = new_class_types,
    fill = cat_colors,
    border = NA,
    bty = "n",
    cex = 0.8
  )
  
  dev.off()
  
  invisible(NULL)
}

main()