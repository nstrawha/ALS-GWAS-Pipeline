#!/usr/bin/env Rscript

# rcircos_histogram.R
# Code by Noah Strawhacker
# Jul. 2025
# Script to create rcircos circular histogram of snp distributions


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
library(here)
library(tidyverse)
library(data.table)
library(RCircos)
source(here("r_scripts", "plotting", "rcircos_modified_functions.R"))


# Main --------------------------------------------------------------------

main <- function() {
  
  # perform pval cut or no
  pval_cut <- FALSE
  
  # set params
  num_chromosomes <- 22
  
  in_floc <- "summary_statistics_tmic_ps" # TODO change this
  
  plot_store <- list()
  
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Recording SNP data for RCircos plot for chr ", chrom_idx)
    
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
  
  # separate into named list
  mutation_cats <- list("other", "m6a_lof", "m6a_gof", "m5c_lof")
  snps_list <- list()
  
  for (cat in mutation_cats) {
    snps_list[[cat]] <- plot_store[plot_store$mut_cat == cat, ]
  }
  
  # perform pval cut
  if (pval_cut) snps_list <- lapply(snps_list, function(df) {df[df$p < 0.05, ]})
  
  # set up genes to plot
  genes_to_plot <- data.frame(
    Chromosome = c("chr9",    "chr21",  "chr1",   "chr16",  "chr17",  "chr6",    "chr19"),
    chromStart = c(27546546,  31659693, 11012654, 31180139, 45784320, 31111223,  797452),
    chromEnd   = c(27573481,  31668931, 11025492, 31191605, 45835828, 31112575,  812312),
    Gene       = c("C9orf72", "SOD1",   "TARDBP", "FUS",    "CRHR1",  "C6orf15", "PTBP1")
    #              All        All       All       All       m6A LOF   m6A GOF    m5C LOF
  )
  

  # Setting up RCircos ----------------------------------------------------
  
  # set up rcircos
  pdf(here("outputs", "rcircos_plot.pdf"), width = 7, height = 7)
  
  # load data for rcircos plot
  data(UCSC.HG38.Human.CytoBandIdeogram)
  
  # load ideogram
  cyto <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(
    cyto.info = cyto, 
    chr.exclude = c("chrX"), # show chrY so we can put freq labels above it
    tracks.inside = 0,
    tracks.outside = 4
  )
  
  # reset plot params
  params <- RCircos.Get.Plot.Parameters()
  params$chr.ideo.pos <- 0.7          # shrink ideogram
  params$hist.width <- 20             # change bar width
  params$track.out.start <- 0.9       # shrink hists around ideogram
  params$track.background <- "white"  # transparent hist bgs
  params$grid.line.color <- "white"   # transparent hist bgs
  params$track.height <- 0.35         # taller histograms
  RCircos.Reset.Plot.Parameters.Custom(new.params = params)
  
  # reset ideogram params
  ideo <- RCircos.Get.Plot.Ideogram()
  
  RCircos.Set.Plot.Area()
  
  # plot genes
  RCircos.Gene.Name.Plot(
    gene.data = genes_to_plot, 
    name.col = 4, 
    track.num = 1,    
    side = "in"
  )
  
  RCircos.Gene.Connector.Plot.Custom(
    genomic.data = genes_to_plot, 
    track.num = 1,
    side = "out"
  )
  
  RCircos.Chromosome.Ideogram.Plot.Custom()
  
  text(0, 0.05, "ALS SNP Distribution", cex = 0.8, font = 2)
  text(0, -0.10, "(Bin = 1 MB)", cex = 0.6, font = 1)
  
  bin_size <- 1e6
  track_idxs <- setNames(seq(from = 1, to = 4), mutation_cats)
  
  # named list of hist colors
  col_vec <- c(
    "purple", 
    "red", 
    "darkgreen", 
    "blue"
    )
  hist_colors <- setNames(as.list(col_vec), mutation_cats)
  
  message("Ideogram generated")

  
# Plot the data -----------------------------------------------------------
  
  max_counts <- list()
  
  for (cat in mutation_cats) {
    current_mat <- snps_list[[cat]]
    
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
    binned_counts_trimmed$PlotColor <- hist_colors[[cat]]
    
    # plot data
    RCircos.Histogram.Plot.Custom(
      hist.data = binned_counts_trimmed,
      data.col = 4,
      track.num = track_idxs[[cat]],
      side = "out"
    )
    
    # record the max bin counts
    max_counts[[cat]] <- max(binned_counts_trimmed$Data)

  }
  
  # add legend for mutation categories
  legend(
    "bottomleft", 
    legend = c(
      "Other SNPs", 
      "m6A LOF SNPs", 
      "m6A GOF SNPs", 
      "m5C LOF SNPs"
      ),
    fill = col_vec,
    border = NA,
    bty = "n",
    cex = 0.8
  )
  

# Frequency axis label ----------------------------------------------------

  # long vertical bar
  vbar_height = length(mutation_cats) * 0.371 - 0.046
  
  segments(
    x0 = -0.025, y0 = 0.9, 
    x1 = -0.025, y1 = 0.9 + vbar_height, 
    col = "black", lwd = 1
    )
  
  # ticks
  for (idx in 1:length(mutation_cats)) {
    segments(
      x0 = -0.025, y0 = 0.9 + 0.371 * (idx - 1), 
      x1 = -0.050, y1 = 0.9 + 0.371 * (idx - 1), 
      col = "black", lwd = 1
    )
    
    segments(
      x0 = -0.025, y0 = 1.225 + 0.371 * (idx - 1), 
      x1 = -0.050, y1 = 1.225 + 0.371 * (idx - 1), 
      col = "black", lwd = 1
    )
  }
  
  # freq labels
  for (idx in 1:length(mutation_cats)) {
    text(
      x = 0, y = 0.9 + 0.371 * (idx - 1),
      labels = "0", 
      cex = 0.30, 
      pos = 2
    )
    
    text(
      x = 0, y = 1.225 + 0.371 * (idx - 1),
      labels = max_counts[[idx]], 
      cex = 0.30, 
      pos = 2
    )
  }
  
  # "freq" label
  text(
    x = 0.065, y = 1 + length(mutation_cats) * 0.371 - 0.046, 
    labels = "Frequency", 
    cex = 0.60, 
    pos = 2
    )
  
  dev.off()
  
  message("Plotting complete")
  
  invisible(NULL)
}

main()