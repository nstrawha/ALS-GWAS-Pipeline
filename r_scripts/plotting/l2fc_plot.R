#!/usr/bin/env Rscript

# l2fc_plot.R
# Code by Noah Strawhacker
# Sep. 2025
# Script to plot log2 fold changes across mutation categories


# Setup and functions -----------------------------------------------------

rm(list = ls())
options(scipen = 0)

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(coin))


# Main --------------------------------------------------------------------

main <- function() {
  
  num_chromosomes <- 22
  
  in_floc <- "summary_statistics_tmic_ps" # TODO change this
  
  plot_store <- list()
  
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Recording transcriptomic data to plot for chr ", chrom_idx)
    
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
  
  # remove nones
  plot_store <- plot_store[!(plot_store$tmic_p == "none") & !(plot_store$tmic_l2fc == "none"), ]
  plot_store$tmic_l2fc <- as.numeric(plot_store$tmic_l2fc)
  
  # deduplicate
  plot_store <- plot_store[!duplicated(plot_store$snp), ]
  
  # cumulative frequency plot
  cats_to_test <- c(
    "m6a_lof", 
    "m6a_gof", 
    "m5c_lof"
  )
  
  colors_to_plot <- setNames(list(
    "red", 
    "darkgreen", 
    "purple"
  ), cats_to_test)
  
  labels_to_plot <- setNames(list(
    "m6A LOF SNPs", 
    "m6A GOF SNPs", 
    "m5C LOF SNPs"
  ), cats_to_test)
  
  ctrl_counts <- plot_store %>%
    filter(mut_cat == "other") %>% 
    dplyr::count(tmic_l2fc) %>% # something masks dplyr::count
    mutate(
      cum_l2fc = cumsum(n),
      cum_prop = cum_l2fc / sum(n), 
      mut_cat = "other"
      )
  
  for (cat in cats_to_test) {
    
    # compute cumulative sum
    obs_counts <- plot_store %>%
      filter(mut_cat == cat) %>% 
      dplyr::count(tmic_l2fc) %>% # something masks dplyr::count
      mutate(
        cum_l2fc = cumsum(n),
        cum_prop = cum_l2fc / sum(n), 
        mut_cat = cat
        )
    
    all_counts <- rbind(ctrl_counts, obs_counts)
    
    # plot
    ggplot(all_counts, aes(x = tmic_l2fc, y = cum_prop, color = mut_cat)) + 
      geom_line(linewidth = 1, alpha = 0.7, show.legend = TRUE) +
      labs(
        x = "Log2 Fold Change (Control vs. ALS)", 
        y = "Cumulative Proportion", 
        title = "Cumulative Proportion of Log2 Fold Changes for ALS vs. Control"
      ) + 
      scale_color_manual(
        name = "SNP Category", 
        values = c(colors_to_plot[[cat]], "grey"), 
        labels = c(labels_to_plot[[cat]], "Other SNPs")
        ) + 
      theme_minimal()
    
    ggsave(here("outputs", paste0("fold_change_plot_", cat, ".png")), height = 5, width = 8)
  }
  
  # compute statistics
  comparison_stats <- list()
  
  for (cat in cats_to_test) {
    l2fcs_ctrl <- plot_store$tmic_l2fc[plot_store$mut_cat == "other"]
    l2fcs_obs <- plot_store$tmic_l2fc[plot_store$mut_cat == cat]
    
    # compute mean difference
    m1 <- mean(l2fcs_obs)
    m2 <- mean(l2fcs_ctrl)
    
    meandiff <- m1 - m2
    
    # compute effect size
    
    # compute p (wilcoxon rank sum)
    test_results <- wilcox.test(
      x = l2fcs_obs, 
      y = l2fcs_ctrl, 
      alternative = "two.sided"
    )
    
    p <- test_results$p.value
    
    # effect size
    temp_df_ctrl <- data.frame(
      l2fc = l2fcs_ctrl, 
      mut_cat = "other"
    )
    
    temp_df_obs <- data.frame(
      l2fc = l2fcs_obs, 
      mut_cat = cat
    )
    
    temp_df <- rbindlist(list(temp_df_ctrl, temp_df_obs))
    
    effsize <- wilcox_effsize(
      temp_df,
      formula = l2fc ~ mut_cat,
      alternative = "two.sided"
    )$effsize
    
    # store stats
    current_stats_df <- data.frame(
      mean_difference = meandiff,
      effect_size = effsize,
      padj = p
      )
    
    comparison_stats[[cat]] <- current_stats_df
  }
  
  comparison_stats <- rbindlist(comparison_stats)
  rownames(comparison_stats) <- cats_to_test
  
  comparison_stats$padj <- p.adjust(comparison_stats$padj, method = "BH")
  
  fwrite(
    comparison_stats, 
    file = here("outputs", "transcriptomic_analysis.csv"), 
    col.names = TRUE, 
    row.names = TRUE, 
    quote = TRUE
  )
  
  invisible(NULL)
}

main()