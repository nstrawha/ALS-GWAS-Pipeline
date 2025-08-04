#!/usr/bin/env Rscript

# peak_calling_functions.R
# Code by Noah Strawhacker
# Jul. 2025
# Function to call peaks of one distribution relative to another w/ bootstrapping
# Assumes pop and obs each have columns "count" and "bin_start"

bin_data <- function(df, binwidth) {
  snp_data <- data.frame(
    bp = df$bp
  )
    
  total_snps <- nrow(snp_data)
  
  snp_data$bin_start <- as.integer(floor(snp_data$bp / binwidth) * binwidth)
  snp_data$bin_end <- as.integer(snp_data$bin_start + binwidth)
  
  binned_data <- snp_data %>%
    group_by(bin_start, bin_end) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      rel_freq = count / total_snps,
      bin_mid = bin_start + binwidth / 2
    )
  
  return(binned_data)
}


call_relative_peaks <- function(pop, obs, nboots, p_cutoff) {
  
  # set up data to sample from
  sample_size <- sum(obs$count)
  bin_starts <- pop$bin_start
  bin_freqs <- pop$count
  
  df_to_boot <- rep(bin_starts, bin_freqs)
  
  # initialize output storage
  booted_df <- matrix(0, nrow = nboots, ncol = length(bin_starts))
  colnames(booted_df) <- as.character(bin_starts)
  
  message("Bootstrap initialized")
  
  # carry out sampling
  for (sample in 1:nboots) {
    
    if (sample %% 50 == 0) message("Sample ", sample, "/", nboots)
    
    boot_sample <- sample(df_to_boot, size = sample_size, replace = TRUE)
    idxs <- match(boot_sample, bin_starts)
    counts <- tabulate(idxs, nbins = length(bin_starts))
    
    booted_df[sample, ] <- counts
  }
  
  message("Sampling complete. Conducting analysis")
  
  # obtain pvals for bins
  obs_counts <- setNames(obs$count, obs$bin_start)
  observed_bins <- as.character(obs$bin_start)
  
  pvals <- sapply(observed_bins, function(bin) {
    i <- which(colnames(booted_df) == bin)
    mean(booted_df[, i] >= obs_counts[bin])
  })
  
  # adjust pvals for multiple comparisons
  adj_pvals <- p.adjust(pvals, method = "BH")
  
  significant_pvals <- adj_pvals[adj_pvals <= p_cutoff]
  sig_bin_names <- names(significant_pvals)
  
  peak_regions <- pop %>%
    filter(bin_start %in% as.numeric(sig_bin_names)) %>%
    mutate(pval = significant_pvals[as.character(bin_start)])
  
  return(peak_regions)
}