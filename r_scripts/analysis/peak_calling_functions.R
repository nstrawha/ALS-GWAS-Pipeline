#!/usr/bin/env Rscript

# peak_calling_functions.R
# Code by Noah Strawhacker
# Jul. 2025
# Function to call peaks of one distribution relative to another w/ bootstrapping
# Assumes pop and obs each have columns "count" and "bin_start"

bin_data_reg <- function(df, binwidth, chr, class, type) {
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
      chr = chr, 
      rel_freq = count / total_snps,
      bin_mid = bin_start + binwidth / 2, 
      mut_class = class, 
      data_type = type
    )
  
  return(binned_data)
}


bin_data_gene <- function(df, ranges, chr, class, type) {
  
  # restrict gene ranges to this chromosome
  gene_chr <- ranges[seqnames(ranges) == paste0("chr", chr), ]
  
  # sort genes to define igen regions
  gene_chr <- sort(gene_chr)
  
  # get intergenic as the complement of gene ranges
  igen <- gaps(gene_chr)
  
  # remove seqnames and empty ranges outside the chromosome
  igen <- igen[strand(igen) == "*" & width(igen) > 0]
  
  # label them
  igen$gene_id <- paste0("igen_", seq_along(igen))
  igen$gene_name <- "igen"
  
  # build gr from df
  snp_gr <- GRanges(
    seqnames = paste0("chr", chr),
    ranges = IRanges(start = df$bp, end = df$bp)
  )
  
  # combine genes + igen
  bins <- c(gene_chr, igen)
  
  # overlap snps with bins
  hits <- findOverlaps(snp_gr, bins)
  
  # assign SNPs to bins
  snp_bin <- tibble(
    bp = start(snp_gr)[queryHits(hits)],
    gene_id = bins$gene_id[subjectHits(hits)],
    gene_name = bins$gene_name[subjectHits(hits)]
  )
  
  # count SNPs per bin
  total_snps <- nrow(df)
  binned_data <- snp_bin %>%
    group_by(gene_id, gene_name) %>%
    summarise(
      count = n(),
      bin_start = min(bp),
      bin_end = max(bp),
      .groups = "drop"
    ) %>%
    mutate(
      chr = chr,
      rel_freq = count / total_snps,
      bin_mid = (bin_start + bin_end) / 2,
      mut_class = class,
      data_type = type
    )
  
  return(binned_data)
}


call_relative_peaks <- function(pop, obs, nboots, p_cutoff, chr) {
  
  # set up data to sample from
  sample_size <- sum(obs$count)
  bin_starts <- pop$bin_start
  bin_freqs <- pop$count
  
  # set up relative frequencies
  obs <- obs %>%
    mutate(rel_freq = count / sum(count))
  
  pop <- pop %>%
    mutate(rel_freq = count / sum(count))
  
  bin_idxs_to_boot <- rep(bin_starts, bin_freqs)
  
  # initialize output storage
  booted_mat <- matrix(0, nrow = nboots, ncol = length(bin_starts))
  colnames(booted_mat) <- as.character(bin_starts)
  
  message("Bootstrap initialized")
  
  # carry out sampling
  for (sample in 1:nboots) {
    
    if (sample %% 100 == 0) message("Sample ", sample, "/", nboots, " for chr ", chr)
    
    boot_sample <- sample(bin_idxs_to_boot, size = sample_size, replace = TRUE)
    idxs <- match(boot_sample, bin_starts)
    counts <- tabulate(idxs, nbins = length(bin_starts))
    
    booted_mat[sample, ] <- counts
  }
  
  # mean and sd of bootstrap counts per bin (columns)
  boot_means <- colMeans(booted_mat)
  boot_sds <- apply(booted_mat, 2, sd)
  
  message("Sampling complete. Conducting analysis")
  
  # obtain pvals for bins
  obs_counts <- setNames(obs$count, obs$bin_start)
  observed_bins <- as.character(obs$bin_start)
  
  pvals <- sapply(observed_bins, function(bin) {
    i <- which(colnames(booted_mat) == bin)
    mean(booted_mat[, i] >= obs_counts[bin])
  })
  
  # adjust pvals for multiple comparisons
  adj_pvals <- p.adjust(pvals, method = "BH")
  
  significant_pvals <- adj_pvals[adj_pvals <= p_cutoff]
  
  # check for no significant peaks
  if (length(significant_pvals) == 0) {
    message("No significant peaks found. Returning empty data table")
    
    empty_dt <- data.table(
      bin_start	= numeric(), 
      bin_end = numeric(), 
      count = integer(), 
      chr	= integer(), 
      rel_freq = numeric(),	
      bin_mid	= integer(), 
      mut_class	= character(), 
      pval = numeric(), 
      pop_rel_freq = numeric(), 
      relfreq_pct_diff = numeric(), 
      conf_lvl = numeric()
    )
    
    return(empty_dt)
  }
  
  sig_bin_names <- names(significant_pvals)
  
  # filter obs and append pval + relative frequency
  peak_regions <- obs %>%
    filter(bin_start %in% as.numeric(sig_bin_names)) %>%
    mutate(
      pval = significant_pvals[as.character(bin_start)]
    ) %>%
    left_join(
      as.data.table(pop)[, .(bin_start, pop_rel_freq = rel_freq)],
      by = "bin_start"
    ) %>%
    mutate(
      relfreq_pct_diff = ifelse(pop_rel_freq == 0, NA_real_,
                                ((rel_freq - pop_rel_freq) / pop_rel_freq) * 100)
    )
  
  peak_regions <- peak_regions %>% 
    arrange(desc(relfreq_pct_diff))
  
  return(peak_regions)
}