#!/usr/bin/env Rscript

# peak_calling_functions.R
# Code by Noah Strawhacker
# Jul. 2025
# Function to call peaks of one distribution relative to another w/ bootstrapping
# Assumes pop and obs each have columns "count" and "bin_start"


bin_data_gene <- function(df, ranges, chr, class, type, chunk_size = 1e5) {
  gene_chr <- ranges[seqnames(ranges) == paste0("chr", chr), ]
  gene_chr <- sort(gene_chr)
  igen <- gaps(gene_chr)
  igen <- igen[strand(igen) == "*" & width(igen) > 0]
  igen$gene_id <- paste0("igen_", seq_along(igen))
  igen$gene_name <- "igen"
  bins <- c(gene_chr, igen)
  
  total_snps <- nrow(df)
  message("Total SNPs: ", total_snps)
  
  # Split SNPs into chunks
  snp_chunks <- split(df, ceiling(seq_len(nrow(df)) / chunk_size))
  message("Processing ", length(snp_chunks), " chunks...")
  
  binned_list <- vector("list", length(snp_chunks))
  
  for (i in seq_along(snp_chunks)) {
    message("Chunk ", i, "/", length(snp_chunks))
    
    snp_df <- snp_chunks[[i]]
    snp_gr <- GRanges(
      seqnames = paste0("chr", chr),
      ranges = IRanges(start = snp_df$bp, end = snp_df$bp)
    )
    
    hits <- findOverlaps(snp_gr, bins)
    
    if (length(hits) == 0L) next
    
    snp_bin <- tibble(
      bp = start(snp_gr)[queryHits(hits)],
      gene_id = bins$gene_id[subjectHits(hits)],
      gene_name = bins$gene_name[subjectHits(hits)]
    )
    
    binned_chunk <- snp_bin %>%
      group_by(gene_id, gene_name) %>%
      summarise(
        count = n(),
        bin_start = min(bp),
        bin_end = max(bp),
        .groups = "drop"
      )
    
    binned_list[[i]] <- binned_chunk
  }
  
  binned_data <- bind_rows(binned_list) %>%
    group_by(gene_id, gene_name) %>%
    summarise(
      count = sum(count),
      bin_start = min(bin_start),
      bin_end = max(bin_end),
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



call_relative_peaks <- function(pop, obs, nboots, chr) {
  
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
  
  # force counts to same bins
  all_bins <- union(pop$bin_start, obs$bin_start)
  pop <- pop %>% complete(bin_start = all_bins, fill = list(count = 0))
  obs <- obs %>% complete(bin_start = all_bins, fill = list(count = 0))
  
  
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
  significant_pvals <- adj_pvals
  
  # check for no significant peaks
  if (length(significant_pvals) == 0) {
    message("No significant peaks found. Returning empty data table")
    
    empty_dt <- data.table(
      bin_start  = numeric(), 
      bin_end = numeric(), 
      count = integer(), 
      chr  = integer(), 
      rel_freq = numeric(),  
      bin_mid  = integer(), 
      mut_class  = character(), 
      pval = numeric(), 
      pop_rel_freq = numeric(), 
      relfreq_pct_diff = numeric(), 
      conf_lvl = numeric()
    )
    
    return(empty_dt)
  }
  
  sig_bin_names <- names(significant_pvals)
  
  # Aggregate pop by bin_start to avoid many-to-many joins
  pop_summary <- as.data.table(pop)[, .(pop_rel_freq = mean(rel_freq)), by = bin_start]
  
  # filter obs and append pval + relative frequency
  peak_regions <- obs %>%
    filter(bin_start %in% as.numeric(sig_bin_names)) %>%
    mutate(pval = significant_pvals[as.character(bin_start)]) %>%
    left_join(pop_summary, by = "bin_start") %>%
    mutate(
      relfreq_pct_diff = ifelse(pop_rel_freq == 0, NA_real_,
                                ((rel_freq - pop_rel_freq) / pop_rel_freq) * 100)
    )
  
  peak_regions <- peak_regions %>% 
    arrange(desc(relfreq_pct_diff))
  
  return(peak_regions)
}
