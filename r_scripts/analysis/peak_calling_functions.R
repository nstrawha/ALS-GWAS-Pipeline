#!/usr/bin/env Rscript

# peak_calling_functions.R
# Code by Noah Strawhacker
# Jul. 2025
# Function to call peaks of one distribution relative to another w/ bootstrapping
# Assumes pop and obs each have columns "count" and "bin_start"


bin_data_gene <- function(df, ranges, chr, class, type, chunk_size = 1e5) {
  gene_chr <- ranges[seqnames(ranges) == paste0("chr", chr), ]
  gene_chr <- sort(gene_chr)
  igen <- gaps(reduce(gene_chr, ignore.strand = TRUE))
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
    message("n hits:", length(hits))
    message(head(bins))
    
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
      mut_class = class,
      data_type = type
    )
  
  return(binned_data)
}



call_relative_peaks <- function(pop, obs, nboots, chr) {
  
  sample_size <- sum(obs$count)
  
  # make sure bin_start are numeric
  pop$bin_start <- as.numeric(pop$bin_start)
  obs$bin_start <- as.numeric(obs$bin_start)
  
  num_pop_in_obs <- sum(pop$gene_id %in% obs$gene_id)
  
  obs$gene_id <- unlist(obs$gene_id)
  
  pop_to_merge <- data.frame(
    gene_id = unlist(pop$gene_id),
    count = pop$count,
    rel_freq = pop$rel_freq
  )
  
  # materials for bootstrap
  df_to_boot <- obs %>%
    left_join(pop_to_merge, by = "gene_id", suffix = c(".als", ".ctrl"))
  
  sample_size <- sum(df_to_boot$count.als)
  
  boot_record <- data.frame(
    below_count = rep.int(0, nrow(df_to_boot)), 
    above_count = rep.int(0, nrow(df_to_boot))
  )
  
  # na check
  df_to_boot$count.ctrl[is.na(df_to_boot$count.ctrl)] <- 0
  
  # perform bootstrap
  for (sample_idx in 1:nboots) {
    
    # update message
    if (sample_idx %% 100 == 0) message("chr ", chr, ", sample ", sample_idx, "/", nboots)
    
    # sample rows with replacement, weighted by control counts
    boot_df <- df_to_boot[sample(seq_len(nrow(df_to_boot)),
                                 size = sample_size,
                                 replace = TRUE,
                                 prob = df_to_boot$count.ctrl), ]
    
    # aggregate sampled counts by gene_id
    boot_summary <- boot_df %>%
      count(gene_id, name = "boot_count")
    
    # join back to full df_to_boot to align rows
    df_merged <- df_to_boot %>%
      left_join(boot_summary, by = "gene_id") %>%
      mutate(boot_count = replace_na(boot_count, 0))
  
    is_smaller <- df_merged$count.als <= df_merged$boot_count
    is_bigger <- df_merged$count.als > df_merged$boot_count
    
    boot_record$below_count <- boot_record$below_count + is_smaller
    boot_record$above_count <- boot_record$above_count + is_bigger
  }
  
  # calculate ps
  pvals <- boot_record$above_count / (boot_record$below_count + boot_record$above_count)
  pvals <- p.adjust(pvals, method = "BH")
  
  df_to_boot$gene_p <- pvals
  peak_regions <- df_to_boot
  
  return(peak_regions)
}