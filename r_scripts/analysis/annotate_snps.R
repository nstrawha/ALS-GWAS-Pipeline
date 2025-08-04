#!/usr/bin/env Rscript

# annotate_snps.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to annotate snp location purposes after finding duplicates

# load packages
library(here)
library(tidyverse)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)

# function to categorize snps based on overlaps
categorize_snps <- function(flist, snps_gr, snps_df, fgrs) {
  snps_list <- list()
  
  for (f in flist) {
    current_fgr <- fgrs[[f]]
    
    hits <- findOverlaps(snps_gr, current_fgr, type = "within")
    snps <- unique(snps_df[queryHits(hits), ])
    snps_list[[f]] <- snps
  }
  
  # get unannotated snps
  categorized_snps <- rbindlist(snps_list, use.names = TRUE, fill = TRUE)
  snps_list[["unannotated"]] <- fsetdiff(snps_df, categorized_snps)
  
  return(snps_list)
}


main <- function() {
  
  # load marts
  snp_mart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
  ensembl_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # set params
  num_chromosomes <- 22
  
  functions_list <- c(
    "intron", 
    "exon", 
    "cds", 
    "promoter", 
    "terminator", 
    "utr3", 
    "utr5", 
    "unannotated"
    )
  
  # load hg38 gtf as txdb
  annotation_address  <- here("original_data", "hg38.ensGene.gtf")
  txdb_hg38 <- makeTxDbFromGFF(annotation_address, format = "gtf")
  
  # get gene info
  gene_info <- genes(txdb_hg38)
  gene_ids <- names(gene_info)
  
  # get linked ensembl ids and hgnc names
  gene_annotations <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = ensembl_mart
  )
  
  # load separate annotation types
  intron_gr      <- intronsByTranscript(txdb_hg38) # no strand info yet
  exon_gr        <- exons(txdb_hg38)
  # promoter_gr    <- promoters(txdb_hg38) # BROKEN
  # terminator_gr  <- terminators(txdb_hg38) # BROKEN
  cds_gr         <- cds(txdb_hg38)
  utr3_gr        <- threeUTRsByTranscript(txdb_hg38)
  utr5_gr        <- fiveUTRsByTranscript(txdb_hg38)
  unannotated_gr <- GRanges() # dummy
  
  intron_gr <- unlist(intron_gr, use.names = FALSE)
  utr3_gr   <- unlist(utr3_gr, use.names = FALSE)
  utr5_gr   <- unlist(utr5_gr, use.names = FALSE)
  
  # PROMOTERS AND TERMINATORS FUNCS NOT INCLUDED IN PIPELINE R VER
  # build manually
  tx <- transcripts(txdb_hg38)
  tss <- ifelse(strand(tx) == "+", start(tx), end(tx))
  tes <- ifelse(strand(tx) == "+", end(tx), start(tx))
  
  promoter_gr <- GRanges(
    seqnames = seqnames(tx),
    ranges = IRanges(
      start = ifelse(strand(tx) == "+", tss - 2000, tss - 200),
      end = ifelse(strand(tx) == "+", tss + 200, tss + 2000)
    ),
    strand = strand(tx)
  )
  
  terminator_gr <- GRanges(
    seqnames = seqnames(tx),
    ranges = IRanges(
      start = ifelse(strand(tx) == "+", tes - 2000, tes - 200),
      end = ifelse(strand(tx) == "+", tes + 200, tes + 2000)
    ),
    strand = strand(tx)
  )
  
  # package all functional region grs
  func_grs_list <- setNames(
    as.list(mget(paste0(functions_list, "_gr"))),
    functions_list
  )
  
  # set output file locs
  output_flocs <- setNames(
    lapply(
      functions_list, 
      function(f) here("summary_statistics_annotated", f)
    ), 
    functions_list
  )
  
  # initialize freqs of each annotation type
  n_annotations <- length(functions_list)
  ffreqs_store <- setNames(
    numeric(n_annotations),
    functions_list
  )
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Processing chr ", chrom_idx)
    
    # access original data
    old_faddress <- here(
      "summary_statistics_hg38", 
      "duplicate_snps", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".hg38.dups.txt")
      )
    
    old_data <- fread(file = old_faddress)
    
    # annotate snp info
    rs_vec <- old_data$snp
    snp_info <- getBM(
      attributes = c(
        "refsnp_id",
        "phenotype_description", 
        "clinical_significance"
        ), 
      filters = "snp_filter", 
      values = rs_vec, 
      mart = snp_mart
    )
    
    old_data <- old_data %>%
      left_join(snp_info, by = c("snp" = "refsnp_id"))
    
    # set up granges
    data_gr <- GRanges(
      seqnames = paste0("chr", chrom_idx), 
      ranges = IRanges(start = old_data$bp, end = old_data$bp)
      )
    
    col_names <- colnames(old_data)
    
    # get ensembl ids and gene names
    gene_ovlps <- findOverlaps(data_gr, gene_info)
    snp_hits <- data_gr[queryHits(gene_ovlps)]
    gene_hits <- gene_info[subjectHits(gene_ovlps)]
    
    snp_gene_df <- data.frame(
      snp = old_data$snp[queryHits(gene_ovlps)],
      ensembl_gene_id = names(gene_hits)
    )
    
    # merge back into old data
    snp_gene_annotated <- left_join(snp_gene_df, gene_annotations, by = c("ensembl_gene_id" = "ensembl_gene_id"))
    snp_gene_annotated <- snp_gene_annotated %>%
      mutate(
        hgnc_symbol = ifelse(is.na(hgnc_symbol), "NA", hgnc_symbol)
      )
    
    old_data <- left_join(old_data, snp_gene_annotated, by = "snp")
    
    # set all na entries to "none"
    old_data <- old_data %>%
      mutate(across(everything(), ~ ifelse(. == "" | is.na(.), "none", .)))
    
    # categorize all snps for each function
    all_snps <- categorize_snps(
      flist = functions_list, 
      snps_gr = data_gr, 
      snps_df = old_data, 
      fgrs = func_grs_list
      )
    
    # set output file names and addresses
    output_fnames <- setNames(
      lapply(
        functions_list, 
        function(f) paste0("als.sumstats.lmm.chr", chrom_idx, ".", f, ".txt")
        ),
      functions_list
    )
    
    output_faddresses <- setNames(
      mapply(file.path, output_flocs, output_fnames, SIMPLIFY = FALSE), 
      functions_list
    )
    
    # write output files
    for (f in functions_list) {
      
      current_snps <- all_snps[[f]]
      current_snps$func <- f
      
      # ensure directory is created
      dir.create(
        dirname(output_faddresses[[f]]), 
        recursive = TRUE, 
        showWarnings = FALSE
      )
      
      # write output
      fwrite(
        current_snps, 
        file = output_faddresses[[f]], 
        sep = "\t", 
        col.names = TRUE, 
        quote = FALSE)
    }
    
    # record freqs of each annotation type
    for (f in functions_list) {
      ffreqs_store[[f]] <- ffreqs_store[[f]] + nrow(unique(all_snps[[f]]))
    }
  }
  
  # create bar plot
  ffreqs_store <- unlist(ffreqs_store)
  rel_ffreqs <- ffreqs_store / sum(ffreqs_store)
  
  labels <- c(
    "Intron", 
    "Exon", 
    "CDS", 
    "Promoter", 
    "Terminator", 
    "3'UTR", 
    "5'UTR", 
    "N/A"
    )
  
  # structure data for plotting
  df_to_plot <- data.frame(
    category = labels,
    count = ffreqs_store,
    rel_freq = rel_ffreqs
  )
  
  # plot freq data
  hist <- ggplot(df_to_plot, aes(x = category, y = rel_freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0("n = ", count)), vjust = -0.5, size = 4) +
    ylab("Relative Frequency") +
    xlab("Functional Region")
  
  # save resulting freq plot
  ggsave(
    filename = here("outputs", "snp_function_barplot.png"), 
    plot = hist, 
    width = 10, height = 5
  )
  
  invisible(NULL)
}

main()
