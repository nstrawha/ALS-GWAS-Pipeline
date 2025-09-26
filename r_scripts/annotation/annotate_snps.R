#!/usr/bin/env Rscript

# annotate_snps.R
# Code by Noah Strawhacker
# Jun. 2025
# Script to annotate snp location purposes after finding duplicates


# Setup and functions -----------------------------------------------------

rm(list = ls())

# load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(biomaRt))

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


# Main -------------------------------------------------------------------

main <- function() {
  

  # Loading the annotation info -------------------------------------------

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
  cds_gr         <- cds(txdb_hg38)
  utr3_gr        <- threeUTRsByTranscript(txdb_hg38)
  utr5_gr        <- fiveUTRsByTranscript(txdb_hg38)
  unannotated_gr <- GRanges() # dummy
  
  intron_gr <- unlist(intron_gr, use.names = FALSE)
  utr3_gr   <- unlist(utr3_gr, use.names = FALSE)
  utr5_gr   <- unlist(utr5_gr, use.names = FALSE)
  
  # promoters and terminators functions are not included in this R version
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
  
  # set output file loc
  output_floc <- "summary_statistics_annotated"
  
  # initialize freqs of each annotation type
  n_annotations <- length(functions_list)
  ffreqs_store <- setNames(
    numeric(n_annotations),
    functions_list
  )
  

  # Carrying out annotation -------------------------------------------------
  
  # iterate through chromosomes
  for (chrom_idx in 1:num_chromosomes) {
    
    message("Annotating chr ", chrom_idx)
    
    # access original data
    in_faddress <- here(
      "summary_statistics_hg38", 
      paste0("als.sumstats.lmm.chr", chrom_idx, ".combined.hg38.csv")
      )
    
    in_data <- fread(file = in_faddress)
    
    # annotate snp info
    rs_vec <- in_data$snp
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
    
    in_data <- in_data %>%
      left_join(snp_info, by = c("snp" = "refsnp_id"))
    
    # set up granges
    data_gr <- GRanges(
      seqnames = paste0("chr", chrom_idx), 
      ranges = IRanges(start = in_data$bp, end = in_data$bp)
      )
    
    col_names <- colnames(in_data)
    
    # get ensembl ids and gene names
    gene_ovlps <- findOverlaps(data_gr, gene_info)
    snp_hits <- data_gr[queryHits(gene_ovlps)]
    gene_hits <- gene_info[subjectHits(gene_ovlps)]
    
    snp_gene_df <- data.frame(
      snp = in_data$snp[queryHits(gene_ovlps)],
      ensembl_gene_id = names(gene_hits)
    )
    
    # merge back into old data
    snp_gene_annotated <- left_join(snp_gene_df, gene_annotations, by = c("ensembl_gene_id" = "ensembl_gene_id"))
    snp_gene_annotated <- snp_gene_annotated %>%
      mutate(
        hgnc_symbol = ifelse(is.na(hgnc_symbol), "NA", hgnc_symbol)
      )
    
    in_data <- left_join(in_data, snp_gene_annotated, by = "snp")
    
    # set all na entries to "none"
    in_data <- in_data %>%
      mutate(across(everything(), ~ ifelse(. == "" | is.na(.), "none", .)))
    
    # categorize all snps for each function
    all_snps <- categorize_snps(
      flist = functions_list, 
      snps_gr = data_gr, 
      snps_df = in_data, 
      fgrs = func_grs_list
      )
    
    # set output file names and addresses
    output_fname <- paste0("als.sumstats.lmm.chr", chrom_idx, ".annotated.csv")
    output_faddress <- here(output_floc, output_fname)
    
    # add col name to identify functional regions
    for (f in functions_list) {
      current_snps <- all_snps[[f]]
      current_snps$func <- f
      all_snps[[f]] <- current_snps
    }
    
    all_snps <- rbindlist(all_snps)
    
    # write to output file
    dir.create(
      dirname(output_faddress), 
      recursive = TRUE, 
      showWarnings = FALSE
      )
    
    # write output
    fwrite(
      all_snps, 
      file = output_faddress, 
      col.names = TRUE, 
      quote = TRUE
      )
  }
  
  invisible(NULL)
}

main()