### ALS GWAS m6A/m5C Analysis Pipeline

The pipeline presented here is intended to identify high-confidence m6A loss of function, m6A gain of function, and m5C loss of function SNPs based on data from Project MinE, the University of Massachusetts, and AnswerALS, as well as a control SNP distribution. The pipeline may be run to completion with the Bash script run_pipeline.sh.

The original patient data is unable to be provided due to privacy concerns.

Dependencies:
- tidyverse (R)
- data.table (R)
- here (R)
- GenomicFeatures (Bioconductor)
- GenomicRanges (Bioconductor)
- rtracklayer (Bioconductor)
- Biostrings (Bioconductor)
- biomaRt (Bioconductor)
- Python v. 3.9
- CrossMap

Contact nstrawha@purdue.edu with questions.

## Bash Script

The full pipeline may be run with ```bash run_pipeline.sh``` in Linux. All directory creation and deletion is handled dynamically in this script. Pipeline figure outputs may be observed in outputs/.

## R Scripts

R Scripts are divided into four subfolders of r_scripts/: aggregation/, formatting/, analysis/, and plotting/.

# Aggregation

R Scripts in this folder are intended to combine SNP data from Project MinE, the University of Massachusetts, and AnswerALS in a consistent format.

# Formatting

R Scripts in this folder are intended to manipulate the data in some way prior to analysis. 

[`lmm_to_bed_convert.R`](formatting/lmm_to_bed_convert.R) converts all MinE, UMass, and AnswerALS data from the LMM format into BED format for realignment with hg38 with CrossMap.

[`bed_to_lmm_convert.R`](formatting/bed_to_lmm_convert.R) converts all MinE, UMass, and AnswerALS data from the BED format into LMM format after realignment with hg38 with CrossMap.

[`find_duplicate_snps.R`](formatting/find_duplicate_snps.R) filters for SNPs that are present in two or more datasets, averaging effect size, standard error, and p-value between them. 

# Analysis

R Scripts in this folder expand on the SNP info presented in the original data and draw conclusions.

[`annotate_snps.R`](analysis/annotate_snps.R) uses the information in biomaRt to annotate each SNP with corresponding Ensembl genes, gene symbols, phenotype, and clinical significance (if known).

[`get_lof_snps.R`](analysis/get_lof_snps.R) uses BED format peak data to find high-confidence m6A and m5C loss of function SNPs.

[`get_gof_snps.R`](analysis/get_gof_snps.R) uses the DRACH motif to find high-confidence m6A gain of function SNPs.

[`find_als_peaks.R`](analysis/find_als_peaks.R) uses a control SNP distribution to identify bins of increased mutation activity in the ALS SNP distribution (for all SNPs, m6A LOF, m6A GOF, and m5C LOF).

[`find_high_conf_snps.R`](analysis/find_high_conf_snps.R) calls SNPs within the previously identified peaks and filters by original LMM p-value to identify high-confidence SNPs for experimental validation.