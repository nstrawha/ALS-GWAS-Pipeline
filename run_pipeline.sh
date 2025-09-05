#!/bin/bash
set -euo pipefail

# -------------------------
# Description
# -------------------------
# Pipeline to format and analyze ALS SNP data (Project MinE, UMass, AnswerALS)
# Author: Noah Strawhacker
# Date: Jul. 2025
# Version: 1.0

# Inputs:
# - summary_statistics_2016/: MinE SNP data from 2016 (hg19)
# - summary_statistics_2018/: MinE SNP data from 2018 (hg19)
# - umass_data/: UMass SNP data (hg19)
# - m6as_m5cs/: peak data for m6A and m5C mutations (BED, hg38)

# Outputs:
# - summary_statistics_categorized/: fully-categorized SNP data into mutation type and function
# - functions_barplot.png: frequency barplot of SNPs in each type of functional region
# - mutation_categories_barplot.png: frequency barplot of each mutation category of SNP

# -------------------------
# Folder and function configuration
# -------------------------

touch .here

SCRIPT_DIR="r_scripts"                                      # Original
AGG_DIR="aggregation"                                       # .
FORM_DIR="formatting"                                       # .
ANA_DIR="analysis"                                          # .
PLT_DIR="plotting"                                          # .
DATA_DIR="original_data"                                    # . 
COMBINED_DIR="summary_statistics_combined"                  # Intermediate
BED_HG19_DIR="summary_statistics_combined_bed_format"       # .
BED_HG38_DIR="summary_statistics_combined_bed_format_hg38"  # .
LMM_HG38_DIR="summary_statistics_hg38_orig"                 # .
LMM_HG19_DIR="summary_statistics_hg19_orig"                 # .
HG38_DIR="summary_statistics_hg38"                          # .
HG19_DIR="summary_statistics_hg19"                          # .
DUP_DIR="duplicate_snps"                                    # .
UNIQ_DIR="unique_snps"                                      # .
ANN_DIR="hg38_annotations_by_type"                          # .
SNP_ANN_DIR="summary_statistics_annotated"                  # .
OUT_DIR="outputs"                                           # Output
CAT_DIR="summary_statistics_categorized"                    # . 
INTRON_DIR="intron"                                         # .
EXON_DIR="exon"                                             # .
CDS_DIR="cds"                                               # .
UTR3_DIR="utr3"                                             # .
UTR5_DIR="utr5"                                             # .
PROM_DIR="promoter"                                         # .
TERM_DIR="terminator"                                       # .
UNAN_DIR="unannotated"                                      # .
M6A_LOF_DIR="m6a_lof"                                       # .
M6A_GOF_DIR="m6a_gof"                                       # .
M5C_LOF_DIR="m5c_lof"                                       # .
OTH_MUTS_DIR="other"                                        # .
ALL_MUTS_DIR="all"                                          # .
BIN_DIR="binned_data"                                       # .
SIG_BIN_DIR="significant_bins"                              # .
SIG_SNP_DIR="significant_snps"                              # .
INFSUM_DIR="info_summary"                                   # .
INF_DIR="pheno_or_clinsig_info"                             # .
NOINF_DIR="no_info"                                         # .

# package certain folders for later iteration
DUP_DIRS=(
    "${HG38_DIR}/${DUP_DIR}"
    "${HG38_DIR}/${UNIQ_DIR}"
    "${HG19_DIR}/${DUP_DIR}"
    "${HG19_DIR}/${UNIQ_DIR}"
)

ANN_DIRS=(
    $INTRON_DIR 
    $EXON_DIR 
    $CDS_DIR 
    $UTR3_DIR 
    $UTR5_DIR 
    $PROM_DIR 
    $TERM_DIR 
    $UNAN_DIR
)

CAT_DIRS=(
    $M6A_LOF_DIR
    $M6A_GOF_DIR
    $M5C_LOF_DIR
    $OTH_MUTS_DIR
)

INF_DIRS=(
    $INFSUM_DIR
    $INF_DIR
    $NOINF_DIR
)

# function to clean intermediate directories
clean_up() {
    local KEEP_OUTPUTS=$1
    echo "Cleaning up..."

    if [[ "$KEEP_OUTPUTS" = true ]]; then
        for dir in */; do
            dirname="${dir%/}"
            if [[ "$dirname" != "$OUT_DIR" && "$dirname" != "$SCRIPT_DIR" && "$dirname" != "$DATA_DIR" && "$dirname" != "$CAT_DIR" && "$dirname" != "$BIN_DIR" && "$dirname" != "$SIG_BIN_DIR" && "$dirname" != "$SIG_SNP_DIR" ]]; then
            echo "Removing $dirname..."
            rm -rf "$dirname"
            fi
        done
    else
        for dir in */; do
            dirname="${dir%/}"
            if [[ "$dirname" != "$SCRIPT_DIR" && "$dirname" != "$DATA_DIR" ]]; then
            echo "Removing $dirname..."
            rm -rf "$dirname"
            fi
        done

        mkdir -p "$OUT_DIR"
    fi
}

# -------------------------
# Main Script
# -------------------------

main () {

    # clean files from previous runs, remake outputs folder
    clean_up false

    # combine and deduplicate MinE data
    echo "Aggregating MinE data..."
    mkdir -p "$COMBINED_DIR"
    Rscript "${SCRIPT_DIR}/${AGG_DIR}/proj_mine_aggregation.R"

    # reformat and combine UMass data
    echo "Aggregating UMass data..."
    Rscript "${SCRIPT_DIR}/${AGG_DIR}/umass_aggregation.R"

    # reformat and combine AnswerALS data
    echo "Aggregating AnswerALS data..."
    Rscript "${SCRIPT_DIR}/${AGG_DIR}/answer_als_aggregation.R"

    # convert files to BED for later realignment
    echo "Converting to BED for realignment..."
    mkdir -p "$BED_HG19_DIR"
    Rscript "${SCRIPT_DIR}/${FORM_DIR}/lmm_to_bed_convert.R"

    # realign files to hg38
    echo "Realigning BED files..."
    mkdir -p "$BED_HG38_DIR"

    for CHR in {1..22}; do
        IN_FILE="${BED_HG19_DIR}/als.sumstats.chr${CHR}.combined.bed"
        OUT_FILE="${BED_HG38_DIR}/als.sumstats.chr${CHR}.combined.bed"
        CHAIN_FILE="${DATA_DIR}/hg19ToHg38.over.chain.gz"

        CrossMap bed $CHAIN_FILE $IN_FILE $OUT_FILE
    done

    # combine with old LMM data
    echo "Data realigned. Converting back to LMM..."
    mkdir -p "$LMM_HG38_DIR"
    mkdir -p "$LMM_HG19_DIR"
    Rscript "${SCRIPT_DIR}/${FORM_DIR}/bed_to_lmm_convert.R"

    echo "Data realigned and formatted"

    # remove intermediates
    rm -r "$COMBINED_DIR"
    rm -r "$BED_HG19_DIR"
    rm -r "$BED_HG38_DIR"

    # find and categorize duplicates
    echo "Categorizing duplicates..."
    for DIR in "${DUP_DIRS[@]}"; do
        mkdir -p "$DIR"
    done
    Rscript "${SCRIPT_DIR}/${FORM_DIR}/find_duplicate_snps.R"

    # remove intermediates
    rm -r "$LMM_HG38_DIR"
    rm -r "$LMM_HG19_DIR"

    # annotate duplicates
    echo "Annotating SNPs..."
    for DIR in "${ANN_DIRS[@]}"; do
        mkdir -p "${SNP_ANN_DIR}/${DIR}"
    done
    mkdir -p "$ANN_DIR"
    Rscript "${SCRIPT_DIR}/${ANA_DIR}/annotate_snps.R"

    # remove intermediates
    for DIR in "${DUP_DIRS[@]}"; do
        rm -r $DIR
    done
    rm -r "$HG38_DIR"
    rm -r "$HG19_DIR"

    # categorize mutation types
    echo "Categorizing LOF mutations..."
    for CDIR in "${CAT_DIRS[@]}"; do
        for ADIR in "${ANN_DIRS[@]}"; do
            mkdir -p "${CAT_DIR}/${CDIR}/${ADIR}"
        done
    done
    Rscript "${SCRIPT_DIR}/${ANA_DIR}/get_lof_snps.R" # LOFs

    echo "Categorizing GOF mutations..."
    Rscript "${SCRIPT_DIR}/${ANA_DIR}/get_gof_snps.R" # m6A GOFS

    # plot mutation categories
    echo "Creating barplot..."
    Rscript "${SCRIPT_DIR}/${PLT_DIR}/mutation_classes_by_function.R"

    CAT_DIRS=(
    $M6A_LOF_DIR
    $M6A_GOF_DIR
    $M5C_LOF_DIR
    $ALL_MUTS_DIR
    )

    # find highest confidence peaks
    echo "Finding high-confidence peaks..."
    for CDIR in "${CAT_DIRS[@]}"; do
        mkdir -p "${BIN_DIR}/${CDIR}"
        mkdir -p "${SIG_BIN_DIR}/${CDIR}"
    done
    Rscript "${SCRIPT_DIR}/${ANA_DIR}/find_als_peaks.R"

    # make rcircos plot
    echo "Creating RCircos plot..."
    Rscript "${SCRIPT_DIR}/${PLT_DIR}/rcircos_histogram.R"

    # make linear histograms
    echo "Creating linear histograms..."
    Rscript "${SCRIPT_DIR}/${PLT_DIR}/snp_dist_histograms.R"

    # call out significant snps
    echo "Identifying significant SNPs"
    for IDIR in "${INF_DIRS[@]}"; do
        mkdir -p "${SIG_SNP_DIR}/${IDIR}"
    done
    Rscript "${SCRIPT_DIR}/${ANA_DIR}/find_high_conf_snps.R"

    # clean all folders except outputs/
    clean_up true
    rm -f .here

    echo "Script finished. Outputs written in ${OUT_DIR}/."
}

main "$@"