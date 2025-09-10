#!/bin/bash
# Download File
# set once, reuse below
ACCESSION="SRR24339408"

# reuse your target folders, set one time and reuse for the same session
GEOROOT="/Users/rusher/Desktop/Li_Lab/GEO_data"
mkdir -p "$GEOROOT/sra_cache" "$GEOROOT/fastq"

# 1 download the .sra
prefetch "$ACCESSION" -O "$GEOROOT/sra_cache"

# 2 convert to FASTQ (paired-end split, 8 threads, pipeline mode)
fasterq-dump "$ACCESSION" --split-files -e 8 -p -O "$GEOROOT/fastq"

# 3 (optional) compress
gzip "$GEOROOT/fastq"/${ACCESSION}*.fastq

# 4  check
ls -lh "$GEOROOT/fastq" | grep "$ACCESSION"

------------------------------------------
## Trimming Notes (bash terminal)

# whole path 
/Users/rusher/Desktop/Trimmomatic-0.40/adapters/TruSeq3-PE.fa

# set variables once
ACCESSION="SRR24339408"
GEOROOT="/Users/rusher/Desktop/Li_Lab/GEO_data/fastq"
ADAPTERS="$HOME/Desktop/Trimmomatic-0.40/adapters"

# (optional) sanity checks
cd "$GEOROOT"
ls ${ACCESSION}_1.fastq.gz ${ACCESSION}_2.fastq.gz
ls "$ADAPTERS/TruSeq3-PE.fa"

# trim (paired-end) – note: the command starts on this line and every line ends with "\" (no trailing spaces)
java -jar ~/Desktop/Trimmomatic-0.40/trimmomatic-0.40.jar PE -threads 8 \
  ${ACCESSION}_1.fastq.gz ${ACCESSION}_2.fastq.gz \
  ${ACCESSION}_1_paired.fq.gz ${ACCESSION}_1_unpaired.fq.gz \
  ${ACCESSION}_2_paired.fq.gz ${ACCESSION}_2_unpaired.fq.gz \
  ILLUMINACLIP:${ADAPTERS}/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  -summary ${ACCESSION}_trim_summary.txt

# save in the same folder
mkdir "${ACCESSION}_trim_data"
mv ${ACCESSION}_1_paired.fq.gz ${ACCESSION}_1_unpaired.fq.gz ${ACCESSION}_2_paired.fq.gz ${ACCESSION}_2_unpaired.fq.gz ${ACCESSION}_trim_summary.txt "${ACCESSION}_trim_data/"


# Parameters:
#	ILLUMINACLIP → removes Illumina adapters (uses Trimmomatic’s TruSeq3-PE.fa).
#	LEADING:3 / TRAILING:3 → trim bases below quality 3 at ends.
#	SLIDINGWINDOW:4:15 → cut when avg quality < 15 in a 4-base window.
#	MINLEN:36 → discard reads shorter than 36 bp.


------------------------------
# run FastQC on trimmed data
ACCESSION="SRR24339409"
GEOROOT="/Users/rusher/Desktop/Li_Lab/GEO_data/fastq/${ACCESSION}_trim_data"

# make output folder
mkdir -p "${GEOROOT}/fastqc"

# run FastQC
~/Desktop/FastQC/fastqc -t 8 \
  ${GEOROOT}/${ACCESSION}_1_paired.fq.gz \
  ${GEOROOT}/${ACCESSION}_2_paired.fq.gz \
  -o "${GEOROOT}/fastqc"