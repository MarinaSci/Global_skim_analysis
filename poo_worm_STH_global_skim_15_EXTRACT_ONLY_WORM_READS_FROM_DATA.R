#9 June 2025

#Pulling the worm reads from all the 150 samples we sequenced as part of this study: 

/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING/05_KEEPING_ONLY_WORM_READS

#Steps:
#List of the 150 files only 
#Ask the script to only look at those 150
#Extract reads names from bam files when samples were mapped to mitochondrial genomes
#Extract reads names from bam files when samples were mapped to NUCLEAR genomes
#Concatenate both files with all read names
#Extract the reads from the RAW FILES (before trimming)
#zip the files
#ready for upload


#!/bin/bash

source ~/conda/etc/profile.d/conda.sh

# Define paths
SAMPLE_LIST="150_SAMPLES_NAME_LIST"
RAW_DIR="/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/01.RAW_DATA/NOVOGENE_DATA_TO_PULL_WORM_READS"
ANALYSIS_DIR="/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING/05_KEEPING_ONLY_WORM_READS"
MT_MAPPING_DIR="/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING/"
NUC_MAPPING_DIR="/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/02.NUCLEAR_MAPPING"

# Loop through each sample in the list
while read -r sample_prefix; do
echo "🔄 Processing $sample_prefix"

mito_bam="$MT_MAPPING_DIR/${sample_prefix}_trimmed.bam"
nuc_bam="$NUC_MAPPING_DIR/${sample_prefix}_trimmed.bam"

# Check BAM files exist
if [[ ! -f "$mito_bam" || ! -f "$nuc_bam" ]]; then
echo "⚠️  Missing BAM file(s) for $sample_prefix, skipping."
continue
fi

# Activate samtools
conda activate samtools
# Extract read names and clean suffixes
samtools view "$mito_bam" | cut -f1 | sed 's/\/[12]$//' | sort | uniq > "$ANALYSIS_DIR/${sample_prefix}_mito_readnames.txt"
samtools view "$nuc_bam" | cut -f1 | sed 's/\/[12]$//' | sort | uniq > "$ANALYSIS_DIR/${sample_prefix}_nuclear_readnames.txt"

conda deactivate

# Combine and deduplicate
cat "$ANALYSIS_DIR/${sample_prefix}_mito_readnames.txt" "$ANALYSIS_DIR/${sample_prefix}_nuclear_readnames.txt" \
| sort | uniq > "$ANALYSIS_DIR/${sample_prefix}_readnames.txt"

# Check that readnames file is not empty
if [[ ! -s "$ANALYSIS_DIR/${sample_prefix}_readnames.txt" ]]; then
echo "❌ No read names found for $sample_prefix — skipping read extraction."
continue
fi

# Extract matching reads using seqtk
conda activate seqtk

seqtk subseq "$RAW_DIR/${sample_prefix}_1.fq.gz" "$ANALYSIS_DIR/${sample_prefix}_readnames.txt" | gzip > "$ANALYSIS_DIR/${sample_prefix}_WORM_READS_1.fq.gz"
seqtk subseq "$RAW_DIR/${sample_prefix}_2.fq.gz" "$ANALYSIS_DIR/${sample_prefix}_readnames.txt" | gzip > "$ANALYSIS_DIR/${sample_prefix}_WORM_READS_2.fq.gz"

conda deactivate
# Step 3: Verify read counts match
reads1=$(zcat "$ANALYSIS_DIR/${sample_prefix}_WORM_READS_1.fq.gz" | wc -l)
reads2=$(zcat "$ANALYSIS_DIR/${sample_prefix}_WORM_READS_2.fq.gz" | wc -l)

count1=$((reads1 / 4))
count2=$((reads2 / 4))

if [ "$count1" -eq "$count2" ]; then
echo "✅ $sample_prefix: Read count check passed ($count1 reads)"
else
  echo "❌ $sample_prefix: Read count mismatch! R1: $count1, R2: $count2"
fi

echo ""


#done 






done < "$SAMPLE_LIST"