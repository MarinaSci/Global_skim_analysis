# Indexing and mapping 
Author: Marina Papaiakovou, mpapaiakovou[at]gmail.com 

## Contents: 
- Indexing reference files 
- Mapping of trimmed paired reads to indexed mtDNA reference files

```bash

# FULL SCRIPT: INDEX
#------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_index_file_%j.out
#SBATCH --error=job_index_file_%j.err
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

bwa index human_mito_ref.fasta

#------------------------------------------------------------------------------

#FULL SCRIPT: BWA MAPPING ----- 
sbatch <scriptname>.sh
#if you want to run job one from your array, the do 
sbatch --array=1 <scriptname>.sh
#------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --partition=day 
#SBATCH --output=job_bwa_%j_%a.out
#SBATCH --error=job_bwa_%j_%a.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-974

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

#TESTING NOW TO RUN THE SCRIPT FROM THE SAME FOLDER ΑS THE TRIMMED FILES - IT WORKS !!!

#make a list of samples present
#ls -1 *_1.fq.gz | sed 's/_1.fq.gz//' > sample_list.txt # sample list needs to only have sample_name_trimmed

#REF_DIR=/mbl/share/workspaces/groups/genome-skimming/03.GLOBAL_SKIM/04.ANALYSIS
TRIM_DIR=/mbl/share/workspaces/groups/genome-skimming/03.GLOBAL_SKIM/02.TRIMMED_DATA
MAPPING_DIR=/mbl/share/workspaces/groups/genome-skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING


SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt) #I got this working !!!!
echo "I will now run the bwa mem"

bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" human_mito_ref.fasta ${TRIM_DIR}/${SAMPLE}_1.fq.gz ${TRIM_DIR}/${SAMPLE}_2.fq.gz > ${MAPPING_DIR}/${SAMPLE}.sam
#added the -R flag with RG and SM tag so I can easily remove duplicates from the samples down the line

#doublechck with samtools view -H your.bam #and look for the sample ID
echo "get rid of unmapped reads"
samtools view -q 30 -F 4 -S -h ${MAPPING_DIR}/${SAMPLE}.sam > ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam #keep the header -h

echo "filter them by length"
samtools view -h ${MAPPING_DIR}/${SAMPLE}_onlymapped.sam | awk 'length($10) > 80 || $1 ~ /^@/' > ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam

#if not done so already, index the fasta to be used with samclip
echo "index the file for the samclip"

samtools faidx human_mito_ref.fasta

#echo "then  filter them by hard clipping"  #adding this back to remove non specific hits to other worms
samclip --ref human_mito_ref.fasta.fai ${MAPPING_DIR}/${SAMPLE}_onlymapped_filtered.sam > ${MAPPING_DIR}/${SAMPLE}_samclip.sam

echo "convert to bam"
samtools view -S -b ${MAPPING_DIR}/${SAMPLE}_samclip.sam > ${MAPPING_DIR}/${SAMPLE}.bam
#sort
samtools sort -o ${MAPPING_DIR}/${SAMPLE}_sorted.bam ${MAPPING_DIR}/${SAMPLE}.bam
#remove duplicates

rm ${MAPPING_DIR}/${SAMPLE}.sam ${MAPPING_DIR}/${SAMPLE}_*.sam

sambamba view -t 12 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${MAPPING_DIR}/${SAMPLE}_sorted.bam -o ${MAPPING_DIR}/${SAMPLE}_uniq.bam

samtools sort -o ${MAPPING_DIR}/${SAMPLE}_uniq_sorted.bam ${MAPPING_DIR}/${SAMPLE}_uniq.bam

echo "tag the duplicate reads"
#tag the duplicates here with picard. It will give a warning/error but it does not stpo picard from running
picard -Xmx4G MarkDuplicates I=${MAPPING_DIR}/${SAMPLE}_uniq_sorted.bam REMOVE_DUPLICATES=TRUE O=${MAPPING_DIR}/${SAMPLE}_uniq_sorted_dup.bam M=${MAPPING_DIR}/${SAMPLE}_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT
#you might need to add -Xmx2G or -Xmx4G with picard, as now the files are larger (mapped to the nuclear/whole genomes),so picard needs more memory now

samtools index  ${MAPPING_DIR}/${SAMPLE}_uniq_sorted_dup.bam  #you need to sort before you index and you need the sorting before calling duplicates too

samtools view -h ${MAPPING_DIR}/${SAMPLE}_uniq_sorted_dup.bam | awk '$6 !~ /H|S/ || $1 ~ /@/' | samtools view -bS - > ${MAPPING_DIR}/${SAMPLE}_filtered_CIGAR_final.bam

samtools index ${MAPPING_DIR}/${SAMPLE}_filtered_CIGAR_final.bam #bedtools multicov needs indexing to calculate stats

echo "calculate some basic stats"

#samtools idxstats ${MAPPING_DIR}/${SAMPLE}_uniq_sorted_dup.bam > ${MAPPING_DIR}/${SAMPLE}_uniq_sorted_dup.txt #or the post CIGAR filtered bam

samtools idxstats ${MAPPING_DIR}/${SAMPLE}_filtered_CIGAR_final.bam > ${MAPPING_DIR}/${SAMPLE}_final.txt

#calculate depth per base - no longer needed
#samtools depth ${MAPPING_DIR}/${SAMPLE}_uniq_sorted.bam > ${MAPPING_DIR}/deduped_${SAMPLE}.coverage #don't need this any more, as we are not working with conse>

mv *.err  ${MAPPING_DIR}/ERR_FILES
mv *.out ${MAPPING_DIR}/OUT_FILES

echo  "I am DONE !!!!!!!"
#------------------------------------------------------------------------------

```