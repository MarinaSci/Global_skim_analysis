# Trimming of raw sequencing reads 
Author: Marina Papaiakovou, mpapaiakovou[at]gmail.com 

## Contents: 
- Trimming of raw sequencing data on an HPC 

```bash

#TRIMMOMATIC ON HPC -----
#global skim update: 4 samples, did not have the Phred33 scores, but Phred64, so I had to add the --phred64 flag on trimmomatic 

#CODE: 
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_trim_%j_%a.out
#SBATCH --error=job_trim_%j_%a.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-4

export PATH=/home/marip3/miniconda3/bin/:$PATH
source activate mapping-env

#make a list of samples present
ls -1 *_1.fq.gz | sed 's/_1.fq.gz//' > sample_list.txt

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_list.txt)

echo "I will now run trimmomatic"

trimmomatic PE -phred33 -threads 8 ${SAMPLE}_1.fq.gz  ${SAMPLE}_2.fq.gz \
${SAMPLE}_trimmed_1.fq.gz ${SAMPLE}_unpaired_1.fq.gz ${SAMPLE}_trimmed_2.fq.gz \
${SAMPLE}_unpaired_2.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#adapters file contains ALL adapters out there 

echo "I AM DONE TRIMMING!"
#test

```