
##########################################
##### NATURE COMMS - REVIEWER COMMENTS 
##########################################
#Assessing reviewer's comments on the filters of coverage selected for Grenedalf in the pools 
#####11 April 2025

conda activate grenedalf

/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING/04_ASCARIS_SPP_MAPPING/03_ASCARIS_SUUM_82_SAMPLES_VARS/01_Pi_diversity_output/04_COVERAGE_TESTS_FOR_NATCOMMS


#Tested coverage of 1, 5, 10, 20 
###############################
### Pooled Pi diversity ----
###############################

##### SCRIPT TO RUN THE CONDITIONS 
#!/bin/bash
# Running a loop to generate pi diversity per bam file per 'kept sites' file,
# with varying filter-sample-min-coverage values: 1, 5, 10, 20
# 17 Oct 2023 (modified version)

for k in *HIGHQUAL*.kept.sites # Steve and I decided to take forward only the HIGHQUAL sites
do
for i in *_filtered_CIGAR_final.bam
do
nickname=$(basename "$k" | awk -F '_' '{print $1"_"$4}')
name=$(basename "$i" | sed 's/_trimmed.*//')

for coverage in 1 5 10 20
do
/home/marip3/apps/grenedalf/bin/grenedalf diversity \
--sam-path "$i" \
--filter-sample-min-count 2 \
--filter-sample-min-coverage "$coverage" \
--filter-region-list "$k" \
--window-type chromosomes \
--pool-sizes 1000 \
--file-prefix "${nickname}_${name}_pi_cov${coverage}_"
done
done
done
#Pi_GRENEDALF_BAMS_FILTER_COVERAGE.sh (END)

###############
#Where are you on the server: 
#we are here: 
/home/marip3/mbl_genome_skimming/03.GLOBAL_SKIM/04.ANALYSIS/01.MTDNA_MAPPING/04_ASCARIS_SPP_MAPPING/03_ASCARIS_SUUM_82_SAMPLES_VARS/01_Pi_diversity_output/04_COVERAGE_TESTS_FOR_NATCOMMS

#and on home directory you are here: 
/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/07.Global_genome_skim_2023/02_DATA/02_TRIMMED_DATA/04_VARIANT_CALLING/01_MITOGENOME_VARS/09_ASCARIS_SPP_VARS/03_ASUUM_US_82SAMPLES/01_Pi_diversity_output/01_CONDITIONAL_COVERAGE_NATCOMMS/

#need to merge the pi files now, and add 'condition' 'column based on the name of the file and then plot them 

library(tidyverse)

setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/07.Global_genome_skim_2023/02_DATA/02_TRIMMED_DATA/04_VARIANT_CALLING/01_MITOGENOME_VARS/09_ASCARIS_SPP_VARS/03_ASUUM_US_82SAMPLES/01_Pi_diversity_output/01_CONDITIONAL_COVERAGE_NATCOMMS/")

# Set your CSV directory
ASUUM_WITH_COVERAGE_CONDITIONS_csv_dir <- "/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/07.Global_genome_skim_2023/02_DATA/02_TRIMMED_DATA/04_VARIANT_CALLING/01_MITOGENOME_VARS/09_ASCARIS_SPP_VARS/03_ASUUM_US_82SAMPLES/01_Pi_diversity_output/01_CONDITIONAL_COVERAGE_NATCOMMS/"

# List CSV files
ASUUM_WITH_COVERAGE_CONDITIONS_csv_files <- list.files(ASUUM_WITH_COVERAGE_CONDITIONS_csv_dir, pattern = ".*_pi_cov\\d+_diversity\\.csv", full.names = TRUE)

# Define country map
country_map <- c(
  "BEN" = "Benin","CMR" = "Cameroon", "TZA" = "Tanzania","ETH" = "Ethiopia",
  "UGA" = "Uganda", "IND" = "India", "NHD" = "Honduras","MWI" = "Malawi",
  "MMR" = "Myanmar","NGA" = "Nigeria","ARG" = "Argentina", "BGD" = "Bangladesh",
  "CHN" = "China","ECU" = "Ecuador", "GLP" = "Guadeloupe","FJI" = "Fiji",
  "LKA" = "Sri Lanka","MOZ" = "Mozambique", "ITA" = "Italy", "PR" = "Puerto Rico",
  "COD" = "DRC","SEN" = "Senegal","THA" = "Thailand","USA" = "U.S.A",
  "ZAF" = "South Africa","MYS" = "Malaysia", "KEN" = "Kenya"
)

# Init list
all_data <- list()

# Loop through each file
for (ASUUM_WITH_COVERAGE_CONDITIONS_csv_file in ASUUM_WITH_COVERAGE_CONDITIONS_csv_files) {
  
  filename <- basename(ASUUM_WITH_COVERAGE_CONDITIONS_csv_file)
  
  # Extract sample ID properly (between known prefix and "_pi")
  sample_id <- str_extract(filename, "(?<=ASUUM_bcftools_)(.*?)(?=_pi)")
  
  if (is.na(sample_id)) {
    warning(paste("Could not extract sample_id from:", filename))
    next
  }
  
  # Country code = first 3 letters of sample_id
  country_code <- substr(sample_id, 1, 3)
  
  # Extract coverage
  coverage <- str_extract(filename, "(?<=_cov)\\d+")
  
  # Read file
  df <- read_csv(ASUUM_WITH_COVERAGE_CONDITIONS_csv_file, show_col_types = FALSE)
  
  # Add metadata
  df$sample.id <- sample_id
  df$country <- country_map[country_code]
  df$coverage <- as.numeric(coverage)
  df$sample_type <- "pools"
  
  # Rename first 10 columns
  colnames(df)[1:10] <- c(
    "species", "start", "end", "snp_count", "coverage_fraction",
    "theta_pi_abs", "theta_pi_rel", "theta_watterson_abs",
    "theta_watterson_rel", "tajimas_d"
  )
  
  # Add to list
  all_data[[length(all_data) + 1]] <- df
}

# Combine all into one dataframe
final_df <- bind_rows(all_data)

#isolate only coverage 1 from the above to match it to the mean depth below 
FINAL_DF_Pi_simplified <- final_df %>%
  select(1,4,7, 11,12,13) %>%
  filter(coverage=="1")



library(ggplot2)

# Basic boxplot of theta_pi_rel per coverage level
ggplot(final_df, aes(x = factor(coverage), y = theta_pi_rel)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  geom_jitter(width=0.1)+
  labs(
    title = "Theta_pi (rel) per Coverage Level \n 16 samples (BAM) \n Grenedalf",
    x = "Coverage",
    y = expression(theta[pi]~"(relative)")
  ) +
  theme_bw()
ggsave("XX_SUPPLEMENTARY_FIGX_POPGEN_Nucleotide_diversity_REVIEWER_1_COVERAGE_FILTERS.pdf", width=400, height=250, units="mm")


############################################################################
### MEAN SAMPLE DEPTH AND THETA RELATIVE Pi PER SAMPLE -----
############################################################################
#18 April 2025

#CALCULATING MEAN DEPTH FOR THE BAM FILES (ASCARIS SUUM) FOR NATCOMMS 
##############################################################################################
#!/bin/bash

# Output header
echo -e "sample_id\tmean_depth\tchrom" > ASUUM_16_SAMPLES_collated_MEAN_depth_results.tsv

# Target chromosome
chrom="NC_001327_Ascaris_suum_mitochondrion_genome_USA"
species="Ascaris_suum"

# Loop through each BAM file
for bam in *.bam; do
sample=$(basename "$bam" .bam)

# Create a depth file for each sample (optional to keep)
samtools depth -a "$bam" > "${sample}_depth.txt"

# Calculate mean depth for the specific chromosome
mean_depth=$(awk -v chr="$chrom" '$1 == chr { sum += $3; count++ } END { print (count > 0 ? sum/count : 0) }' "${sample}_depth.txt")

# Append results to output table
echo -e "$sample\t$mean_depth\t$species" >> ASUUM_16_SAMPLES_collated_MEAN_depth_results.tsv
done

#############################################################################################

#import the mean depth 
setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/07.Global_genome_skim_2023/02_DATA/02_TRIMMED_DATA/04_VARIANT_CALLING/01_MITOGENOME_VARS/09_ASCARIS_SPP_VARS/03_ASUUM_US_82SAMPLES/01_Pi_diversity_output/02_MEAN_DEPTH_NATCOMMS/")

MEAN_DEPTH <- read.table("ASUUM_16_SAMPLES_collated_MEAN_depth_results.tsv", header =T)

MEAN_DEPTH$sample_id <- sub("_trimmed.*", "", MEAN_DEPTH$sample_id)

colnames(MEAN_DEPTH) <- c('sample.id', 'mean_depth','chrom')

#now need to combine the two 
COMBINED_Pi_MEAN_DEPTH <- merge(MEAN_DEPTH, FINAL_DF_Pi_simplified, by = "sample.id")

#pLot
# Plot
ggplot(COMBINED_Pi_MEAN_DEPTH, aes(x = mean_depth, y = theta_pi_rel, color = country)) +
  geom_point(size = 4) +  # Use points to represent each data entry
  scale_x_log10() +  # Log scale for x-axis (mean_depth)
  #scale_y_log10() +  # Log scale for y-axis (theta_pi_rel)
  labs(
    title = "Mean Depth vs Theta Pi Relative by Country",
    x = "Mean Depth",
    y = "Theta Pi Relative",
    color = "Country"
  ) +
  theme_bw() +  # A clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity
ggsave("XX_SUPPLEMENTARY_FIGX_POPGEN_Nucleotide_diversity_Vs_Mean_Coverage_REVIEWER_1.pdf", width=200, height=150, units="mm")



#####################

#Tested coverage of 1, 5, 10, 20 
###############################
### Pooled Dxy  ----
###############################

#Using Grenedalf again, the code will be like: 
/home/marip3/apps/grenedalf/bin/grenedalf fst --write-pi-tables --sam-path . --filter-sample-min-count 2 --filter-sample-min-coverage 1 --filter-region-list ALUM_bcftools_mtDNA_HIGHQUAL_n500_FORMAT_AD_removeddups.recode_MAX_MISS_0.7_pools.kept.sites --window-type chromosomes --pool-sizes 1000 --method unbiased-hudson --file-prefix ALUM_pi_table_hudson_

#where: 
#--filter-sample-min-coverage 1, I will change to 1, 5, 10, 20 like I did with the pools and then take a single table 
#also added 40 and 80 to see if the number of SNPs does change after some point. 
#the above gives me number of SNPs so I can calculate: 
#(i) does the number of SNPs change? No until coverage of 20 
#(ii) does the between pi change? (maybe plot it as a heatmap like I did before or in boxplots)

#####################################################################
#!/bin/bash
# Loop to run fst with different filter-sample-min-coverage values

for coverage in 1 5 10 20 40 80
do
/home/marip3/apps/grenedalf/bin/grenedalf fst \
--write-pi-tables \
--sam-path . \
--filter-sample-min-count 2 \
--filter-sample-min-coverage "$coverage" \
--filter-region-list ASUUM_US_82SAMPLES_bcftools_mtDNA_HIGHQUAL_n500_NOMINIMUMALLELEFREQ_nodups.recode_MAX_MISS_0.7_w_filtered_Indiv_pools.kept.sites \
--window-type chromosomes \
--pool-sizes 1000 \
--method unbiased-hudson \
--file-prefix cov${coverage}_ALUM_pi_table_hudson_
done
#Between_pi_Dxy_GRENEDALF_BAMS.sh (END)
################################################################


setwd("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/07.Global_genome_skim_2023/02_DATA/02_TRIMMED_DATA/04_VARIANT_CALLING/01_MITOGENOME_VARS/09_ASCARIS_SPP_VARS/03_ASUUM_US_82SAMPLES/04_Grenedalf_Dxy_output/01_Dxy_CONDITIONAL_COVERAGE_NATCOMMS/")

dxy_data_transformation <- function(file) {
  pi_data2 <- read.table(file, sep=",", header=T)
  columns <- colnames(pi_data2)
  # Use the sub function to replace the specified substring
  new_columns <- gsub("_filtered_CIGAR_final.1", "", columns)
  # Rename the columns in your data frame
  colnames(pi_data2) <- new_columns #this works!
  
  #turn it into long format
  pi_data2_l <- pi_data2 %>% 
    pivot_longer(
      cols = 5:ncol(pi_data2), 
      names_to = "sample_combos",
      values_to = "dxy" )
  #split sample combos
  pi_data2_l_split <- pi_data2_l %>% separate_wider_delim(col = sample_combos, ".", names = c("sample_1", "sample_2")) #works 
  
  
  country_map <- c("BEN" = "Benin","CMR" = "Cameroon", "TZA" = "Tanzania","ETH" = "Ethiopia","UGA" = "Uganda", "IND" = "India", "HND" ="Honduras",
                   "MWI" = "Malawi","MMR" = "Myanmar","NGA" = "Nigeria","ARG" = "Argentina", "BGD" = "Bangladesh","CHN" = "China",
                   "ECU" = "Ecuador", "GLP" = "Guadeloupe","FJI" = "Fiji","LKA" = "Sri Lanka","MOZ" = "Mozambique", "ITA" = "Italy",
                   "PR" = "Puerto Rico","COD" = "DRC","SEN" = "Senegal","THA" = "Thailand","USA" = "U.S.A","ZAF" = "South Africa",
                   "MYS" = "Malaysia", "KEN" ="Kenya")
  
  pi_data2_l_split_country <- pi_data2_l_split %>%
    mutate(country_1 = country_map[substr(sample_1, 1, 3)],
           country_2=country_map[substr(sample_2, 1,3)]) #YAAAAAAS THAT WORKS !!!!!
  
  
  #write_csv(xx, "xxx.csv")
  #Calculate now mean Dxy per country?
  pi_data2_l_split_country_mean_dxy <- pi_data2_l_split_country %>%
    group_by(country_1, country_2) %>%
    mutate(mean_dxy_country = mean(dxy))
  
}

#Call the function
colnames(cov1_dxy_pi_data_ASUUM_82SAMPLES_pools)
cov1_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov1_ALUM_pi_table_hudson_pi-between.csv")
cov1_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "1"

cov5_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov5_ALUM_pi_table_hudson_pi-between.csv")
cov5_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "5"

cov10_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov10_ALUM_pi_table_hudson_pi-between.csv")
cov10_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "10"

cov20_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov20_ALUM_pi_table_hudson_pi-between.csv")
cov20_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "20"

cov40_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov40_ALUM_pi_table_hudson_pi-between.csv")
cov40_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "40"

cov80_dxy_pi_data_ASUUM_82SAMPLES_pools <- dxy_data_transformation("cov80_ALUM_pi_table_hudson_pi-between.csv")
cov80_dxy_pi_data_ASUUM_82SAMPLES_pools$coverage <- "80"


#keeping unique country pairs and values per dataset 
library(dplyr)

# Get all the relevant data frame names from the global environment
dxy_datasets <- ls(pattern = "^cov\\d+_dxy_pi_data_ASUUM_82SAMPLES_pools$")

# Loop through each and apply the transformation
for (df_name in dxy_datasets) {
  # Get the dataframe
  df <- get(df_name)
  
  # Apply the logic: keep only one row per unordered country pair
  df_clean <- df %>%
    rowwise() %>%
    mutate(
      country_min = min(country_1, country_2),
      country_max = max(country_1, country_2)
    ) %>%
    ungroup() %>%
    distinct(country_min, country_max, .keep_all = TRUE) %>%
    select(chrom, snps, country_1, country_2, mean_dxy_country, coverage)
  
  # Save cleaned version back with a "_unique" suffix
  assign(paste0(df_name, "_unique"), df_clean)
}


#now combine all the datasets and plot them/facet them by coverage

library(dplyr)
library(ggplot2)

# Combine all cleaned *_unique dataframes into one
dxy_heatmap_data <- bind_rows(
  cov1_dxy_pi_data_ASUUM_82SAMPLES_pools_unique,
  cov5_dxy_pi_data_ASUUM_82SAMPLES_pools_unique,
  cov10_dxy_pi_data_ASUUM_82SAMPLES_pools_unique,
  cov20_dxy_pi_data_ASUUM_82SAMPLES_pools_unique,
  cov40_dxy_pi_data_ASUUM_82SAMPLES_pools_unique,
  cov80_dxy_pi_data_ASUUM_82SAMPLES_pools_unique
)

# Create the base row for Mozambique with coverage set to "1"
Mozambique_base <- data.frame(
  chrom = "NC_001327_Ascaris_suum_mitochondrion_genome_USA",
  snps = "308",
  country_1 = "Mozambique",
  country_2 = "Mozambique",
  mean_dxy_country = NA, 
  coverage = "1"
)

# Create a vector of coverage values
coverage_values <- c(1, 5, 10, 20, 40, 80)

# Use bind_rows to create identical rows with different coverage values
Mozambique_expanded <- do.call(rbind, lapply(coverage_values, function(cov) {
  Mozambique_base$coverage <- as.character(cov)
  return(Mozambique_base)
}))

#bind now Mozambique with the rest 
dxy_heatmap_data_MOZ <- rbind (dxy_heatmap_data, Mozambique_expanded)

#doing some magic so the heatmap is symmetrical 

ASUUM_ALL_unique_countries_1_pools_coverage_filters <- dxy_heatmap_data_MOZ  %>%
  distinct(country_1) %>%
  pull()

ASUUM_ALL_unique_countries_2_pools_coverage_filters <- dxy_heatmap_data_MOZ %>%
  distinct(country_2) %>%
  pull()
ASUUM_ALL_unique_countries_pools_coverage_filters <- union(ASUUM_ALL_unique_countries_1_pools_coverage_filters, ASUUM_ALL_unique_countries_2_pools_coverage_filters)



##PLOT
#ggplot(
#  dxy_heatmap_data_MOZ,
#  aes(
#    x = factor(country_2, levels = ASUUM_ALL_unique_countries_pools_coverage_filters),
#    y = factor(country_1, levels = ASUUM_ALL_unique_countries_pools_coverage_filters),
#    fill = mean_dxy_country,
#    label = round(mean_dxy_country, 4)
#  )
#) +
#  geom_tile(color = "white") +
#  geom_text(
#    aes(label = ifelse(is.na(mean_dxy_country), "NA", round(mean_dxy_country, 4))),
#    color = "white",
#    size = 2.5
#  ) +
#  scale_fill_gradient(
#    low = "purple4",
#    high = "yellow2",
#    limits = c(0, 0.5),
#    name = "Mean DXY",
#    na.value = "grey"
#  ) +
# # theme_minimal(base_size = 13) +
#  labs(
#    title = "Mean DXY per Country Pair (ASUUM Pools \n 16 samples based on BAM \n Grenedalf)",
#    x = "Country 2",
#    y = "Country 1"
#  ) +
#  facet_wrap(~ coverage, ncol = 3) +
#  theme(
#    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
# 
#   # strip.background = element_rect(fill = "grey90")
#  )
#ggsave("XX_SUPPLEMENTARY_FIGX_POPGEN_DXY_REVIEWER_1_COVERAGE_FILTERS.pdf", width=400, height=250, units="mm")
#
#iF ONLY WANTING THE 1,5, 10, 20 to be consistent with the pi diversity 

# Define which coverages to keep and how to label them
selected_coverages <- c("1", "5", "10", "20")
coverage_labels <- setNames(
  paste("Coverage of", selected_coverages),
  selected_coverages
)

# Filter data
filtered_heatmap_data <- dxy_heatmap_data_MOZ %>%
  filter(coverage %in% selected_coverages)

#PLOT 
filtered_heatmap_data$coverage <- factor(
  filtered_heatmap_data$coverage,
  levels = c("1", "5", "10", "20", "40", "80"),
  labels = c("Coverage of 1", "Coverage of 5", "Coverage of 10", 
             "Coverage of 20", "Coverage of 40", "Coverage of 80")
)

#PLOT
ggplot(
  filtered_heatmap_data,
  aes(
    x = factor(country_2, levels = ASUUM_ALL_unique_countries_pools_coverage_filters),
    y = factor(country_1, levels = ASUUM_ALL_unique_countries_pools_coverage_filters),
    fill = mean_dxy_country,
    label = round(mean_dxy_country, 4)
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(is.na(mean_dxy_country), "NA", round(mean_dxy_country, 4))),
    color = "white",
    size = 2.5
  ) +
  scale_fill_gradient(
    low = "purple4",
    high = "yellow2",
    limits = c(0, 0.5),
    name = "Mean DXY",
    na.value = "grey"
  ) +
  # theme_minimal(base_size = 13) +
  labs(
    title = "Mean DXY per Country Pair (ASUUM Pools \n 16 samples based on BAM \n Grenedalf)",
    x = "Country 2",
    y = "Country 1"
  ) +
  facet_wrap(~ coverage, ncol = 3) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    
    # strip.background = element_rect(fill = "grey90")
  )
ggsave("XX_SUPPLEMENTARY_FIGX_POPGEN_DXY_REVIEWER_1_COVERAGE_FILTERS_ONLY_1_5_10_20.pdf", width=400, height=250, units="mm")
