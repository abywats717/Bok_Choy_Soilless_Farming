# Bok choy hydroponic study
#Differential abundance analysis at ASV level
#Auja Bywater
#Last updated: 1/23/26

#Load packages
library(ggplot2)
library(BiocManager)
#remotes::install_git("https://git.bioconductor.org/packages/ALDEx2")
library(ALDEx2, lib.loc = "/storage/home/akb6471/R/")

library(viridis)
library(dendextend)
library(readxl)
library(gplots)
library(BiocParallel)
library(phyloseq)
library(ape)
library(compositions)
library(dplyr)
library(tidyr)
library(psych)
library(cowplot)
library(reshape2)
library(svglite)
library(tibble)
library(stringr)
library(multtest)
library(ggplot2)
library(forcats)
library(randomcoloR)
library(stringr)

#Set working directory
setwd('/storage/home/akb6471/work/hydroponics_rerun/Reads_forSRA/')

# Import Data

# Load in filtered raw 16s data
asv_16s <- read.csv("asv_16s_filtered.csv", header=TRUE)
colnames(asv_16s)[colnames(asv_16s) == "X"] <- "SampleID"

taxon_16s <- as.matrix(read.csv("Taxon_16s_clean.csv", header=TRUE, row.names=1))

metadata_16s <- read.csv('metadata_16s_updated99.csv', header=TRUE, row.names=1)
colnames(metadata_16s)[colnames(metadata_16s) == "X"] <- "SampleID"

# ALDeX2 can only be used for two groups, and I have five systems I need to compare. 
#I'm creating a for loop to handle pairwise comparisons between each system

#Create a list of the hydroponic systems to compare
unique_groups <- c("di", "dwc", "ef", "kr", "nf")

# The ASV table needs to have samples as columns, not rows
asv_16s_T_transposed <- t(asv_16s)
colnames(asv_16s_T_transposed) <- metadata_16s$SampleID
asv_16s_T_transposed <- asv_16s_T_transposed[-1, ]

# Create a list for the results
results_list <- list()

# Pairwise comparisons
for(i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    #Define the two groups for comparison
    group1 <- unique_groups[i]
    group2 <- unique_groups[j]
    
    # Print which groups are being compared
    print(paste("Comparing:", group1, "vs", group2))
    
    #Subset metadata and ASV table for the two groups
    selected_samples <- metadata_16s$system %in% c(group1, group2)
    metadata_subset <- metadata_16s[selected_samples, ]
    metadata_subset$system <- as.factor(metadata_subset$system)
    asv_subset <- asv_16s_T_transposed[, selected_samples]
    
    # Convert asv_subset to a data frame and then to numeric
    asv_subset_df <- as.data.frame(asv_subset)
    asv_subset_df[] <- lapply(asv_subset_df, as.numeric)
    asv_subset <- as.matrix(asv_subset_df)
    
    # Run ALDEx2 for the two groups
    aldex_clr <- aldex.clr(asv_subset, mc.samples = 128, conds = metadata_subset$system)
    aldex_text <- aldex.ttest(aldex_clr, paired.test = FALSE)
    
    #Add group names to the results
    aldex_text$group1 <- group1
    aldex_text$group2 <- group2
    
    # Store the results in the results list
    results_list[[paste(group1, "_vs_", group2)]] <- aldex_text
     }
}

# Combine all results into a single data frame
final_results <- do.call(rbind, results_list)

# Filter for significant interactions based on the p-value
significant_results <- final_results[final_results$we.eBH < 0.05, ]
   
# Restructure significant_results so the and system interactions asvs are in their own columns
library(dplyr)
library(stringr)

# Convert the first row to a proper column
significant_results <- significant_results %>%
  rownames_to_column(var = "ASV")

# Now extract the ASV identifier
significant_results <- significant_results %>%
  mutate(ASV = str_extract(ASV, "(?<=\\.)ASV\\d+"))

# Combine group1 and group2
significant_results <- significant_results %>% 
  mutate(systems = paste(group1, group2, sep = "-"))

# Make relative abundance plots
# Make the column (SampleId) rownames
asv_16s_r <- rownames(asv_16s) <- asv_16s[[1]]
asv_16s_r <- asv_16s[, -1]
asv_16s_r[] <- lapply(asv_16s_r, function(x) as.numeric(as.character(x)))

# Convert to relative abundance (per sample)
asv_16s_rel <- sweep(asv_16s_r, 1, rowSums(asv_16s_r), FUN = "/")
# Convert all columns of asv_16s_rel to numeric
asv_16s_rel[] <- lapply(asv_16s_rel, as.numeric)

# Save the relative abundance data
write.csv(asv_16s_rel, file="asv_16s_rel.csv")

# Add SampleID column from rownames
asv_16s_rel$SampleID <- rownames(asv_16s_rel)

# Reshape to long format
asv_16s_long <- pivot_longer(
  asv_16s_rel,
  cols = -SampleID,
  names_to = "ASV",
  values_to = "Abundance"
)

# Extract System
asv_16s_long$System <- str_extract(asv_16s_long$SampleID, "(?<=_)(di|dwc|dw|ef|kr|nft|nf)(?=_)")

# Replace 'dw' with 'dwc' only in the rows where 'dw' appears
asv_16s_long$System[asv_16s_long$System == "dw"] <- "dwc"
# Replace 'nf' with 'nft' only in the rows where 'nf' appears
asv_16s_long$System[asv_16s_long$System == "nf"] <- "nft"

# Extract Season
asv_16s_long$Season <- str_extract(asv_16s_long$SampleID, "_f_|_s_") %>%
  gsub("_", "", .)

# Merge taxonomy into the abundance data
taxon_16s_df <- as.data.frame(taxon_16s)
taxon_16s_df$ASV <- rownames(taxon_16s_df)
taxon_16s_df <- taxon_16s_df[, c(ncol(taxon_16s_df), 1:(ncol(taxon_16s_df) - 1))]
asv_16s_long <- left_join(asv_16s_long, taxon_16s_df, by = "ASV")

# STEP 1: Filter out low-abundance genera (>1%)
asv_16s_long_filtered <- asv_16s_long %>%
  filter(Abundance > 0.025)

# STEP 2: Calculate total observed abundance per sample
observed_abundance <- asv_16s_long_filtered %>%
  group_by(SampleID) %>%
  summarise(observed = sum(Abundance), .groups = "drop")

# STEP 3: Create the "Other" category
other_data <- observed_abundance %>%
  mutate(Abundance = 1 - observed) %>%
  filter(Abundance > 0) %>%
  transmute(SampleID, Abundance, Genus = "Other")

# STEP 4: Add System and Season metadata back to "Other"
other_data <- left_join(
  other_data,
  asv_16s_long_filtered %>% dplyr::select(SampleID, System, Season) %>% distinct(),
  by = "SampleID"
)

# STEP 5: Combine "Other" with filtered data
asv_16s_long_filled <- bind_rows(asv_16s_long_filtered, other_data)

# STEP 6: Generate distinct colors for genera (including "Other")
unique_genera <- unique(asv_16s_long_filled$Genus)
genus_colors <- distinctColorPalette(length(unique_genera))
names(genus_colors) <- unique_genera
genus_colors["Other"] <- "#CCCCCC"  # Optional: set "Other" to grey

# STEP 7: Create and save plots per system
systems <- unique(asv_16s_long_filled$System)

for (sys in systems) {
  
  # Filter for current system
  system_data <- asv_16s_long_filled %>% filter(System == sys)
  
  # Extract numeric day from SampleID (after the last underscore)
  system_data$Day <- as.numeric(str_extract(system_data$SampleID, "[^_]+$"))
  
  # Reorder SampleID and Genus
  system_data <- system_data %>%
    arrange(Season, Day) %>%
    mutate(
      SampleID = factor(SampleID, levels = unique(SampleID)),
      Genus = factor(Genus, levels = c("Other", setdiff(unique(Genus), "Other")))
    )
  
  # Plot
  p <- ggplot(system_data, aes(x = SampleID, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = genus_colors) +
    labs(
      title = paste("Relative Abundance in", toupper(sys), "System"),
      x = "Sample ID",
      y = "Relative Abundance",
      fill = "Genus"
    ) +
    facet_wrap(~Season, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
      plot.title = element_text(hjust = 0.5)
    )
  
  print(p) 
  
  # Optional: save to file
  ggsave(
    filename = paste0("relative_abundance_", sys, ".png"),
    plot = p,
    width = 10, height = 6, dpi = 300
  )
}
