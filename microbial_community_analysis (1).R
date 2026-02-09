# Bok choy hydroponic study
#Microbial community analysis at genus level
#Auja Bywater
#Last updated: 8/26/25
#Based on: https://github.com/LauRolon/Apple-Year2/blob/main/Code/Core%20and%20network%20FINAL.R

#Set working directory
setwd('/storage/home/akb6471/work/hydroponics_rerun/Reads_forSRA/')

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
library(phyloseq)
library(vegan)
library(minpack.lm)
library(Hmisc)
library(stats4)
library(zCompositions)
library(compositions)
library(vegan)

# Load in data
asv_16s <- read.csv("asv_16s_filtered.csv", header=TRUE)
colnames(asv_16s)[1] <- "SampleID"
rownames(asv_16s) <- asv_16s$SampleID
asv_16s <- asv_16s[, -which(colnames(asv_16s) == "SampleID")]
taxon_16s <- as.matrix(read.csv("Taxon_16s_clean.csv", header=TRUE, row.names=1))

# Metadata 
metadata_16s <- read.csv("metadata_16s_updated99.csv", header = TRUE)
# Set the first column as rownames
rownames(metadata_16s) <- metadata_16s[[1]]
# Remove the first column (optional, but avoids duplication)
metadata_16s <- metadata_16s[, -1]

#### Prepare data for analyses ####
#Make phyloseq
phyloseq16s <- phyloseq(otu_table(asv_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))

# Subset phyloseq by season and system
physeq_16s_spring_dwc <- subset_samples(phyloseq16s, season == "s" & system == "dwc")
physeq_16s_fall_dwc <- subset_samples(phyloseq16s, season == "f" & system == "dwc")
physeq_16s_spring_kr <- subset_samples(phyloseq16s, season == "s" & system == "kr")
physeq_16s_fall_kr <- subset_samples(phyloseq16s, season == "f" & system == "kr")
physeq_16s_spring_nft <- subset_samples(phyloseq16s, season == "s" & system == "nf")
physeq_16s_fall_nft <- subset_samples(phyloseq16s, season == "f" & system == "nf")
physeq_16s_spring_ef <- subset_samples(phyloseq16s, season == "s" & system == "ef")
physeq_16s_fall_ef <- subset_samples(phyloseq16s, season == "f" & system == "ef")
physeq_16s_spring_di <- subset_samples(phyloseq16s, season == "s" & system == "di")
physeq_16s_fall_di <- subset_samples(phyloseq16s, season == "f" & system == "di")

# Extract ASV tables for each subset
asv_16s_spring_dwc <- t(as.data.frame(as(otu_table(physeq_16s_spring_dwc), "matrix")))
asv_16s_fall_dwc <- t(as.data.frame(as(otu_table(physeq_16s_fall_dwc), "matrix")))
asv_16s_spring_kr <- t(as.data.frame(as(otu_table(physeq_16s_spring_kr), "matrix")))
asv_16s_fall_kr <- t(as.data.frame(as(otu_table(physeq_16s_fall_kr), "matrix")))
asv_16s_spring_nft <- t(as.data.frame(as(otu_table(physeq_16s_spring_nft), "matrix")))
asv_16s_fall_nft <- t(as.data.frame(as(otu_table(physeq_16s_fall_nft), "matrix")))
asv_16s_spring_ef <- t(as.data.frame(as(otu_table(physeq_16s_spring_ef), "matrix")))
asv_16s_fall_ef <- t(as.data.frame(as(otu_table(physeq_16s_fall_ef), "matrix")))
asv_16s_spring_di <- t(as.data.frame(as(otu_table(physeq_16s_spring_di), "matrix")))
asv_16s_fall_di <- t(as.data.frame(as(otu_table(physeq_16s_fall_di), "matrix")))

# Extract metadata for each subset
metadata_16s_spring_dwc <- subset(metadata_16s, season == "s" & system == "dwc")
metadata_16s_fall_dwc <- subset(metadata_16s, season == "f" & system == "dwc")
metadata_16s_spring_kr <- subset(metadata_16s, season == "s" & system == "kr")
metadata_16s_fall_kr <- subset(metadata_16s, season == "f" & system == "kr")
metadata_16s_spring_nft <- subset(metadata_16s, season == "s" & system == "nf")
metadata_16s_fall_nft <- subset(metadata_16s, season == "f" & system == "nf")
metadata_16s_spring_ef <- subset(metadata_16s, season == "s" & system == "ef")
metadata_16s_fall_ef <- subset(metadata_16s, season == "f" & system == "ef")
metadata_16s_spring_di <- subset(metadata_16s, season == "s" & system == "di")
metadata_16s_fall_di <- subset(metadata_16s, season == "f" & system == "di")

# Remove ASVs that have count zero in all samples
asv_16s_spring_dwc <- asv_16s_spring_dwc[ which(rowSums(asv_16s_spring_dwc)>0),]
asv_16s_fall_dwc <- asv_16s_fall_dwc[ which(rowSums(asv_16s_fall_dwc)>0),]
asv_16s_spring_kr <- asv_16s_spring_kr[ which(rowSums(asv_16s_spring_kr)>0),]
asv_16s_fall_kr <- asv_16s_fall_kr[ which(rowSums(asv_16s_fall_kr)>0),]
asv_16s_spring_nft <- asv_16s_spring_nft[ which(rowSums(asv_16s_spring_nft)>0),]
asv_16s_fall_nft <- asv_16s_fall_nft[ which(rowSums(asv_16s_fall_nft)>0),]
asv_16s_spring_ef <- asv_16s_spring_ef[ which(rowSums(asv_16s_spring_ef)>0),]
asv_16s_fall_ef <- asv_16s_fall_ef[ which(rowSums(asv_16s_fall_ef)>0),]
asv_16s_spring_di <- asv_16s_spring_di[ which(rowSums(asv_16s_spring_di)>0),]
asv_16s_fall_di <- asv_16s_fall_di[ which(rowSums(asv_16s_fall_di)>0),]

#### Occupancy abundance curves for each subset ####
# Calculate presence absence for 16s data
asv_16s_spring_dwc_PA <- 1*((asv_16s_spring_dwc > 0) == 1)
asv_16s_fall_dwc_PA <- 1*((asv_16s_fall_dwc > 0) == 1)
asv_16s_spring_kr_PA <- 1*((asv_16s_spring_kr > 0) == 1)
asv_16s_fall_kr_PA <- 1*((asv_16s_fall_kr > 0) == 1)
asv_16s_spring_nft_PA <- 1*((asv_16s_spring_nft > 0) == 1)
asv_16s_fall_nft_PA <- 1*((asv_16s_fall_nft > 0) == 1)
asv_16s_spring_ef_PA <- 1*((asv_16s_spring_ef > 0) == 1)
asv_16s_fall_ef_PA <- 1*((asv_16s_fall_ef > 0) == 1)
asv_16s_spring_di_PA <- 1*((asv_16s_spring_di > 0) == 1)
asv_16s_fall_di_PA <- 1*((asv_16s_fall_di > 0) == 1)

# Calculate mean occupancy for each subset and rename columns 
asv_16s_spring_dwc_occ <- as.data.frame(rowSums(asv_16s_spring_dwc_PA) / ncol(asv_16s_spring_dwc_PA))
colnames(asv_16s_spring_dwc_occ) <- "MeanOcc"
asv_16s_fall_dwc_occ <- as.data.frame(rowSums(asv_16s_fall_dwc_PA) / ncol(asv_16s_fall_dwc_PA))
colnames(asv_16s_fall_dwc_occ) <- "MeanOcc"
asv_16s_spring_kr_occ <- as.data.frame(rowSums(asv_16s_spring_kr_PA) / ncol(asv_16s_spring_kr_PA))
colnames(asv_16s_spring_kr_occ) <- "MeanOcc"
asv_16s_fall_kr_occ <- as.data.frame(rowSums(asv_16s_fall_kr_PA) / ncol(asv_16s_fall_kr_PA))
colnames(asv_16s_fall_kr_occ) <- "MeanOcc"
asv_16s_spring_nft_occ <- as.data.frame(rowSums(asv_16s_spring_nft_PA) / ncol(asv_16s_spring_nft_PA))
colnames(asv_16s_spring_nft_occ) <- "MeanOcc"
asv_16s_fall_nft_occ <- as.data.frame(rowSums(asv_16s_fall_nft_PA) / ncol(asv_16s_fall_nft_PA))
colnames(asv_16s_fall_nft_occ) <- "MeanOcc"
asv_16s_spring_ef_occ <- as.data.frame(rowSums(asv_16s_spring_ef_PA) / ncol(asv_16s_spring_ef_PA))
colnames(asv_16s_spring_ef_occ) <- "MeanOcc"
asv_16s_fall_ef_occ <- as.data.frame(rowSums(asv_16s_fall_ef_PA) / ncol(asv_16s_fall_ef_PA))
colnames(asv_16s_fall_ef_occ) <- "MeanOcc"
asv_16s_spring_di_occ <- as.data.frame(rowSums(asv_16s_spring_di_PA) / ncol(asv_16s_spring_di_PA))
colnames(asv_16s_spring_di_occ) <- "MeanOcc"
asv_16s_fall_di_occ <- as.data.frame(rowSums(asv_16s_fall_di_PA) / ncol(asv_16s_fall_di_PA))
colnames(asv_16s_fall_di_occ) <- "MeanOcc"

# Calculate relative abundance for each ASV by System
# First replace zero values with a small non-zero value before clr transformation. 
# Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16s_spring_dwc <- t(cmultRepl(t(asv_16s_spring_dwc), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_fall_dwc <- t(cmultRepl(t(asv_16s_fall_dwc), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_spring_kr <- t(cmultRepl(t(asv_16s_spring_kr), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_fall_kr <- t(cmultRepl(t(asv_16s_fall_kr), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_spring_nft <- t(cmultRepl(t(asv_16s_spring_nft), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_fall_nft <- t(cmultRepl(t(asv_16s_fall_nft), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_spring_ef <- t(cmultRepl(t(asv_16s_spring_ef), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_fall_ef <- t(cmultRepl(t(asv_16s_fall_ef), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_spring_di <- t(cmultRepl(t(asv_16s_spring_di), label=0, method="CZM", output="p-counts", z.warning = 0.99))
asv.n0_16s_fall_di <- t(cmultRepl(t(asv_16s_fall_di), label=0, method="CZM", output="p-counts", z.warning = 0.99))

# Then convert to compositions (relative abundance) with Aitchinson
asv.n0.acomp_16s_spring_dwc <- as.data.frame(acomp(t(asv.n0_16s_spring_dwc)), total=1)
asv.n0.acomp_16s_fall_dwc <- as.data.frame(acomp(t(asv.n0_16s_fall_dwc)), total=1)
asv.n0.acomp_16s_spring_kr <- as.data.frame(acomp(t(asv.n0_16s_spring_kr)), total=1)
asv.n0.acomp_16s_fall_kr <- as.data.frame(acomp(t(asv.n0_16s_fall_kr)), total=1)
asv.n0.acomp_16s_spring_nft <- as.data.frame(acomp(t(asv.n0_16s_spring_nft)), total=1)
asv.n0.acomp_16s_fall_nft <- as.data.frame(acomp(t(asv.n0_16s_fall_nft)), total=1)
asv.n0.acomp_16s_spring_ef <- as.data.frame(acomp(t(asv.n0_16s_spring_ef)), total=1)
asv.n0.acomp_16s_fall_ef <- as.data.frame(acomp(t(asv.n0_16s_fall_ef)), total=1)
asv.n0.acomp_16s_spring_di <- as.data.frame(acomp(t(asv.n0_16s_spring_di)), total=1)
asv.n0.acomp_16s_fall_di <- as.data.frame(acomp(t(asv.n0_16s_fall_di)), total=1)

# Calculate the mean relative abundance for each ASV in each facility and rename column head
asv_16s_spring_dwc_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_spring_dwc)) / ncol(t(asv.n0.acomp_16s_spring_dwc)))
colnames(asv_16s_spring_dwc_abun) <- "MeanAbun"
asv_16s_fall_dwc_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_fall_dwc)) / ncol(t(asv.n0.acomp_16s_fall_dwc)))
colnames(asv_16s_fall_dwc_abun) <- "MeanAbun"
asv_16s_spring_kr_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_spring_kr)) / ncol(t(asv.n0.acomp_16s_spring_kr)))
colnames(asv_16s_spring_kr_abun) <- "MeanAbun"
asv_16s_fall_kr_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_fall_kr)) / ncol(t(asv.n0.acomp_16s_fall_kr)))
colnames(asv_16s_fall_kr_abun) <- "MeanAbun"
asv_16s_spring_nft_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_spring_nft)) / ncol(t(asv.n0.acomp_16s_spring_nft)))
colnames(asv_16s_spring_nft_abun) <- "MeanAbun"
asv_16s_fall_nft_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_fall_nft)) / ncol(t(asv.n0.acomp_16s_fall_nft)))
colnames(asv_16s_fall_nft_abun) <- "MeanAbun"
asv_16s_spring_ef_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_spring_ef)) / ncol(t(asv.n0.acomp_16s_spring_ef)))
colnames(asv_16s_spring_ef_abun) <- "MeanAbun"
asv_16s_fall_ef_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_fall_ef)) / ncol(t(asv.n0.acomp_16s_fall_ef)))
colnames(asv_16s_fall_ef_abun) <- "MeanAbun"
asv_16s_spring_di_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_spring_di)) / ncol(t(asv.n0.acomp_16s_spring_di)))
colnames(asv_16s_spring_di_abun) <- "MeanAbun"
asv_16s_fall_di_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16s_fall_di)) / ncol(t(asv.n0.acomp_16s_fall_di)))
colnames(asv_16s_fall_di_abun) <- "MeanAbun"

# Bind columns for each subset
OccAbun_16s_spring_dwc <- bind_cols(asv_16s_spring_dwc_occ, asv_16s_spring_dwc_abun)
OccAbun_16s_fall_dwc <- bind_cols(asv_16s_fall_dwc_occ, asv_16s_fall_dwc_abun)
OccAbun_16s_spring_ef <- bind_cols(asv_16s_spring_ef_occ, asv_16s_spring_ef_abun)
OccAbun_16s_fall_ef <- bind_cols(asv_16s_fall_ef_occ, asv_16s_fall_ef_abun)
OccAbun_16s_spring_kr <- bind_cols(asv_16s_spring_kr_occ, asv_16s_spring_kr_abun)
OccAbun_16s_fall_kr <- bind_cols(asv_16s_fall_kr_occ, asv_16s_fall_kr_abun)
OccAbun_16s_spring_nft <- bind_cols(asv_16s_spring_nft_occ, asv_16s_spring_nft_abun)
OccAbun_16s_fall_nft <- bind_cols(asv_16s_fall_nft_occ, asv_16s_fall_nft_abun)
OccAbun_16s_spring_di <- bind_cols(asv_16s_spring_di_occ, asv_16s_spring_di_abun)
OccAbun_16s_fall_di <- bind_cols(asv_16s_fall_di_occ, asv_16s_fall_di_abun)

# Preserve sample-level relative abundance for all system-season combinations for alpha diversity
# Create a list of all your phyloseq subsets
physeq_list <- list(
  "spring_dwc" = physeq_16s_spring_dwc,
  "fall_dwc"   = physeq_16s_fall_dwc,
  "spring_ef"  = physeq_16s_spring_ef,
  "fall_ef"    = physeq_16s_fall_ef,
  "spring_kr"  = physeq_16s_spring_kr,
  "fall_kr"    = physeq_16s_fall_kr,
  "spring_nft" = physeq_16s_spring_nft,
  "fall_nft"   = physeq_16s_fall_nft,
  "spring_di"  = physeq_16s_spring_di,
  "fall_di"    = physeq_16s_fall_di
)

# Function to get sample-level relative abundance
get_sample_level_rel <- function(physeq_obj, system_name, season_name){
  asv_table <- t(as.data.frame(as(otu_table(physeq_obj), "matrix")))  # rows=samples
  # Replace zeros using Count Zero Multiplicative method
  asv_table_n0 <- t(cmultRepl(t(asv_table), label=0, method="CZM", output="p-counts", z.warning = 0.99))
  # Convert to compositions (relative abundance)
  asv_acomp <- as.data.frame(acomp(t(asv_table_n0)), total=1)
  # Add sample metadata columns
  asv_acomp$SampleID <- rownames(asv_acomp)
  asv_acomp$System   <- system_name
  asv_acomp$Season   <- season_name
  return(asv_acomp)
}

# Apply function to all phyloseq subsets
sample_level_rel_list <- mapply(
  get_sample_level_rel,
  physeq_list,
  system_name = c("DWC","DWC","EF","EF","KR","KR","NFT","NFT","DI","DI"),
  season_name = c("Spring","Fall","Spring","Fall","Spring","Fall","Spring","Fall","Spring","Fall"),
  SIMPLIFY = FALSE
)

# Combine all system-season sample-level relative abundance tables
all_samples_rel <- bind_rows(sample_level_rel_list)

#Reorder columns 
all_samples_rel <- all_samples_rel %>%
  dplyr::select(SampleID, System, Season, dplyr::everything())

# Reformat the data
all_samples_long <- all_samples_rel %>%
  pivot_longer(
    cols = starts_with("ASV"),
    names_to = "ASV",
    values_to = "RelAbundance"
  )
  
# Join taxonomy information
all_samples_long <- all_samples_long %>%
  left_join(taxon_16s_df, by = c("ASV" = "ASV"))

# Add a column with day
all_samples_long <- all_samples_long %>%
  mutate(
    day = as.numeric(sub(".*_(\\d+)$", "\\1", SampleID))
  )

# Select only ASV abundance columns (exclude metadata)
asv_cols <- setdiff(colnames(all_samples_rel), c("SampleID", "System", "Season"))
asv_table <- all_samples_rel[, asv_cols]

# convert to matrix 
asv_matrix <- as.matrix(sapply(asv_table, as.numeric))# Replace NA with 0
asv_matrix[is.na(asv_matrix)] <- 0

# Calculate alpha diversity metrics
alpha_div <- data.frame(
  SampleID = all_samples_rel$SampleID,
  System   = all_samples_rel$System,
  Season   = all_samples_rel$Season,
  Observed = specnumber(asv_matrix),      # observed ASVs
  Shannon  = diversity(asv_matrix, index = "shannon"),
  Simpson  = diversity(asv_matrix, index = "simpson")
)

alpha_div <- alpha_div %>%
  mutate(Day = as.numeric(sub(".*_(\\d+)$", "\\1", SampleID)))
  
#Plot occupancy-abundance
# Graph occupancy abundance data
# Prepare taxon data frame - add asv as a column name
taxon_16s <- as.data.frame(taxon_16s)
taxon_16s$asv <- rownames(taxon_16s)

  library(ggrepel)

  # Make the graph for dwc spring
  # Add asv as the column name
  OccAbun_16s_spring_dwc$asv <- rownames(OccAbun_16s_spring_dwc)
  # Join with taxon data
  OccAbun_taxa_spring_dwc <- left_join(OccAbun_16s_spring_dwc, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_spring_dwc_genus <- OccAbun_taxa_spring_dwc %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_spring_dwc <- ggplot(OccAbun_taxa_spring_dwc_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - DWC Spring") +
    theme(plot.title = element_text(hjust = 0.45))
    theme_minimal()
  Plot_OccAbun_spring_dwc
  ggsave(
    "Plot_OccAbun_spring_dwc.pdf",
    Plot_OccAbun_spring_dwc,
    width = 7,
    height = 5,
    dpi = 600
  )

  # Make the graph for dwc fall
  # Add asv as the column name
  OccAbun_16s_fall_dwc$asv <- rownames(OccAbun_16s_fall_dwc)
  # Join with taxon data
  OccAbun_taxa_fall_dwc <- left_join(OccAbun_16s_fall_dwc, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_fall_dwc_genus <- OccAbun_taxa_fall_dwc %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_fall_dwc <- ggplot(OccAbun_taxa_fall_dwc_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - DWC Fall") +
    theme(plot.title = element_text(hjust = 0.45))
    theme_minimal()
  Plot_OccAbun_fall_dwc
  ggsave(
    "Plot_OccAbun_fall_dwc.pdf",
    Plot_OccAbun_fall_dwc,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for KR spring
  # Add asv as the column name
  OccAbun_16s_spring_kr$asv <- rownames(OccAbun_16s_spring_kr)
  # Join with taxon data
  OccAbun_taxa_spring_kr <- left_join(OccAbun_16s_spring_kr, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_spring_kr_genus <- OccAbun_taxa_spring_kr %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_spring_kr <- ggplot(OccAbun_taxa_spring_kr_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - KR Spring") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_spring_kr
  ggsave(
    "Plot_OccAbun_spring_kr.pdf",
    Plot_OccAbun_spring_kr,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for kr fall
  # Add asv as the column name
  OccAbun_16s_fall_kr$asv <- rownames(OccAbun_16s_fall_kr)
  # Join with taxon data
  OccAbun_taxa_fall_kr <- left_join(OccAbun_16s_fall_kr, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_fall_kr_genus <- OccAbun_taxa_fall_kr %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_fall_kr <- ggplot(OccAbun_taxa_fall_kr_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - KR Fall") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_fall_kr
  ggsave(
    "Plot_OccAbun_fall_kr.pdf",
    Plot_OccAbun_fall_kr,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for nft spring
  # Add asv as the column name
  OccAbun_16s_spring_nft$asv <- rownames(OccAbun_16s_spring_nft)
  # Join with taxon data
  OccAbun_taxa_spring_nft <- left_join(OccAbun_16s_spring_nft, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_spring_nft_genus <- OccAbun_taxa_spring_nft %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_spring_nft <- ggplot(OccAbun_taxa_spring_nft_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf, segment.color = "grey60", segment.alpha = 0.4,) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundace of Bacterial Taxa - NFT Spring") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_spring_nft
  ggsave(
    "Plot_OccAbun_spring_nft.pdf",
    Plot_OccAbun_spring_nft,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for nft fall
  # Add asv as the column name
  OccAbun_16s_fall_nft$asv <- rownames(OccAbun_16s_fall_nft)
  # Join with taxon data
  OccAbun_taxa_fall_nft <- left_join(OccAbun_16s_fall_nft, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_fall_nft_genus <- OccAbun_taxa_fall_nft %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_fall_nft <- ggplot(OccAbun_taxa_fall_nft_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - NFT Fall") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_fall_nft
  ggsave(
    "Plot_OccAbun_fall_nft.pdf",
    Plot_OccAbun_fall_nft,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for ef spring
  # Add asv as the column name
  OccAbun_16s_spring_ef$asv <- rownames(OccAbun_16s_spring_ef)
  # Join with taxon data
  OccAbun_taxa_spring_ef <- left_join(OccAbun_16s_spring_ef, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_spring_ef_genus <- OccAbun_taxa_spring_ef %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_spring_ef <- ggplot(OccAbun_taxa_spring_ef_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - EF Spring") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_spring_ef
  ggsave(
    "Plot_OccAbun_spring_ef.pdf",
    Plot_OccAbun_spring_ef,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for ef fall
  # Add asv as the column name
  OccAbun_16s_fall_ef$asv <- rownames(OccAbun_16s_fall_ef)
  # Join with taxon data
  OccAbun_taxa_fall_ef <- left_join(OccAbun_16s_fall_ef, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_fall_ef_genus <- OccAbun_taxa_fall_ef %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_fall_ef <- ggplot(OccAbun_taxa_fall_ef_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - EF Fall") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_fall_ef
  ggsave(
    "Plot_OccAbun_fall_ef.pdf",
    Plot_OccAbun_fall_ef,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for di spring
  # Add asv as the column name
  OccAbun_16s_spring_di$asv <- rownames(OccAbun_16s_spring_di)
  # Join with taxon data
  OccAbun_taxa_spring_di <- left_join(OccAbun_16s_spring_di, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_spring_di_genus <- OccAbun_taxa_spring_di %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_spring_di <- ggplot(OccAbun_taxa_spring_di_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf, box.padding = 1, segment.color = "grey60", segment.alpha = 0.4, force = 2, force_pull = 0.5) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - DI Spring") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_spring_di
  ggsave(
    "Plot_OccAbun_spring_di.pdf",
    Plot_OccAbun_spring_di,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Make the graph for di fall
  # Add asv as the column name
  OccAbun_16s_fall_di$asv <- rownames(OccAbun_16s_fall_di)
  # Join with taxon data
  OccAbun_taxa_fall_di <- left_join(OccAbun_16s_fall_di, taxon_16s, by = "asv")
  # Aggregate by genus and graph
  OccAbun_taxa_fall_di_genus <- OccAbun_taxa_fall_di %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarize(
      MeanAbun = mean(MeanAbun),
      MeanOcc = mean(MeanOcc),
      .groups = "drop"
    ) %>%
    dplyr::filter(MeanAbun > 0.005, MeanOcc > 0.1)
  # Create the plot
  Plot_OccAbun_fall_di <- ggplot(OccAbun_taxa_fall_di_genus, aes(x = MeanAbun, y = MeanOcc)) +
    geom_point() +
    geom_text_repel(aes(label = Genus), size = 3, max.overlaps = Inf, box.padding = 1, segment.color = "grey60", segment.alpha = 0.4, force = 2, force_pull = 0.5) +
    scale_x_log10() +
    xlab("log10 Mean Relative Abundance") +
    ylab("Occupancy") +
    ggtitle("Occupancy Abundance of Bacterial Taxa - DI Fall") +
    theme(plot.title = element_text(hjust = 0.45)) 
    theme_minimal()
  Plot_OccAbun_fall_di
  ggsave(
    "Plot_OccAbun_fall_di.pdf",
    Plot_OccAbun_fall_di,
    width = 7,
    height = 5,
    dpi = 600
  )
  
  # Plot abundance with a heat map
  # Only including asvs with an occupancy over 0.5 and abundance over 0.005 as filtered during occupancy-abundance plotting
  # List of your occupancy-abundance dataframes
  filtered_genus_dfs <- list(
    "DI Fall"     = OccAbun_taxa_fall_di_genus,
    "DWC Fall"    = OccAbun_taxa_fall_dwc_genus,
    "EF Fall"     = OccAbun_taxa_fall_ef_genus,
    "KR Fall"     = OccAbun_taxa_fall_kr_genus,
    "NFT Fall"    = OccAbun_taxa_fall_nft_genus,
    "DI Spring"   = OccAbun_taxa_spring_di_genus,
    "DWC Spring"  = OccAbun_taxa_spring_dwc_genus,
    "EF Spring"   = OccAbun_taxa_spring_ef_genus,
    "KR Spring"   = OccAbun_taxa_spring_kr_genus,
    "NFT Spring"  = OccAbun_taxa_spring_nft_genus
  )
  
  # Merge all datasets
  merged_genus <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Genus"),
                         lapply(names(filtered_genus_dfs), function(name) {
                           df <- filtered_genus_dfs[[name]] %>%
                             dplyr::select(Genus, MeanAbun, MeanOcc) %>%
                             dplyr::rename(
                               !!paste0(name, "_Abun") := MeanAbun,
                               !!paste0(name, "_Occ")  := MeanOcc
                             )
                           df
                         })
  )
  
  # Set genus names as rownames, then remove the Genus column
  merged_genus <- as.data.frame(merged_genus)
  rownames(merged_genus) <- merged_genus$Genus
  merged_genus <- merged_genus[, -which(colnames(merged_genus) == "Genus")]
  
  merged_genus[is.na(merged_genus)] <- 0
  
  install.packages("pheatmap")
  library(pheatmap)
  
  # Specify desired column order
  col_order <- c(
    "DI Fall", "DI Spring",
    "DWC Fall", "DWC Spring",
    "EF Fall", "EF Spring",
    "KR Fall", "KR Spring",
    "NFT Fall", "NFT Spring"
  )
  
 # The existing occupancy and abundance cut-off is too low
  # STEP 1: Make a subset of each dataframe for taxa with >1% mean abundance and >0.7 mean occupancy 
    filtered_genus_dfs_1 <- lapply(filtered_genus_dfs, function(df) {
    df %>%
      dplyr::filter(
        MeanAbun > 0.01,
        MeanOcc  >= 0.8
      )
  })
  
  # STEP 2: Merge all datasets for the >1% subset
  merged_genus_1 <- Reduce(function(x, y) full_join(x, y, by = "Genus"),
                           lapply(names(filtered_genus_dfs_1), function(name) {
                             df <- filtered_genus_dfs_1[[name]]
                             df <- df[, c("Genus", "MeanAbun")]
                             colnames(df)[2] <- name
                             return(df)
                           }))
  
  # STEP 3: Set genus names as rownames, then remove the Genus column
  merged_genus_1 <- as.data.frame(merged_genus_1)
  rownames(merged_genus_1) <- merged_genus_1$Genus
  merged_genus_1 <- merged_genus_1[, -which(colnames(merged_genus_1) == "Genus")]
  
  # STEP 4: Replace NAs with 0 for plotting
  merged_genus_1[is.na(merged_genus_1)] <- 0
  
  # Step 5: Make a mapping to ensure column names match col_order
  # Step 5: Fix column names to match col_order
  colnames(merged_genus_1) <- sub("^Di", "DI", colnames(merged_genus_1))
  colnames(merged_genus_1) <- sub("^Dwc", "DWC", colnames(merged_genus_1))
  colnames(merged_genus_1) <- sub("^Ef", "EF", colnames(merged_genus_1))
  colnames(merged_genus_1) <- sub("^Kr", "KR", colnames(merged_genus_1))
  colnames(merged_genus_1) <- sub("^Nft", "NFT", colnames(merged_genus_1))
  
  # STEP 5b: desired column order
  col_order <- c(
    "DI Fall",  "DI Spring",
    "DWC Fall", "DWC Spring",
    "EF Fall",  "EF Spring",
    "KR Fall",  "KR Spring",
    "NFT Fall", "NFT Spring"
  )
  
  # STEP 6: Reorder columns
  missing_cols <- setdiff(col_order, colnames(merged_genus_1))
  if (length(missing_cols) > 0) {
    stop("These columns from col_order are missing in merged_genus_1: ",
         paste(missing_cols, collapse = ", "))
  }
  merged_genus_1 <- merged_genus_1[, col_order]
  # STEP 8: Create the heatmap and save with Cairo
  library(ComplexHeatmap)
  library(circlize)
  Cairo::CairoPNG(
    filename = "genus_heatmap_1pct.png",
    width = 8, height = 10, units = "in",
    dpi = 600
  )
  
  ht <- Heatmap(
    as.matrix(merged_genus_1),
    name = "Mean Relative Abundance",
    row_names_gp = gpar(fontface = "italic", fontsize = 10.5),
    column_names_gp = gpar(fontsize = 10),
    col = colorRamp2(
      c(0,
        max(merged_genus_1) * 0.5,
        max(merged_genus_1)),
      c("#440154", "#21908C", "#FDE725")
    ),
    column_names_side = "bottom",
    column_names_rot = 90,
    column_names_centered = TRUE,
    width = unit(ncol(merged_genus_1) * 8, "mm"),
    height = unit(nrow(merged_genus_1) * 4, "mm"),
    rect_gp = gpar(col = "black"),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter"
    )
  )
  
  draw(ht)
  dev.off()
  
  
  # Get the most abundant asv for each system-season
  # Get the names of the most abundant ASV per column
  top_asv_names <- apply(merged_genus, 2, function(col) {
    rownames(merged_genus)[which.max(col)]
  })
  
  # Get the max value (i.e., highest relative abundance) per column
  top_asv_values <- apply(merged_genus_1, 2, max)
  
  # Convert to percent
  top_asv_percents <- apply(merged_genus_1, 2, function(col) {
    round(max(col) * 100, 2)
  })
  
  # Combine into a data frame
  top_asvs_df <- data.frame(
    System = colnames(merged_genus_1),
    Top_ASV = top_asv_names,
    Percent = top_asv_percents
  )
  
  print(top_asvs_df)
  
# Subset asv with Occupancy=1
  asv_16sdwc_s_occ1 <- subset(OccAbun_16s_spring_dwc, MeanOcc == 1) %>% arrange(rownames(.))
  asv_16sdwc_f_occ1 <- subset(OccAbun_16s_fall_dwc,   MeanOcc == 1) %>% arrange(rownames(.))
  
  asv_16sef_s_occ1  <- subset(OccAbun_16s_spring_ef,  MeanOcc == 1) %>% arrange(rownames(.))
  asv_16sef_f_occ1  <- subset(OccAbun_16s_fall_ef,    MeanOcc == 1) %>% arrange(rownames(.))
  
  asv_16skr_s_occ1  <- subset(OccAbun_16s_spring_kr,  MeanOcc == 1) %>% arrange(rownames(.))
  asv_16skr_f_occ1  <- subset(OccAbun_16s_fall_kr,    MeanOcc == 1) %>% arrange(rownames(.))
  
  asv_16snft_s_occ1 <- subset(OccAbun_16s_spring_nft, MeanOcc == 1) %>% arrange(rownames(.))
  asv_16snft_f_occ1 <- subset(OccAbun_16s_fall_nft,   MeanOcc == 1) %>% arrange(rownames(.))
  
  asv_16sdi_s_occ1  <- subset(OccAbun_16s_spring_di,  MeanOcc == 1) %>% arrange(rownames(.))
  asv_16sdi_f_occ1  <- subset(OccAbun_16s_fall_di,    MeanOcc == 1) %>% arrange(rownames(.))

# Add taxonomic information
  asv_16sdwc_s_occ1_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdwc_s_occ1))
  asv_16sdwc_f_occ1_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdwc_f_occ1))
  
  asv_16sef_s_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sef_s_occ1))
  asv_16sef_f_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sef_f_occ1))
  
  asv_16skr_s_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16skr_s_occ1))
  asv_16skr_f_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16skr_f_occ1))
  
  asv_16snft_s_occ1_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16snft_s_occ1))
  asv_16snft_f_occ1_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16snft_f_occ1))
  
  asv_16sdi_s_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdi_s_occ1))
  asv_16sdi_f_occ1_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdi_f_occ1))
 
# Bind columns
  asv_16sdwc_s_occ1 <- bind_cols(asv_16sdwc_s_occ1, asv_16sdwc_s_occ1_tax)
  asv_16sdwc_f_occ1 <- bind_cols(asv_16sdwc_f_occ1, asv_16sdwc_f_occ1_tax)
  
  asv_16sef_s_occ1  <- bind_cols(asv_16sef_s_occ1,  asv_16sef_s_occ1_tax)
  asv_16sef_f_occ1  <- bind_cols(asv_16sef_f_occ1,  asv_16sef_f_occ1_tax)
  
  asv_16skr_s_occ1  <- bind_cols(asv_16skr_s_occ1,  asv_16skr_s_occ1_tax)
  asv_16skr_f_occ1  <- bind_cols(asv_16skr_f_occ1,  asv_16skr_f_occ1_tax)
  
  asv_16snft_s_occ1 <- bind_cols(asv_16snft_s_occ1, asv_16snft_s_occ1_tax)
  asv_16snft_f_occ1 <- bind_cols(asv_16snft_f_occ1, asv_16snft_f_occ1_tax)
  
  asv_16sdi_s_occ1  <- bind_cols(asv_16sdi_s_occ1,  asv_16sdi_s_occ1_tax)
  asv_16sdi_f_occ1  <- bind_cols(asv_16sdi_f_occ1,  asv_16sdi_f_occ1_tax)
  
# Combine all core taxa into one dataframe
  # Extract core ASV IDs (rownames) from each core dataset
  core_16sdwc_s <- rownames(asv_16sdwc_s_occ1)
  core_16sdwc_f <- rownames(asv_16sdwc_f_occ1)
  
  core_16sef_s <- rownames(asv_16sef_s_occ1)
  core_16sef_f <- rownames(asv_16sef_f_occ1)
  
  core_16skr_s <- rownames(asv_16skr_s_occ1)
  core_16skr_f <- rownames(asv_16skr_f_occ1)
  
  core_16snft_s <- rownames(asv_16snft_s_occ1)
  core_16snft_f <- rownames(asv_16snft_f_occ1)
  
  core_16sdi_s <- rownames(asv_16sdi_s_occ1)
  core_16sdi_f <- rownames(asv_16sdi_f_occ1)
  # Combine
  core_16s_all <- data.frame(
    OTU = c(core_16sdwc_s, core_16sdwc_f,
            core_16sef_s, core_16sef_f,
            core_16skr_s, core_16skr_f,
            core_16snft_s, core_16snft_f,
            core_16sdi_s, core_16sdi_f),
    System = c(
      rep("DWC", length(core_16sdwc_s)),
      rep("DWC", length(core_16sdwc_f)),
      
      rep("EF", length(core_16sef_s)),
      rep("EF", length(core_16sef_f)),
      
      rep("KR", length(core_16skr_s)),
      rep("KR", length(core_16skr_f)),
      
      rep("NFT", length(core_16snft_s)),
      rep("NFT", length(core_16snft_f)),
      
      rep("DI", length(core_16sdi_s)),
      rep("DI", length(core_16sdi_f))
    ),
    Season = c(
      rep("Spring", length(core_16sdwc_s)),
      rep("Fall", length(core_16sdwc_f)),
      
      rep("Spring", length(core_16sef_s)),
      rep("Fall", length(core_16sef_f)),
      
      rep("Spring", length(core_16skr_s)),
      rep("Fall", length(core_16skr_f)),
      
      rep("Spring", length(core_16snft_s)),
      rep("Fall", length(core_16snft_f)),
      
      rep("Spring", length(core_16sdi_s)),
      rep("Fall", length(core_16sdi_f))
    )
  )
  
  # Add taxonomy information
  core_16s_all_taxa <- subset(taxon_16s, rownames(taxon_16s) %in% core_16s_all$OTU)
  core_16s_all <- core_16s_all %>%
    left_join(
      core_16s_all_taxa %>%
        rownames_to_column(var = "OTU"),
      by = "OTU"
    )
  
  #Write csv file
  write.csv(core_16s_all_taxa, file="Core 16s.csv")  
  
# Visualize core taxa by Season
library(ggvenn)
  coretaxa_16s_system <- core_16s_all %>%
    group_by(System) %>%
    summarise(OTUs = list(unique(OTU))) %>%   # list of unique OTUs per system
    deframe()  # converts two-column tibble into a named list System -> OTU vector
  
 # A max of four groups can be compared, there are five systems
  ggvenn(coretaxa_16s_system[c("DWC", "KR", "NFT", "EF")],
         show_percentage = FALSE,
         fill_color = c('#2A9D8F', '#420A68','#FCA50A', '#E56F00'),
         set_name_size = 4,
         stroke_size = 0.5) +
    ggtitle("EF and DWC") +
    theme(plot.title = element_text(hjust = 0.5))

  # Visualize core taxa by season
    coretaxa_16s_season <- core_16s_all %>% 
      group_by(Season) %>% 
      summarise(OTUs = list(unique(OTU))) %>% 
      deframe()
  
        library(Cairo)
    Cairo::CairoPNG("core_venn.png", width = 6, height = 5, units = "in", dpi = 600)
    print(
      ggvenn(coretaxa_16s_season,
             show_percentage = FALSE,
             fill_color = c('#2A9D8F', "#420A68"),
             set_name_size = 6,
             stroke_size = 0.5) +
        theme(
          text = element_text(size = 80, face = "bold")  # bigger and bold
        )
    )
    dev.off()
   
    # Identify and modify only the number labels (not the set names)
    for (i in seq_along(p$layers)) {
      layer <- p$layers[[i]]
      if ("GeomText" %in% class(layer$geom) &&
          "label" %in% names(layer$mapping)) {
        
        # Apply size/bold only to the intersection number labels
        p$layers[[i]]$aes_params$size <- 15    
      }
    }
    
    # Step 3: Display the plot with no clipping so text outside the panel shows
    p +
      coord_cartesian(clip = "off") +  # allow labels outside the plot area
      theme(
        plot.margin = margin(t = 30, r = 10, b = 10, l = 10)  # increase top margin
      )
    
  # Identify if any genus is in all systems 
    # Create binary presence table: one row per OTU, one column per System
    otu_presence_matrix <- core_16s_all %>%
      distinct(System, OTU) %>%              # Keep only System-OTU pairs
      mutate(present = 1) %>%
      pivot_wider(names_from = System, values_from = present, values_fill = 0)
    
    # Get OTUs found in all systems (i.e., all columns = 1)
    otu_in_all_systems <- otu_presence_matrix %>%
      filter(if_all(-OTU, ~ . == 1))
    # View results - the only asv in all systems is ASV14 or Acidovorax
    nrow(otu_in_all_systems)
    otu_in_all_systems$OTU 
    
  # Visualize system overlaps with an UpSet
    library(UpSetR)
    # Convert to data.frame and remove OTU column
    otu_mat <- otu_presence_matrix %>%
      dplyr::select(-OTU) %>%
      as.data.frame()
    
    install.packages("Cairo")
    library(Cairo)
    # Save and create the graph
    Cairo::CairoPNG(
      filename = "upset.png",
      width = 18, height = 9, units = "in",
      dpi = 600
    )
    
    UpSetR::upset(
      otu_mat,
      sets = colnames(otu_mat),
      order.by = "freq",
      mainbar.y.label = "Shared Temporal Core ASVs",
      sets.x.label = "Total Temporal Core ASVs",
      main.bar.color = "#420A68",
      sets.bar.color = "#E56F00",
      text.scale = 2.4,
      keep.order = TRUE,
      mb.ratio = c(0.7, 0.3)  # shrink left-side bars
    )
    
    dev.off()
    
## Network analyses ##
  # Subset asv with Occupancy > 0.5
  asv_16sdwc_s_occ0.5 <- subset(OccAbun_16s_spring_dwc, MeanOcc >= 0.5) %>% arrange(rownames(.))
  asv_16sdwc_f_occ0.5 <- subset(OccAbun_16s_fall_dwc,   MeanOcc >= 0.5) %>% arrange(rownames(.))
  
  asv_16sef_s_occ0.5  <- subset(OccAbun_16s_spring_ef,  MeanOcc >= 0.5) %>% arrange(rownames(.))
  asv_16sef_f_occ0.5  <- subset(OccAbun_16s_fall_ef,    MeanOcc >= 0.5) %>% arrange(rownames(.))
  
  asv_16skr_s_occ0.5  <- subset(OccAbun_16s_spring_kr,  MeanOcc >= 0.5) %>% arrange(rownames(.))
  asv_16skr_f_occ0.5  <- subset(OccAbun_16s_fall_kr,    MeanOcc >= 0.5) %>% arrange(rownames(.))
  
  asv_16snft_s_occ0.5 <- subset(OccAbun_16s_spring_nft, MeanOcc >= 0.5) %>% arrange(rownames(.))
  asv_16snft_f_occ0.5 <- subset(OccAbun_16s_fall_nft,   MeanOcc >= 0.5) %>% arrange(rownames(.))
  
  asv_16sdi_s_occ0.5  <- subset(OccAbun_16s_spring_di,  MeanOcc >= 0.5) %>% arrange(rownames(.))
  asv_16sdi_f_occ0.5  <- subset(OccAbun_16s_fall_di,    MeanOcc >= 0.5) %>% arrange(rownames(.))
  
  # Add taxonomic information
  asv_16sdwc_s_occ0.5_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdwc_s_occ0.5))
  asv_16sdwc_f_occ0.5_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdwc_f_occ0.5))
  
  asv_16sef_s_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sef_s_occ0.5))
  asv_16sef_f_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sef_f_occ0.5))
  
  asv_16skr_s_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16skr_s_occ0.5))
  asv_16skr_f_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16skr_f_occ0.5))
  
  asv_16snft_s_occ0.5_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16snft_s_occ0.5))
  asv_16snft_f_occ0.5_tax  <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16snft_f_occ0.5))
  
  asv_16sdi_s_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdi_s_occ0.5))
  asv_16sdi_f_occ0.5_tax   <- subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sdi_f_occ0.5))
  
  # Combine all core taxa into one dataframe
  # Extract core ASV IDs (rownames) from each core dataset
  core_16sdwc_s_0.5 <- rownames(asv_16sdwc_s_occ0.5)
  core_16sdwc_f_0.5 <- rownames(asv_16sdwc_f_occ0.5)
  
  core_16sef_s_0.5 <- rownames(asv_16sef_s_occ0.5)
  core_16sef_f_0.5 <- rownames(asv_16sef_f_occ0.5)
  
  core_16skr_s_0.5 <- rownames(asv_16skr_s_occ0.5)
  core_16skr_f_0.5 <- rownames(asv_16skr_f_occ0.5)
  
  core_16snft_s_0.5 <- rownames(asv_16snft_s_occ0.5)
  core_16snft_f_0.5 <- rownames(asv_16snft_f_occ0.5)
  
  core_16sdi_s_0.5 <- rownames(asv_16sdi_s_occ0.5)
  core_16sdi_f_0.5 <- rownames(asv_16sdi_f_occ0.5)
 
   # Combine
  core_16s_all_0.5 <- data.frame(
    OTU = c(core_16sdwc_s_0.5, core_16sdwc_f_0.5,
            core_16sef_s_0.5, core_16sef_f_0.5,
            core_16skr_s_0.5, core_16skr_f_0.5,
            core_16snft_s_0.5, core_16snft_f_0.5,
            core_16sdi_s_0.5, core_16sdi_f_0.5),
    System = c(
      rep("DWC", length(core_16sdwc_s_0.5)),
      rep("DWC", length(core_16sdwc_f_0.5)),
      
      rep("EF", length(core_16sef_s_0.5)),
      rep("EF", length(core_16sef_f_0.5)),
      
      rep("KR", length(core_16skr_s_0.5)),
      rep("KR", length(core_16skr_f_0.5)),
      
      rep("NFT", length(core_16snft_s_0.5)),
      rep("NFT", length(core_16snft_f_0.5)),
      
      rep("DI", length(core_16sdi_s_0.5)),
      rep("DI", length(core_16sdi_f_0.5))
    ),
    Season = c(
      rep("Spring", length(core_16sdwc_s_0.5)),
      rep("Fall", length(core_16sdwc_f_0.5)),
      
      rep("Spring", length(core_16sef_s_0.5)),
      rep("Fall", length(core_16sef_f_0.5)),
      
      rep("Spring", length(core_16skr_s_0.5)),
      rep("Fall", length(core_16skr_f_0.5)),
      
      rep("Spring", length(core_16snft_s_0.5)),
      rep("Fall", length(core_16snft_f_0.5)),
      
      rep("Spring", length(core_16sdi_s_0.5)),
      rep("Fall", length(core_16sdi_f_0.5))
    )
  )
  
  # Add taxonomy information
  core_16s_all_taxa_0.5 <- subset(taxon_16s, rownames(taxon_16s) %in% core_16s_all_0.5$OTU)
  core_16s_all_0.5 <- core_16s_all_0.5 %>%
    left_join(
      core_16s_all_taxa_0.5 %>%
        rownames_to_column(var = "OTU"),
      by = "OTU"
    )
  
  #Write csv file
  write.csv(core_16s_all_taxa_0.5, file="Core 16s_0.5.csv")  
  
  # Identify if any taxa is in all systems
  # Create binary presence table: one row per OTU, one column per System
  otu_presence_matrix_0.5 <- core_16s_all_0.5 %>%
    distinct(System, OTU) %>%              # Keep only System-OTU pairs
    mutate(present = 1) %>%
    pivot_wider(names_from = System, values_from = present, values_fill = 0)
  # Get OTUs found in all systems (i.e., all columns = 1)
  otu_in_all_systems_0.5 <- otu_presence_matrix_0.5 %>%
    filter(if_all(-OTU, ~ . >= 0.5))
  # View results - the only asv in all systems is ASV14 or Acidovorax
  nrow(otu_in_all_systems_0.5)
  otu_in_all_systems_0.5$OTU
  
  # Combine list of core ASVs greater than 0.5
  nettaxa_16s_list <- unique(c(
    rownames(asv_16sdwc_s_occ0.5),
    rownames(asv_16sdwc_f_occ0.5),
    rownames(asv_16sef_s_occ0.5),
    rownames(asv_16sef_f_occ0.5),
    rownames(asv_16skr_s_occ0.5),
    rownames(asv_16skr_f_occ0.5),
    rownames(asv_16snft_s_occ0.5),
    rownames(asv_16snft_f_occ0.5),
    rownames(asv_16sdi_s_occ0.5),
    rownames(asv_16sdi_f_occ0.5)
  ))
  
  asv_16s.T<-t(asv_16s)
  
  # Subset the full asv table with unique core ASVs
  asvs_16s_net <- subset(asv_16s.T, rownames(asv_16s.T) %in% nettaxa_16s_list)
  
  #Deal with sparcity. Only include taxa that have an occupancy above 0.5 in the three facilities from previous analyses.
  
  #Make phyloseq for Network analyses
  # Subset taxon_16s to match ASVs with over 0.5 Occupancy 
  taxon_16s_net <- taxon_16s[rownames(taxon_16s) %in% rownames(asvs_16s_net), ]
  # Reorder taxonomy to have the same order as your ASV table rows
  taxon_16s_net <- taxon_16s_net[match(rownames(asvs_16s_net), rownames(taxon_16s_net)), ]
  # Convert taxonomy dataframe to matrix
  taxon_16s_net <- as.matrix(taxon_16s_net)
  # Make sure metadata has rownames as sample IDs
  rownames(metadata_16s) <- metadata_16s$SampleID  # replace 'SampleID' with your actual column
  # Subset metadata to match OTU table columns
  metadata_16s_net <- metadata_16s[colnames(asvs_16s_net), , drop = FALSE]
  # Reorder metadata to match OTU table columns
  metadata_16s_net <- metadata_16s[match(colnames(asvs_16s_net), rownames(metadata_16s)), , drop = FALSE]
  
  #Build phyloseq object
  phyloseq_16s_net<-phyloseq(otu_table(asvs_16s_net, taxa_are_rows = TRUE), tax_table(taxon_16s_net), sample_data(metadata_16s_net))  

  #### Networks  ####
  #Install NetCoMi package
  devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE,
                           repos = c("https://cloud.r-project.org/",
                                     BiocManager::repositories()))
  
  library(NetCoMi)
  
  #Install package dependencies
  installNetCoMiPacks()  
  
  #Construct network with spieceasi
  net_16s<-netConstruct(phyloseq_16s_net, 
                        measure='spieceasi',
                        measurePar = list(nlambda=10,lambda.min.ratio=1e-2,pulsar.params=list(rep.num=99)),
                        normMethod = "none",
                        zeroMethod = "none",
                        sparsMethod = "none",
                        verbose = 3,
                        seed=10000)
  
  #Analyze network. Detect hubs by highest betweenness centrality
  net_16s_a<-netAnalyze(net_16s,
                        centrLCC = TRUE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = "betweenness",
                        weightDeg = FALSE, normDeg = FALSE)
  
  #Get network summary
  network_summary_16s<-summary(net_16s_a, numbNodes=5L)
  
  #Plot network - color by cluster, highlight hubs
  cols16s<-c("#F37B59","#D89000","#AFA100","#00BF7D","#00BF4F","#BF80FF","#FF689F", "#009ACD", "#8FBC8F", "#CD5C5C")
  
  # Save network plot as PNG with Cairo
  CairoPNG("net_16s_plot.png",
           width = 8, height = 6, units = "in", dpi = 600)
  
  plot(net_16s_a, 
       nodeColor = "cluster", 
       nodeSize = "clr", 
       rmSingles = TRUE, 
       #labels = NULL,
       labelScale = TRUE,
       nodeSizeSpread = 1, 
       highlightHubs = TRUE,
       showTitle = FALSE, 
       cexLabels = 0,
       repulsion = 1.5,
       nodeTransp = 10, 
       hubTransp = 0,
       edgeTranspLow = 90,
       edgeTranspHigh = 50,
       colorVec = cols16s)
  
  dev.off()

  #Get hub ASV names
  hubs_16s<-net_16s_a$hubs$hubs1
  #Make table with hub taxonomy
  hub_16s_taxa<-subset(taxon_16s, rownames(taxon_16s) %in% hubs_16s)
  
# I am looking at the genera with the highest relative abundances by individual systems over time
 # Identify target genera
  target_genera <- c("Methylophilus", "Eoetvoesia", "Acidovorax")
 # Subset data
  focus_genera <- all_samples_long %>%
    filter(Genus %in% target_genera)
  
  # Order systems and convert day to numeric
  focus_genera <- focus_genera %>%
    mutate(
      System = factor(System, levels = c("DI", "DWC", "EF", "KR", "NFT", "NS")),
      day = as.numeric(day)
    )
  
  # Create combined label for ordered grouping
  focus_genera <- focus_genera %>%
    mutate(System_Day = paste0(System, "_", day))
  
  # Build ordered factor
  focus_genera$System_Day <- factor(
    focus_genera$System_Day,
    levels = focus_genera %>%
      arrange(System, day) %>%
      pull(System_Day) %>%
      unique()
  )
  
  # Plot using facet_wrap so each Season gets its own x-axis
  ggplot(focus_genera, aes(x = System_Day, y = RelAbundance, fill = System)) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_wrap(Season + Genus ~ ., scales = "free_x", ncol = length(unique(focus_genera$Genus))) +
    scale_x_discrete(labels = function(x) sub("^[^_]*_([0-9]+)$", "\\1", x)) +
    scale_fill_manual(values = c(
      "DI"  = "#41afaa",
      "DWC" = "#466eb4",
      "EF"  = "#00a0e1",
      "KR"  = "#e6a532",
      "NFT" = "#d7642c",
      "NS"  = "#af4b91"
    )) +
    labs(
      x = "Day",
      y = "Relative Abundance"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      strip.text = element_text(size = 12),
      panel.grid.major.x = element_blank()
    )
  
  # Save ggplot
  ggsave(
    filename = "focus_genera_plot.png",   
    plot = last_plot(),                    
    width = 18,                            
    height = 10,                            
    dpi = 600                              
  )
  
# Calculate standard deviation for Table 1  
  focus_genera %>%
    filter(System == "KR", Season == "Fall", Genus == "Acidovorax") %>%
    summarise(sd_rel_abundance = sd(RelAbundance, na.rm = TRUE))
  