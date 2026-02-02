#DADA2 pipeline applied on NESARE data
#From https://benjjneb.github.io/dada2/tutorial.html

#Last updated: AB 05/4/25
#This is looking at the data for the 2023-2024 bok choy soilless system experiment. 

#Install DADA2

#Attach libraries
library(dada2)
library(ggplot2)

#Set working directory
setwd('/storage/home/akb6471/work/hydroponics_rerun/Reads_forSRA/')

#Make path for fastq files
path_all <- '/storage/home/akb6471/work/hydroponics_rerun/Reads_forSRA/' 
list.files(path_all)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_all <- sort(list.files(path_all, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs_all <- sort(list.files(path_all, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs_all), "_"), function(x) {paste(x[1:5], collapse="_")})

#Inspect read quality for forward and reverse reads 
qc_F<-plotQualityProfile(fnFs_all[28])
ggsave("qc_F.png", plot=qc_F, device='png')
qc_R<-plotQualityProfile(fnRs_all[28])
ggsave("qc_R.png", plot=qc_R, device='png')

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs_all <- file.path(path_all, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs_all <- file.path(path_all, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs_all) <- sample.names
names(filtRs_all) <- sample.names

#Trim and filter for Forward and Reverse at 240, 240 - Check based on graph
out_all <- filterAndTrim(fnFs_all, filtFs_all, fnRs_all, filtRs_all, truncLen=c(240,240),
                         maxN=0, maxEE=c(2,2), truncQ=3, rm.phix=TRUE,
                         compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out_all

#Data Statistics after Trimming 
sum(out_all[,1]) #total reads in
sum(out_all[,2]) #total reads out
sum(out_all[,1]) - sum(out_all[,2]) #reads lost
sum(out_all[,2])/sum(out_all[,1]) 

#Learn error rates
errF_all <- learnErrors(filtFs_all, multithread=TRUE)
errR_all <- learnErrors(filtRs_all, multithread=TRUE)

error_F<-plotErrors(errF_all, nominalQ=TRUE) 
error_R<-plotErrors(errR_all, nominalQ=TRUE)
ggsave("error_F.png", plot=error_F, device='png')
ggsave("error_R.png", plot=error_R, device='png')

#Sample inference (adjust the sequences based on learned error rates)
dadaFs_all <- dada(filtFs_all, err=errF_all, multithread=TRUE)
dadaRs_all <- dada(filtRs_all, err=errR_all, multithread=TRUE)

dadaFs_all[[1]]

#Merge paired reads
mergers_all <- mergePairs(dadaFs_all, filtFs_all, dadaRs_all, filtRs_all, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_all[[1]])

#construct sequence table
seqtab_all <- makeSequenceTable(mergers_all)
dim(seqtab_all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_all)))

#Remove non-target sequences from the sequence table
seqtab2_all <- seqtab_all[,nchar(colnames(seqtab_all)) %in% 251:256]

#Remove chimeras
seqtab.nochim_all <- removeBimeraDenovo(seqtab2_all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_all)
sum(seqtab.nochim_all)/sum(seqtab_all)

#Track reads through the pipeline - check that there is no great drop in reads throughout any of the previous steps
getN <- function(x) sum(getUniques(x))
track_all <- cbind(out_all, sapply(dadaFs_all, getN), sapply(dadaRs_all, getN), sapply(mergers_all, getN), rowSums(seqtab.nochim_all))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_all) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_all) <- sample.names
track_all

#Assign taxonomy
taxa_all <- assignTaxonomy(seqtab.nochim_all, "/storage/work/akb6471/hydroponics_rerun/Reads_forSRA/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

#Add species column
taxa_all <- addSpecies(taxa_all, "/storage/work/akb6471/hydroponics_rerun/Reads_forSRA/silva_v138.2_assignSpecies.fa.gz")
taxa.print_all <- taxa_all # Removing sequence rownames for display only
rownames(taxa.print_all) <- NULL
taxa.print_all

#Make phyloseq
library(phyloseq)
ps_all <- phyloseq(otu_table(seqtab.nochim_all, taxa_are_rows=FALSE), 
                   tax_table(taxa_all))

#Change name from sequence to ASV#
dna_all <- Biostrings::DNAStringSet(taxa_names(ps_all))
names(dna_all) <- taxa_names(ps_all)
ps_all <- merge_phyloseq(ps_all, dna_all)
taxa_names(ps_all) <- paste0("ASV", seq(ntaxa(ps_all)))
ps_all

#Save ASV table as csv
asvs_all<-as.data.frame(otu_table(ps_all))
Taxon_all<-as.data.frame(tax_table(ps_all))
DNA<-as.data.frame(refseq(ps_all))

write.csv(asvs_all, file = "ASV_all.csv")
write.csv(Taxon_all, file = "Taxon_all.csv")
write.csv(DNA, file="DNA_all.csv")

### Incorporating code from Jasna Kovac's script at 130 for analysis

# Install CRAN packages
install.packages(c("ape", "vegan", "cowplot", "tidyr", "dplyr", "compositions", "zCompositions", "viridis", "readxl", "psych", "svglite", "pairwiseAdonis"))

# Install additional packages
install.packages("remotes")
install.packages("tidyverse")
install.packages("zCompositions")
install.packages("vegan")

# Install pairwiseAdonis from GitHub
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#Load package
library(ggplot2)
library(microshades)
library(Polychrome)
library(colorspace)

#remove NAs from the sequence names in the ASV_all.csv file to match the metadata file names
#Import data - 16s rRNA V4
asvs_16s<-read.csv('ASV_all.csv', header = TRUE, row.names = 1)
taxon_16s<-read.csv('Taxon_all.csv', header = TRUE, row.names = 1)
metadata_16s<-read_excel('metadata.xlsx')
metadata_16s <- as.data.frame(metadata_16s)
rownames(metadata_16s) <- metadata_16s[[1]]
metadata_16s <- metadata_16s[, -1]

#Clean up taxonomy tables
#Add '_unclassified' marker to NAs in taxonomy table
taxon_16s$Phylum<-ifelse(is.na(taxon_16s$Phylum), paste(taxon_16s$Kingdom, "unclassified", sep = '_'), taxon_16s$Phylum)
taxon_16s$Class<-ifelse(is.na(taxon_16s$Class), paste(taxon_16s$Phylum, "unclassified", sep = '_'), taxon_16s$Class)
taxon_16s$Order<-ifelse(is.na(taxon_16s$Order), paste(taxon_16s$Class, "unclassified", sep = '_'), taxon_16s$Order)
taxon_16s$Family<-ifelse(is.na(taxon_16s$Family), paste(taxon_16s$Order, "unclassified", sep = '_'), taxon_16s$Family)
taxon_16s$Genus<-ifelse(is.na(taxon_16s$Genus), paste(taxon_16s$Family, "unclassified", sep = '_'), taxon_16s$Genus)
taxon_16s$Species<-ifelse(is.na(taxon_16s$Species), paste(taxon_16s$Genus, "unclassified", sep = '_'), taxon_16s$Species)

#Remove extra _unclassified
taxon_16s$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Class)
taxon_16s$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Species)

#Convert ASV and taxonomy tables to a matrices
asvs_16s<-as.matrix(asvs_16s)
taxon_16s<-as.matrix(taxon_16s)

#Make phyloseq object
ps_16s<-phyloseq(otu_table(asvs_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))

#Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- ps_16s %>%  subset_taxa( Order!="Chloroplast" | is.na(Order) )
physeq_16s <- physeq_16s %>% subset_taxa( Family!= "Mitochondria" | is.na("Family"))

#Get ASV table from phyloseq object
#this is where asv_16s is introduced
asv_16s<-as.data.frame(t(otu_table(physeq_16s)))
tail(rowSums(asv_16s))

#Get taxonomy table from phyloseq object
taxon.16s<-as.matrix(tax_table(physeq_16s))

#Remove ASVs with zero counts in all samples
asv_16s<-asv_16s[ which(rowSums(asv_16s)>0),]
asv_16s<-t(asv_16s)

#Save ASV, taxon, and metadata files to use in downstream analyses
write.csv(asv_16s, file="ASV_16s_clean.csv")
write.csv(taxon.16s, file="Taxon_16s_clean.csv")
write.csv(metadata_16s, file="metadata_16s_clean.csv")

#### COMPOSITIONAL ANALYSIS AT ASV LEVEL ####
#Based on Microbiome Analysis in R. Chap 10.

#Convert ASV table to appropriate format
#Samples need to be in rows and ASVs in columns 
asv_16s <- read.csv("ASV_16s_clean.csv", header=TRUE)
head(asv_16s) 

# Convert all columns to numeric, ensuring the dataframe structure is preserved
for (i in 2:ncol(asv_16s)) {
  asv_16s[,i] <- as.numeric(asv_16s[,i])
}

# make the first column a label instead
rownames(asv_16s) <- asv_16s[, 1]
asv_16s <- asv_16s[-1]

# CLR transform the data
#Step 1: Replace zero values before clr transformation
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
library(zCompositions)

# Testing with a cut-off of 99% or more 0s (removes 23 samples)
asv.n0_16s_99 <- t(cmultRepl(asv_16s, label=0, method="CZM", output="p-counts", z.warning = 0.99)) # No. of adjusted imputations: 47,907

#Step 2: Convert data to proportions 
asv.n0_16s_prop_99<-apply(asv.n0_16s_99, 2, function(x) {x/sum(x)})

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv.n0_16s_99<-ifelse(asv.n0_16s_99 < 0, asv.no_16s_99*(-1), asv.n0_16s_99)

#Step 3: Perform abundance and sample filtering
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16s_prop_f_99<-asv.n0_16s_99[apply(asv.n0_16s_prop_99, 1, min) > 0.000001, ]

#Check that samples are in columns and ASVs in rows
head(asv.n0_16s_prop_f_99)

#Step 4: perform CLR transformation
asv.n0.clr_16s_99<-t(apply(asv.n0_16s_prop_f_99, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16s_99)

#Remove nutrient solution samples from the data set
asv.n0.clr_16s_99_c <- asv.n0.clr_16s_99[1:(nrow(asv.n0.clr_16s_99) - 12), ]
#Remove bc samples from the data set
asv.n0.clr_16s_99_bc <- asv.n0.clr_16s_99_c[-c(5, 11, 24, 32, 40, 52, 53, 67), ]
# Remove samples with less than 1,000 sequencing reads
# DI, DW, EF, NF Day 35
# DW Day 0, 33
asv.n0.clr_16s_99_cl <- asv.n0.clr_16s_99_bc[-c(9, 10, 17, 20, 34, 59), ]

# Save the CLR transformed data
write.csv(asv.n0.clr_16s_99_cl, file="asv.n0.clr_16s.csv")
#FIX METADATA TABLE: remove samples that were excluded at the zero imputation step
metadata_16s_99 <- read.csv('metadata_16s_updated99.csv', header=TRUE, row.names=1)

#Filter the original asv_16s file to only include samples in asv.n0.clr_16s (raw counts are needed for ALDEx2 analysis)
# Step 1: Get the common ASVs (rows) and samples (columns) between the original ASV table and CLR-transformed table
common_asvs <- intersect(rownames(asvs_16s), rownames(asv.n0.clr_16s_99_cl))  # ASVs in both tables
common_samples <- intersect(colnames(asvs_16s), colnames(asv.n0.clr_16s_99_cl))  # Samples in both tables
# Step 2: Subset the original ASV table to retain only the common ASVs and samples
asv_16s_clean <- asvs_16s[common_asvs, common_samples]
# Save the filtered output
write.csv(asv_16s_clean, file="asv_16s_filtered.csv")

# Step 5: Perform Singular Value Decomposition (PCA)
pc.clr_16s<-prcomp(asv.n0.clr_16s_99_cl) #Doesn't give error in R 4.1, but gives error in R 4.0

#library(ade4) #This library contains the PCA function that doesn't give an error. Alternatively use prcomp()
#pc.clr_16s<-dudi.pca(asv.n0.clr_16s, scannf= FALSE, nf=5)

png("Screeplot - PCA2.png", width = 400, height = 300, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
screeplot(pc.clr_16s, type='barplot', main="Bacteria")
dev.off()

#Calculate total variance in the data (if this doesn't work, reload compositions library)
library(compositions)
mvar.clr_16s<-mvar(asv.n0.clr_16s_99_cl)

#Display results - 16s 
library(dplyr)
row_16s<-rownames(asv.n0.clr_16s_99_cl) #Make vector with sample names
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:2]) #Get PC1 and PC2 

#Combine and reformat PCA output
pc_out_meta_16s<-as.data.frame(bind_cols(pc_out_16s,metadata_16s_99)) #Add metadata information
row.names(pc_out_meta_16s)<-row_16s #Add rownames to dataframe
pc_out_meta_16s$system<-as.factor(pc_out_meta_16s$system)
pc_out_meta_16s$season<-as.factor(pc_out_meta_16s$season)

# Make PCA plot - First 2 axis 
library(ggplot2)

#Fig 1 - color by system

scale_color_manual(values = sys_palette)

# Recode system and season labels
pc_out_meta_16s$system <- recode(pc_out_meta_16s$system,
                                 'nf' = 'NFT',
                                 'di' = 'DI',
                                 'dwc' = 'DWC',
                                 'ef' = 'EF',
                                 'kr' = 'KR')

pc_out_meta_16s$season <- recode(pc_out_meta_16s$season,
                                 'f' = 'Fall',
                                 's' = 'Spring')

PCA_16s <- ggplot(pc_out_meta_16s, aes(x = PC1, y = PC2, color = system)) +
  geom_point(aes(shape = season), size = 3) +
  scale_color_manual(values = sys_palette, name = "Hydroponic System") +
  scale_x_continuous(
    name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2 / mvar.clr_16s * 100, digits = 1), "%", sep = "")
  ) +
  scale_y_continuous(
    name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2 / mvar.clr_16s * 100, digits = 1), "%", sep = "")
  ) +
  labs(
    color = "Hydroponic System",
    shape = "Season"
  ) +
  theme(
    legend.text = element_text(size = 15, color = 'black'),
    legend.title = element_text(size = 15, face = 'italic'),
    legend.position = 'right',
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 15, color = 'black'),
    panel.background = element_rect(fill = 'white', color = NA),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 1)
  )

print(PCA_16s)
ggsave("Fig4.png", plot = PCA_16s, device = "png", width = 8, height = 5, units = "in", dpi = 600)

library(vegan)

### PERMANOVA analysis ###

install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

dist_16s_99 <- dist(asv.n0.clr_16s_99_cl, method = 'euclidean')

# Convert the dist object to a matrix
dist_matrix_99 <- as.matrix(dist_16s_99)

# Ensure pH and temp are numeric
metadata_16s_99$ph <- as.numeric(metadata_16s_99$ph)
metadata_16s_99$temp <- as.numeric(metadata_16s_99$temp)

# Check overall system and season differences
permanova_16s2_99<- adonis2(dist_16s_99~system * season + ph + temp + day, data=metadata_16s, perm = 999, by = "terms")
permanova_16s2_99

adonis2(
  dist_16s_99 ~ (system + season + ph + temp + day)^3,  # includes all main effects + 2-way + 3-way interactions
  data = metadata_16s,
  permutations = 999,
  by = "terms"
)
  
#Check overall season differences
fall_data <- subset(metadata_16s_99, season == "f")
spring_data <- subset(metadata_16s_99, season == "s")

fall_dist <- dist(asv.n0.clr_16s_99[rownames(fall_data), ], method = 'euclidean')
spring_dist <- dist(asv.n0.clr_16s_99[rownames(spring_data), ], method = 'euclidean')

# For fall systems
permanova_fall <- adonis2(fall_dist ~ system, data = fall_data, perm = 999, p.adjust.m = 'bonferroni')
permanova_fall

# For spring systems
permanova_spring <- adonis2(spring_dist ~ system, data = spring_data, perm = 999, p.adjust.m = 'bonferroni')
permanova_spring

  # Loop through each season and perform PERMANOVA for system pairs 0.99 comparison
  seasons <- c("f", "s")  # Fall and Spring
  permanova_results_season_pairwise_99 <- list()  # To store results for each season
  
  # Convert distance object to matrix
  dist_matrix_99 <- as.matrix(dist_16s_99)
  
  # Loop through seasons and test each pair of systems within the season
  for (s in seasons) {
    # Subset the metadata for the current season
    metadata_subset_99 <- metadata_16s_99[metadata_16s_99$season == s, ]
    
    # Get the unique systems in this season
    systems_in_season_99 <- unique(metadata_subset_99$system)
    
    # Create an empty list to store the results for each pair
    season_results_99 <- list()
    
    # Loop through each pair of systems in the current season
    for (i in 1:(length(systems_in_season_99) - 1)) {
      for (j in (i + 1):length(systems_in_season_99)) {
        system_1 <- systems_in_season_99[i]
        system_2 <- systems_in_season_99[j]
        
        # Subset metadata for both systems in the pair
        metadata_pair_99 <- metadata_subset_99[metadata_subset_99$system %in% c(system_1, system_2), ]
        
        # Get the corresponding sample IDs for this pair of systems
        sample_ids_99 <- rownames(metadata_pair_99)
        
        # Check that the sample IDs match the row/column names of dist_matrix
        if (!all(sample_ids_99 %in% rownames(dist_matrix_99))) {
          stop("Some sample IDs from metadata do not match the row names of the distance matrix!")
        }
        
        # Subset the distance matrix by the sample IDs
        dist_subset_pair_99 <- dist_matrix_99[sample_ids_99, sample_ids_99, drop = FALSE]  # Ensure it's a matrix
        
        # Run PERMANOVA for the pair of systems in the current season
        result_99 <- adonis2(dist_subset_pair_99 ~ system, data = metadata_pair_99, permutations = 9999)
        
        # Store the results in the list
        season_results_99[[paste(system_1, "vs", system_2)]] <- result_99
      }
    }
    
    # Store the results for the current season
    permanova_results_season_pairwise_99[[s]] <- season_results_99
  }
  
  # Print the pairwise PERMANOVA results for each season
  for (s in seasons) {
    print(paste("PERMANOVA pairwise results for season", s))
    for (pair in names(permanova_results_season_pairwise_99[[s]])) {
      print(paste("Comparing:", pair))
      print(permanova_results_season_pairwise_99[[s]][[pair]])
    }
  }
  