# Bok Choy Hydroponic Study
# Shannon Diversity Analysis at ASV level
# Last updated: 11/13/2025

# Set working directory 
setwd('/storage/home/akb6471/work/hydroponics_rerun/Reads_forSRA/')

# Load libraries 
library(vegan)
library(dplyr)
library(ggplot2)
library(viridis)

# Upload data
# Load in data
asv_16s <- read.csv("asv_16s_filtered.csv", header=TRUE)
colnames(asv_16s)[1] <- "SampleID"
rownames(asv_16s) <- asv_16s$SampleID
asv_16s <- asv_16s[, -which(colnames(asv_16s) == "SampleID")]
metadata_16s <- read.csv('metadata_16s_updated99.csv', header=TRUE, row.names=1)

# Check sequencing depths 
sample_depths <- rowSums(asv_16s)
summary(sample_depths)
table(sample_depths)               # quick frequency table
which_low <- which(sample_depths < 1500)
length(which_low)                   # how many samples < 1500
sample_depths[which_low]            # depths of low samples
names(which_low) <- NULL            

# Rarefy to 1,500 reads per sample
set.seed(42)
asv_raref <- rrarefy(asv_16s, sample = 1500)

# Calculate shannon diversity
shannon_raref <- diversity(asv_raref, index = "shannon")

# Add shannon values to metadata
metadata_16s$Shannon <- shannon_raref

# Remove day 0 samples for temporal analysis
metadata_time <- metadata_16s %>% filter(day != 0)

# Fit linear model with system, time, pH, and temperature
lm_combined <- lm(Shannon ~ system + day + ph + season + temp, data = metadata_time)
summary(lm_combined)

# ANOVA table
anova(lm_combined)

# Plot shannon diversity
# Keeping day 0 for graphing
# 1. Start from metadata_time and fix labels 
metadata_plot <- metadata_16s %>%
  mutate(
    # Recode system codes to nice labels
    system = recode(system,
                    "di"  = "DI",
                    "dwc" = "DWC",
                    "ef"  = "EF",
                    "kr"  = "KR",
                    "nf"  = "NFT"),
    system = factor(system, levels = c("DI", "DWC", "EF", "KR", "NFT")),
    
    # Recode seasons
    season = recode(season,
                    "f" = "Fall",
                    "s" = "Spring"),
    season = factor(season, levels = c("Fall", "Spring"))
  )

metadata_plot

# 2. Keep the same Viridis palette for the 5 systems ----
sys_palette <- viridis(5, option = "D", end = 0.85)
names(sys_palette) <- c("DI", "DWC", "EF", "KR", "NFT")

# 3. Plot: Shannon on Y, Day on X, facetted by Season ----
gg <- ggplot(metadata_plot,
             aes(x = day, y = Shannon, color = system, group = system)) +
  geom_point(size = 5, alpha = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
  facet_wrap(~ season) +
  scale_color_manual(values = sys_palette, name = "System") +
  theme_bw(base_size = 30) +
  labs(
    x = "Day",
    y = "Shannon Diversity"
  ) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 26),
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 34)
  )

gg

# Save ggplot
ggsave(
  filename = "fig5.png",   
  plot = gg,                    
  width = 18,                            
  height = 10,                            
  dpi = 600                              
)
