# Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Subset the data for negative controls and adults in the acquisition experiment ####
supp_adults2023 <- subset_samples(ps, sample_type %in% c("Adults", "Negative_control"))
supp_adults2023

# Obtain top 20 genera ####
ps_Genus <- tax_glom(supp_adults2023, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Plot the relative abundance ####
(supp_RAAE <- df_Genus %>%
   mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
   ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
   geom_bar(width = 1, stat = "identity") +
   scale_fill_manual(values = my_palette) +
   facet_nested(~ sample_type + Env_exposure + sampleid, scales = "free", space = "free") +
   labs(x = "sample_type", y = "Relative abundance") +
   theme( 
     axis.text.y = element_text(size=16, face = 'bold'),
     axis.title.y = element_text(size=16, face = 'bold'),
     axis.ticks.y = element_line(linewidth = 1),
     axis.ticks.x = element_blank(),
     axis.text.x = element_blank(),
     axis.title.x = element_blank(),
     legend.title = element_blank(),
     legend.text = element_text(size = 18),
     legend.position = "bottom",
     strip.background = element_blank(),
     strip.text = element_textbox_simple(
       padding = margin(5, 0, 5, 0),
       margin = margin(5, 5, 5, 5),
       size = 10,
       face = "bold",
       halign = 0.5,
       fill = "white",
       box.color = "grey",
       linewidth = 1.5,
       linetype = "solid",),
     panel.background = element_blank()
   ))

#Save plot
ggsave("Supp1.png", height=15, width=25)

91 changes: 91 additions & 0 deletions91  
Figure S2.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,91 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

### Acinetobacter plot - Supplementary Figure 2 ###
# Load phyloseq object with only Acinetobacter ASVs
ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps_acinetobacter@tax_table) 
df_temp$asvs <- row.names(ps_acinetobacter@tax_table)
ps_acinetobacter@tax_table <- tax_table(as.matrix(df_temp))

# Subset samples
supp_adults2023 <- subset_samples(ps_acinetobacter, sample_type %in% c("Adults", "Negative_control"))
supp_adults2023

# Determine the number of reads per sample ####
sample_sums(supp_adults2023)

# Sort top 20 ASVs
top20ASV = names(sort(taxa_sums(supp_adults2023), TRUE)[1:20])
taxtab20 = cbind(tax_table(supp_adults2023), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(supp_adults2023)
                                   [top20ASV, "asvs"], "character")

tax_table(supp_adults2023) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(supp_adults2023, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sample_type)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")

# Colour palette to distinguish true Acinetobacter from contam
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ASVList = unique(tax_table(ps_acinetobacter)[,"asvs"])
ASVPalette = getPalette(length(ASVList))
names(ASVPalette) = ASVList
new_color_Acinetobacter <- "#FA26A0"
Acinetobacter_to_change <- "c624f2e4228eea7296b2a77e2d4b7e50"
ASVPalette[Acinetobacter_to_change] <- new_color_Acinetobacter

#Plot relative abundance of Acinetobacter ASVs
Acinetobacter_plot_supp2 <- df_ASV %>%
  filter(Abundance > 0) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = ASV_20)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ASVPalette) +
  facet_nested(~ sample_type + Env_exposure + Sample, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face = 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 16,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid",),
    panel.background = element_blank())

Acinetobacter_plot_supp2

# Save plot
ggsave("FigureS2.png", height=15, width=25)
91 changes: 91 additions & 0 deletions91  
Figure S3.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,91 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load phyloseq object with only Acinetobacter ASVs
ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps_acinetobacter@tax_table) 
df_temp$asvs <- row.names(ps_acinetobacter@tax_table)
ps_acinetobacter@tax_table <- tax_table(as.matrix(df_temp))

# Subset samples
AB_exp <- subset_samples(ps_acinetobacter, AB_treatment %in% c("control", "antibiotic"))
Negs <- subset_samples (ps_acinetobacter, sample_type %in% c("Negative_control", "Food"))
AB_exp_Acine <- merge_phyloseq(AB_exp, Negs) 

# Determine the number of reads per sample ####
sample_sums(AB_exp_Acine)

# Sort top 20 ASVs
top20ASV = names(sort(taxa_sums(AB_exp_Acine), TRUE)[1:20])
taxtab20 = cbind(tax_table(AB_exp_Acine), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(AB_exp_Acine)
                                   [top20ASV, "asvs"], "character")

tax_table(AB_exp_Acine) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(AB_exp_Acine, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sample_type)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")

# Colour palette to distinguish true Acinetobacter from contam
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ASVList = unique(tax_table(ps_acinetobacter)[,"asvs"])
ASVPalette = getPalette(length(ASVList))
names(ASVPalette) = ASVList
new_color_Acinetobacter <- "#FA26A0"
Acinetobacter_to_change <- "c624f2e4228eea7296b2a77e2d4b7e50"
ASVPalette[Acinetobacter_to_change] <- new_color_Acinetobacter

#Plot relative abundance of Acinetobacter ASVs
Acinetobacter_plot_supp3 <- df_ASV %>%
  filter(Abundance > 0) %>%
  ggplot(aes(x = sample_type, y = Abundance, fill = ASV_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = ASVPalette) +
  facet_nested(~ sample_type + AB_treatment + sampleid, scales = "free", space = "free") +
  labs(x = "sample_type", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=16, face = 'bold'),
    axis.title.y = element_text(size=16, face = 'bold'),
    axis.ticks.y = element_line(linewidth = 1),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_textbox_simple(
      padding = margin(5, 0, 5, 0),
      margin = margin(5, 5, 5, 5),
      size = 16,
      face = "bold",
      halign = 0.5,
      fill = "white",
      box.color = "grey",
      linewidth = 1.5,
      linetype = "solid",),
    panel.background = element_blank()
  )

#Save plot
ggsave("Acinetobacter_supp_3.png", height = 15, width = 30)

116 changes: 116 additions & 0 deletions116  
Figure3.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,116 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Rarefy ####
ps_rar <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 1337)

# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps_rar)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_rar <- prune_taxa(asvs_to_keep, ps_rar)
ps_rar

#Subset prepupae from antibiotic experiment only ####
AB_exp <- subset_samples(ps_rar, AB_treatment %in% c("control", "antibiotic"))
AB_exp #821 taxa and 39 samples

# Remove 0 abundsance samples ####
zero_abundance_samples <- sample_sums(AB_exp) == 0
filtered_physeq <- subset_samples(AB_exp, !zero_abundance_samples)
filtered_physeq

# Alpha diversity box plot ####
alpha_diversity <- alpha(filtered_physeq, index = "Shannon")
metadata <- meta(filtered_physeq)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")

# Plot alpha diversity ####
fig5 <- alpha_diversity_metadata %>%
  ggplot(aes(x = AB_treatment, y = diversity_shannon, colour = AB_treatment)) +
  geom_boxplot() +
  geom_point(size = 3)

# Obtain top 20 genera ####
ps_Genus <- tax_glom(filtered_physeq, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# % of reads that make up the top 20 genera ####
mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)

# Plot the relative abundance ####
(fig4 <- df_Genus %>%
   mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
   ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
   geom_bar(width = 1, stat = "identity") +
   scale_fill_manual(values = my_palette) +
   facet_nested(~ AB_treatment + sampleid + Year, scales = "free", space = "free") +
   labs(x = "sample_type", y = "Relative abundance") +
   theme( 
     axis.text.y = element_text(size=16, face = 'bold'),
     axis.title.y = element_text(size=16, face = 'bold'),
     axis.ticks.y = element_line(linewidth = 1),
     axis.ticks.x = element_blank(),
     axis.text.x = element_blank(),
     axis.title.x = element_blank(),
     legend.title = element_blank(),
     legend.text = element_text(size = 12),
     legend.position = "bottom",
     strip.background = element_blank(),
     strip.text = element_textbox_simple(
       padding = margin(5, 0, 5, 0),
       margin = margin(5, 5, 5, 5),
       size = 10,
       face = "bold",
       halign = 0.5,
       fill = "white",
       box.color = "grey",
       linewidth = 1.5,
       linetype = "solid",),
     panel.background = element_blank()
   ))
143 changes: 143 additions & 0 deletions143  
Figure4.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,143 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps #1089 taxa and 208 samples

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Rarefy ####
ps_rar <- rarefy_even_depth(ps, sample.size = 1500, rngseed = 1337)
ps_rar 

# Filter out 0 abundance reads ####
zero_abundance_samples <- sample_sums(ps_rar) == 0
ps_rar <- subset_samples(ps_rar, !zero_abundance_samples)
ps_rar

# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")

chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps_rar)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_rar <- prune_taxa(asvs_to_keep, ps_rar)
ps_rar

# Subset the data for adults in the acquisition experiment ####
adults2023 <- subset_samples(ps_rar, sample_type %in% c("Adults"))
adults2023

# Alpha diversity of pollen provisions as they are consumed ####
alpha_diversity <- alpha(adults2023, index = "Shannon")
metadata <- meta(adults2023)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadata <- merge(alpha_diversity, metadata, by = "name")
alpha_diversity_metadata$sample_type <- factor(alpha_diversity_metadata$sample_type)
alpha_diversity_metadata$Env_exposure <- factor(alpha_diversity_metadata$Env_exposure)

# Box plot for alpha diversity
ggplot(alpha_diversity_metadata, aes(x = sample_type, y = diversity_shannon, fill = Env_exposure)) +
  geom_violin() +
  geom_point(position = position_jitterdodge(), size = 3)

# Beta diversity - PCOA ordination of unweighted UniFrac distances ####
unweighted_unifrac <- ordinate(adults2023, method = "PCoA", distance = "unifrac", weighted=F)

#Axes 1/2 ---- variation, 
(p1 <- plot_ordination(physeq = adults2023, 
                       ordination = unweighted_unifrac, 
                       color = "Env_exposure",
                       shape = "sample_type",
                       axes = c(1, 2)) +
    theme_minimal() +
    geom_point(size = 5, alpha = 0.6))

# Obtain top 20 genera ####
ps_Genus <- tax_glom(adults2023, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")
tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# % of reads that make up the top 20 genera ####
mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)

# Plot the relative abundance ####
(RAAE <- df_Genus %>%
   mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
   ggplot(aes(x = sample_type, y = Abundance, fill = Genus_20)) +
   geom_bar(width = 1, stat = "identity") +
   scale_fill_manual(values = GenusPalette) +
   facet_nested(~ sample_type + Env_exposure + sampleid, scales = "free", space = "free") +
   labs(x = "sample_type", y = "Relative abundance") +
   theme( 
     axis.text.y = element_text(size=16, face = 'bold'),
     axis.title.y = element_text(size=16, face = 'bold'),
     axis.ticks.y = element_line(linewidth = 1),
     axis.ticks.x = element_blank(),
     axis.text.x = element_blank(),
     axis.title.x = element_blank(),
     legend.title = element_blank(),
     legend.text = element_text(size = 18),
     legend.position = "bottom",
     strip.background = element_blank(),
     strip.text = element_textbox_simple(
       padding = margin(5, 0, 5, 0),
       margin = margin(5, 5, 5, 5),
       size = 10,
       face = "bold",
       halign = 0.5,
       fill = "white",
       box.color = "grey",
       linewidth = 1.5,
       linetype = "solid",),
     panel.background = element_blank()
   ))

# Combine the plots ####
AE_plot <- cowplot::plot_grid(RAAE,
                              p1,
                              ncol = 2,
                              rel_heights = c(1,0.6),
                              rel_widths = c(1, 0.5))

# Save plot ####
ggsave("Figure4.png", height=15, width=25)


88 changes: 88 additions & 0 deletions88  
Figure5.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,88 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(extrafont)
library(forcats)
library(ggforce)
library(vegan)
library(patchwork)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Identify zero abundance samples
zero_abundance_samples <- sample_sums(ps) == 0
ps <- subset_samples(ps, !zero_abundance_samples)
ps

# Subset the data for samples that underwent qPCRs (have CT scores)
qPCRs <- subset_samples(ps, complete.cases(AVG_CTscore) &
                          sample_type %in% c("Negative_control", "Food", "Adults", "Prepupae", "Honey_bee"))

# Make the phyloseq into a dataframe
sample_data_df <- as.data.frame(sample_data(qPCRs))

# Custom order for figure
custom_order <- c("Negative_control", "Adults", "Prepupae", "Food", "Honey_bee")
sample_data_df$sample_type <- factor(sample_data_df$sample_type, levels = custom_order)

# Custom colors
custom_colors <- c("control" = "#34C9CD", "antibiotic" = "#F87970", "NA" = "black")

# Convert NA values to strings so that they are still included in the plot
sample_data_df$AB_treatment <- as.character(sample_data_df$AB_treatment)
sample_data_df$AB_treatment[is.na(sample_data_df$AB_treatment)] <- "NA"
sample_data_df$AB_treatment <- factor(sample_data_df$AB_treatment, levels = c("control", "antibiotic", "NA"))

sample_data_df$Env_exposure <- as.character(sample_data_df$Env_exposure)
sample_data_df$Env_exposure[is.na(sample_data_df$Env_exposure)] <- "NA"
sample_data_df$Env_exposure <- factor(sample_data_df$Env_exposure, levels = c("Free-flying", "Nest_only", "None", "NA"))

# Create the ggplot
(qPCRPlot <- ggplot(sample_data_df, aes(x = sampleid, y = logDNA)) +
    geom_hline(aes(linetype = "limit of detection", yintercept = 2.339477 ), color = "green") +
    geom_point(aes(
      color = AB_treatment,
      shape = Env_exposure
    ), size = 3, na.rm = TRUE) +
    scale_shape_manual(values = c(
      "Free-flying" = 1,
      "Nest_only" = 4,
      "None" = 3,
      "NA" = 16
    )) +
    scale_color_manual(values = custom_colors) +
    scale_y_continuous(name = "logDNA") +
    theme(axis.text.x = element_blank(),
          strip.text = element_text(size = 12)) +
    facet_wrap(~sample_type, ncol = 5, scale = "free_x"))


# Save the figure 
ggsave("figures/Figure5.png", height=10, width=15)


192 changes: 192 additions & 0 deletions192  
Figure6.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,192 @@
  # Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(extrafont)
library(forcats)
library(ggforce)
library(vegan)
library(patchwork)

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load custom colour palette ####
colours_df <- file.path(
  "D:/Research/PhD/Manuscripts/Chapter 2/Git/2024_Chapter2_tosticauda_development_microbiome/input_files/",
  "colour_list.csv"
) %>% read_csv

my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Identify zero abundance samples
zero_abundance_samples <- sample_sums(ps) == 0
ps <- subset_samples(ps, !zero_abundance_samples)
ps

# Remove primary isolates
ps <- subset_samples(ps, !sample_type %in% c("Primary_isolate"))

#Collapse to Genus level and pull our relative abundance of top 20 common genuses.
ps_Genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:20])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")

tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")
mean(
  sample_sums(
    prune_taxa(top20Genus, ps_Genus_ra)
  )
)

custom_order <- c("Mitochondria", "Chloroplast", "Pseudomonas", "Lactobacillus", "Staphylococcus", "Escherichia-Shigella", "Other", "Rhodococcus","Brevibacterium", "Bacillus", "Snodgrassella",  "Bartonella", "Erwinia", "Curtobacterium","Streptomyces", "Sphingobium", "Gilliamella", "Tyzzerella", "Acinetobacter")
df_Genus$Genus_20 <- factor(df_Genus$Genus_20, levels = custom_order)

df_Genus <- df_Genus %>%
  dplyr::mutate(sample_type2 = dplyr::if_else(sample_type == "Adults",
                                              stringr::str_c("Tosti_", Env_exposure),
                                              dplyr::if_else(sample_type == "Prepupae",
                                                             stringr::str_c("ABexp_", AB_treatment),
                                                             sample_type)),
                .after = sample_type)



order <- c("Negative_control", "Tosti_Free-flying", "Tosti_Nest_only", "Tosti_None", "ABexp_control", "ABexp_antibiotic", "Food", "Honey_bee")
df_Genus$sample_type2 <- factor(df_Genus$sample_type2, levels=order)


(relativeBar <- df_Genus %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    #mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2 , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=16, face = 'bold'),
      axis.title.y = element_text(size=16, face = 'bold'),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      #  element_textbox_simple(
      # padding = margin(5, 0, 5, 0),
      # margin = margin(5, 5, 5, 5),
      # size = 10,
      # face = "bold",
      # halign = 0.5,
      # fill = "white",
      # box.color = "grey",
      # linewidth = 1.5,
      # linetype = "solid",),
      panel.background = element_blank()
    )
)


(absolutePlot <- df_Genus %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    dplyr::distinct(sampleid, Genus_20, .keep_all = TRUE) %>%
    #mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
    dplyr::mutate(logDNA_weighted = logDNA*(Abundance/100)) %>%
    ggplot(aes(x = sampleid)) +
    geom_col(width = 1, aes(fill = Genus_20, y = logDNA_weighted)) +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2, scales = "free", space = "free") +
    #facet_nested(~ sample_type , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Absolute abundance (ng)") +
    theme( 
      axis.text.y = element_text(size=16, face =, 'bold'),
      axis.title.y = element_text(size=16, face = 'bold'),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_textbox_simple(
        padding = margin(5, 0, 5, 0),
        margin = margin(5, 5, 5, 5),
        size = 10,
        face = "bold",
        halign = 0.5,
        fill = "white",
        box.color = "grey",
        linewidth = 1.5,
        linetype = "solid",),
      panel.background = element_blank()
    )
)

(legendPlot <- df_Genus %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    dplyr::distinct(sampleid, Genus_20, .keep_all = TRUE) %>%
    mutate(Genus_20 = reorder(Genus_20, -Abundance)) %>%
    dplyr::mutate(logDNA_weighted = logDNA*(Abundance/100)) %>%
    ggplot(aes(x = sampleid)) +
    geom_col(width = 1, aes(fill = Genus_20, y = logDNA_weighted)) +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2 , scales = "free", space = "free") +
    #facet_nested(~ sample_type , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=16, face =, 'bold'),
      axis.title.y = element_text(size=16, face = 'bold'),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 18),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_textbox_simple(
        padding = margin(5, 0, 5, 0),
        margin = margin(5, 5, 5, 5),
        size = 10,
        face = "bold",
        halign = 0.5,
        fill = "white",
        box.color = "grey",
        linewidth = 1.5,
        linetype = "solid",),
      panel.background = element_blank()
    )
)


(Abundances_plot <- cowplot::plot_grid(absolutePlot,
                                       relativeBar,
                                       cowplot::get_legend(legendPlot),
                                       cols = 1,
                                       rel_heights = c(1,0.8,0.3)))

ggsave("figures/Figure6.png", height=10, width=15)

225 changes: 225 additions & 0 deletions225  
Table 1.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,225 @@
  #load libraries
  library(phyloseq)
library(tidyverse)
library(qiime2R)
library(decontam)
library(microbiome)
library(ggh4x)
library(ggtext)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)

# Load in sequencing data
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

################################################################################
#####          calculating the % of chloroplasts and mitochondria           ####
################################################################################

# Function to computer percentage
compute_percentage <- function(taxon, ps) {
  counts_data <- as.data.frame(otu_table(ps))  #Extract the counts data from the phyloseq object
  taxon_row <- which(tax_table(ps)[, "Genus"] == taxon)  #Identify the row corresponding to the specific taxon in the taxonomic table
  taxon_counts <- sum(counts_data[taxon_row, ])  #Sum the counts of the specific taxon across all samples
  total_counts <- sum(counts_data)  #Calculate the percentage of reads attributed to the specific taxon
  percentage_taxon <- (taxon_counts / total_counts) * 100
  return(percentage_taxon)
}

# Sublist sample types
Prepupae <- subset_samples(ps, sample_type %in% c("Prepupae")& !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))
Larva <- subset_samples(ps, sample_type %in% c("Larvae"))
Honey_bee <- subset_samples(ps, sample_type %in% c("Honey_bee"))
Food <- subset_samples(ps, sample_type %in% c("Food"))
AdultsFree <- subset_samples(ps, Env_exposure %in% c("Free-flying")) #seperate only natural adult samples

antibiotic <- subset_samples(ps, AB_treatment %in% c("antibiotic"))
control <- subset_samples(ps, AB_treatment %in% c("control"))
AdultsNest <- subset_samples(ps, Env_exposure %in% c("Nest_only"))
AdultsNone <- subset_samples(ps, Env_exposure %in% c("None"))

# Run function for each sublisted sample type
# For Chloroplast
percentage_Chloroplast_Adults <- compute_percentage("Chloroplast", AdultsFree) 
percentage_Chloroplast_Larva <- compute_percentage("Chloroplast", Larva)
percentage_Chloroplast_Prepupae <- compute_percentage("Chloroplast", Prepupae)
percentage_Chloroplast_Food <- compute_percentage("Chloroplast", Food)
percentage_Chloroplast_HB <- compute_percentage("Chloroplast", Honey_bee)

percentage_Chloroplast_Adults
percentage_Chloroplast_Food
percentage_Chloroplast_Larva
percentage_Chloroplast_Prepupae
percentage_Chloroplast_HB

percentage_Chloroplast_antibiotic <- compute_percentage("Chloroplast", antibiotic) 
percentage_Chloroplast_control <- compute_percentage("Chloroplast", control) 
percentage_Chloroplast_AdultsNest <- compute_percentage("Chloroplast", AdultsNest)
percentage_Chloroplast_AdultsNone <- compute_percentage("Chloroplast", AdultsNone)

percentage_Chloroplast_antibiotic
percentage_Chloroplast_control
percentage_Chloroplast_AdultsNest
percentage_Chloroplast_AdultsNone

# For Mitochondria
percentage_Mitochondria_Adults <- compute_percentage("Mitochondria", AdultsFree) 
percentage_Mitochondria_Larva <- compute_percentage("Mitochondria", Larva)
percentage_Mitochondria_Prepupae <- compute_percentage("Mitochondria", Prepupae)
percentage_Mitochondria_Food <- compute_percentage("Mitochondria", Food)
percentage_Mitochondria_HB <- compute_percentage("Mitochondria", Honey_bee)

percentage_Mitochondria_Adults
percentage_Mitochondria_Food
percentage_Mitochondria_Larva
percentage_Mitochondria_Prepupae
percentage_Mitochondria_HB

percentage_Mitochondria_antibiotic <- compute_percentage("Mitochondria", antibiotic) 
percentage_Mitochondria_control <- compute_percentage("Mitochondria", control) 
percentage_Mitochondria_AdultsNest <- compute_percentage("Mitochondria", AdultsNest)
percentage_Mitochondria_AdultsNone <- compute_percentage("Mitochondria", AdultsNone)

percentage_Mitochondria_antibiotic
percentage_Mitochondria_control
percentage_Mitochondria_AdultsNest
percentage_Mitochondria_AdultsNone

### Calculating the proportion of putitive contaminate asvs (called using decontam at a threshold of 0.5)
#vector for contaminant asvs
contam_asvs <- read_delim("input_files/decontam_asvs.txt", 
                          delim = "\n", 
                          col_names = "asv")
contam_asvs <- contam_asvs$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% contam_asvs)]
ps_no_contam <- prune_taxa(asvs_to_keep, ps)
ps_no_contam #pyloseq object without any contaminant asvs

#Adult tostis contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying")))))) #Specifically for free-fying the adults
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Pollen provision contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Food")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Larvae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Larvae")))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic")))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Prepupae" & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Honey bee contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) 
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, sample_type == "Honey_bee"))))) 
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

# Antibiotic prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("antibiotic"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, AB_treatment %in% c("antibiotic"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

# Control prepupae contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("control"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, AB_treatment %in% c("control"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Nest eclosed bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

#Controlled eclosion bees contamination %
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_no_contam <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_no_contam, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_no_contam / total_reads_ps) * 100
100 - percentage_remaining

### Working biomass after removing ALL non-target reads ###
# Remove contaminants, chloroplasts, and mitochondria ####
chloro_mito_decontam_asvs <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                                        delim = "\n", 
                                        col_names = "asv")
chloro_mito_decontam_asvs <- chloro_mito_decontam_asvs$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps_filtered <- prune_taxa(asvs_to_keep, ps)
ps_filtered #phyloseq object with all off-target reads removed

# Free flying adults working biomass 
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Free-flying"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered,Env_exposure %in% c("Free-flying"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Larvae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Larvae")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Larvae"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Prepupae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Prepupae") & !(sample_data(ps)$AB_treatment %in% c("control", "antibiotic"))))))) 
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Prepupae" & !(sample_data(ps_filtered)$AB_treatment %in% c("control", "antibiotic")))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#Food working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Food")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Food"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

#honey bee working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, sample_type == c("Honey_bee")))))) #replace with sample_type in question
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, sample_type == "Honey_bee"))))) #replace with sample_type in question
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# antibiotic prepupae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("antibiotic"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, AB_treatment %in% c("antibiotic"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# control prepupae working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, AB_treatment %in% c("control"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, AB_treatment %in% c("control"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# nest eclosed adults working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("Nest_only"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("Nest_only"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

# controlled eclosion adults working biomass
total_reads_ps <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps, Env_exposure %in% c("None"))))))
total_reads_ps_filtered <- sum(rowSums(as.data.frame(otu_table(subset_samples(ps_filtered, Env_exposure %in% c("None"))))))
percentage_remaining <- (total_reads_ps_filtered / total_reads_ps) * 100
percentage_remaining

109 changes: 109 additions & 0 deletions109  
decontam.R
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,109 @@
