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


