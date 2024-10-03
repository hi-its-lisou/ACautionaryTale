# Load the required libraries ####
library(phyloseq)
library(tidyverse)
library(qiime2R)
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
library(cowplot)
library(gridGraphics)


# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Load custom colour palette ####
colours_df <- read_csv("input_files/colour_list.csv")
my_palette <- colours_df$colours
names(my_palette) <- colours_df$genera
my_palette

# Identify zero abundance samples
zero_abundance_samples <- sample_sums(ps) == 0
ps <- subset_samples(ps, !zero_abundance_samples)
ps

# Remove primary isolates
ps <- subset_samples(ps, sample_type %in% c("Prepupae"))

#Collapse to Genus level and pull our relative abundance of top 20 common genuses.
ps_Genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

if ("Family" %in% colnames(tax_table(ps_Genus))) {
  
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Family"]), "unclassified", sep = "_")
} else {
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- "Unclassified"
}

# Replace "NA" (as a string) in the Genus column with "NA_unclassified"
tax_table(ps_Genus)[tax_table(ps_Genus)[, "Genus"] == "NA", "Genus"] <- "NA_unclassified"
tax_table(ps_Genus)[tax_table(ps_Genus)[, "Genus"] == "<NA>", "Genus"] <- "NA_unclassified"

top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")

tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus$Genus_20))

# % of reads that make up the top 20 genera ####
mean(sample_sums(prune_taxa(top20Genus, ps_Genus_ra)))

custom_order <- c("Chloroplast",
                  "Mitochondria",
                  "Tyzzerella",
                  "Enterobacteriaceae_unclassified",
                  "Other",
                  "Pseudomonas",
                  "Hydrogenophilus",
                  "Acinetobacter",
                  "Aliterella",
                  "Sphingobium",
                  "Enhydrobacter",
                  "Neisseria",
                  "Staphylococcus",
                  "Chryseobacterium",
                  "Comamonadaceae_unclassified",
                  "Curvibacter",
                  "Arsenophonus",
                  "Aquabacterium",
                  "Paenibacillus",
                  "Bacillus")
                   

(unfiltered_prepupae_Bar <- df_Genus %>% 
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ AB_treatment, scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_text(size=5),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size=14, face = "bold"),
      panel.background = element_blank()
    ))


ggsave("figures/Prepupae_unfiltered_nolegend.png", height=6, width=12)



# Remove contaminants, chloroplasts, and mitochondria ####
contam <- read_delim("input_files/chloro_mito_decontam_asvs.txt", 
                     delim = "\n", 
                     col_names = "asv")
chloro_mito_decontam_asvs <- contam$asv
all_asvs <- taxa_names(ps)
asvs_to_keep <- all_asvs[!(all_asvs %in% chloro_mito_decontam_asvs)]
ps <- prune_taxa(asvs_to_keep, ps)
ps




#Collapse to Genus level and pull our relative abundance of top 20 common genuses.
ps_Genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)

if ("Family" %in% colnames(tax_table(ps_Genus))) {
  
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- 
    paste(as.character(tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Family"]), "unclassified", sep = "_")
} else {
  tax_table(ps_Genus)[is.na(tax_table(ps_Genus)[, "Genus"]), "Genus"] <- "Unclassified"
}

# Replace "NA" (as a string) in the Genus column with "NA_unclassified"
tax_table(ps_Genus)[tax_table(ps_Genus)[, "Genus"] == "NA", "Genus"] <- "NA_unclassified"
tax_table(ps_Genus)[tax_table(ps_Genus)[, "Genus"] == "<NA>", "Genus"] <- "NA_unclassified"

top20Genus = names(sort(taxa_sums(ps_Genus), TRUE)[1:19])
taxtab20 = cbind(tax_table(ps_Genus), Genus_20 = NA)
taxtab20[top20Genus, "Genus_20"] <- as(tax_table(ps_Genus)
                                       [top20Genus, "Genus"], "character")

tax_table(ps_Genus) <- tax_table(taxtab20)
ps_Genus_ra <- transform_sample_counts(ps_Genus, function(x) 100 * x/sum(x))
df_Genus <- psmelt(ps_Genus_ra)
df_Genus <- arrange(df_Genus, sample_type)
df_Genus$Genus_20[is.na(df_Genus$Genus_20)] <- c("Other")

# View the result
print(unique(df_Genus$Genus_20))

# % of reads that make up the top 20 genera ####
mean(sample_sums(prune_taxa(top20Genus, ps_Genus_ra)))


custom_order <- c("Tyzzerella",
                  "Enterobacteriaceae_unclassified",
                  "Acinetobacter",
                  "Other",
                  "Hydrogenophilus",
                  "Enhydrobacter",
                  "Aliterella",
                  "Arsenophonus",
                  "Denitratisoma",
                  "Bacillus",
                  "Paenibacillus",
                  "Chryseobacterium",
                  "Comamonadaceae_unclassified",
                  "Pseudomonas",
                  "Pantoea",
                  "Curvibacter",
                  "Blastococcus",
                  "Aquabacterium",
                  "Burkholderia-Caballeronia-Paraburkholderia",
                  "Sodalis")

df_Genus$Genus_20 <- factor(df_Genus$Genus_20, levels = custom_order)

(filtered_prepupae_Bar <- df_Genus %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ AB_treatment, scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_text(size=8),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(size=14, face = "bold"),
      panel.background = element_blank()
    ))


ggsave("figures/Prepupae_filtered_nolegend.png", height=6, width=12)
