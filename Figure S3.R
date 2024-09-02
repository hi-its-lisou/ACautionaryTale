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

