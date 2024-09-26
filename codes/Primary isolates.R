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

#primary isolates genera
# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

ps <- subset_samples(ps, sample_type %in% c("Primary_isolate"))

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

(relativeBar <- df_Genus %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sampleid , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
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
        fill = "white"),
      panel.background = element_blank()
    ))

ggsave("figures/primaryisolates.png", height=8, width=15)



#primary isolates ASVs

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps
#ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps@tax_table) 
df_temp$asvs <- row.names(ps@tax_table)
ps@tax_table <- tax_table(as.matrix(df_temp))


# Determine the number of reads per sample ####
sample_sums(ps)
ps <- subset_samples(ps, sample_type %in% c("Primary_isolate", "Honeybee_Culture", "Tosti_Culture", "Larva_Culture"))

# Sort top 20 ASVs
top20ASV = names(sort(taxa_sums(ps), TRUE)[1:40])
taxtab20 = cbind(tax_table(ps), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(ps)
                                   [top20ASV, "asvs"], "character")

tax_table(ps) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(ps, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sampleid)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")

# Colour palette to distinguish true Acinetobacter from contam
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ASVList = unique(tax_table(ps)[,"asvs"])
ASVPalette = getPalette(length(ASVList))
names(ASVPalette) = ASVList

new_color_Acinetobacter <- "#FA26A0"
Acinetobacter_to_change <- "c624f2e4228eea7296b2a77e2d4b7e50"
ASVPalette[Acinetobacter_to_change] <- new_color_Acinetobacter

new_color_Tyzz <- "#FF0000"
Tyzz_to_change <- "80626c0d45293428d118ce1f05a1ab18"
ASVPalette[Tyzz_to_change] <- new_color_Tyzz

#new_color_Erwinia <- "#3D550C"
new_color_Erwinia <- "#FF0000"
Erwinia_to_change <- "ce7f481ae002dcbd2d89b6b625a75992"
ASVPalette[Erwinia_to_change] <- new_color_Erwinia

new_color_NA_unclassified <- "#777777"
NA_unclassified_to_change <- "1a3617c309356d969a4e4d202a6591e1"
ASVPalette[NA_unclassified_to_change] <- new_color_NA_unclassified

#new_color_Enterobacteriaceae_unclassified <- "#FFFFAF"
new_color_Enterobacteriaceae_unclassified <- "#FF0000"
Enterobacteriaceae_unclassified_to_change <- "4a32f1635f88214c4874eec85df56595"
ASVPalette[Enterobacteriaceae_unclassified_to_change] <- new_color_Enterobacteriaceae_unclassified


#Plot relative abundance of Acinetobacter ASVs
(ASVplot <- df_ASV %>%
  ggplot(aes(x = sampleid, y = Abundance, fill = ASV_20)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = ASVPalette) +
  facet_nested(~ sampleid , scales = "free", space = "free") +
  labs(x = "sampleid", y = "Relative abundance") +
  theme( 
    axis.text.y = element_text(size=14, face = "bold"),
    axis.title.y = element_text(size=14, face = "bold"),
    axis.ticks.y = element_line(linewidth = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
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
      fill = "white"),
    panel.background = element_blank()
  ))

ggsave("figures/primaryisolatesASVs.png", height=15, width=20)






#primary isolates ASVs

# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps
#ps_acinetobacter <- subset_taxa(ps, Genus=="Acinetobacter")
df_temp <- as.data.frame(ps@tax_table) 
df_temp$asvs <- row.names(ps@tax_table)
ps@tax_table <- tax_table(as.matrix(df_temp))


# Determine the number of reads per sample ####
sample_sums(ps)
ps <- subset_samples(ps, sampleid %in% c("EW-185", 
"EW-186", 
"EW-187", 
"EW-188", 
"EW-189", 
"EW-194", 
"EW-195", 
"EW-196", 
"EW-197", 
"EW-200", 
"EW-201"))

# Sort top 20 ASVs
top20ASV = names(sort(taxa_sums(ps), TRUE)[1:40])
taxtab20 = cbind(tax_table(ps), ASV_20 = NA)
taxtab20[top20ASV, "ASV_20"] <- as(tax_table(ps)
                                   [top20ASV, "asvs"], "character")

tax_table(ps) <- tax_table(taxtab20)
ps_ASV_ra <- transform_sample_counts(ps, function(x) 100 * x/sum(x))
df_ASV <- psmelt(ps_ASV_ra)
df_ASV <- arrange(df_ASV, sampleid)
df_ASV$ASV_20[is.na(df_ASV$ASV_20)] <- c("Other")

# Colour palette to distinguish true Acinetobacter from contam
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
ASVList = unique(tax_table(ps)[,"asvs"])
ASVPalette = getPalette(length(ASVList))
names(ASVPalette) = ASVList

new_color_Acinetobacter <- "#FA26A0"
Acinetobacter_to_change <- "c624f2e4228eea7296b2a77e2d4b7e50"
ASVPalette[Acinetobacter_to_change] <- new_color_Acinetobacter

new_color_Tyzz <- "#FF0000"
Tyzz_to_change <- "80626c0d45293428d118ce1f05a1ab18"
ASVPalette[Tyzz_to_change] <- new_color_Tyzz

new_color_Erwinia <- "#3D550C"
Erwinia_to_change <- "ce7f481ae002dcbd2d89b6b625a75992"
ASVPalette[Erwinia_to_change] <- new_color_Erwinia

new_color_NA_unclassified <- "#777777"
NA_unclassified_to_change <- "1a3617c309356d969a4e4d202a6591e1"
ASVPalette[NA_unclassified_to_change] <- new_color_NA_unclassified

new_color_Enterobacteriaceae_unclassified <- "#FFFFAF"
Enterobacteriaceae_unclassified_to_change <- "4a32f1635f88214c4874eec85df56595"
ASVPalette[Enterobacteriaceae_unclassified_to_change] <- new_color_Enterobacteriaceae_unclassified


#Plot relative abundance of Acinetobacter ASVs
(ASVplot <- df_ASV %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = ASV_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = ASVPalette) +
    facet_nested(~ sampleid , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
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
        fill = "white"),
      panel.background = element_blank()
    ))

ggsave("figures/pollenprimaryisolatesASVs.png", height=13, width=18)




# Load phyloseq object ####
ps <- qza_to_phyloseq(
  features = "input_files/merged_table.qza",
  taxonomy = "input_files/merged_taxonomy.qza",
  tree = "input_files/merged_sepp_tree.qza",
  metadata = "input_files/combined_metadata_4.tsv")
ps@sam_data$sampleid = rownames(ps@sam_data)
ps

# Determine the number of reads per sample ####
sample_sums(ps)
ps <- subset_samples(ps, sampleid %in% c("EW-185", 
                                         "EW-186", 
                                         "EW-187", 
                                         "EW-188", 
                                         "EW-189", 
                                         "EW-194", 
                                         "EW-195", 
                                         "EW-196", 
                                         "EW-197", 
                                         "EW-200", 
                                         "EW-201"))

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

(relativeBar <- df_Genus %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sampleid , scales = "free", space = "free") +
    labs(x = "sampleid", y = "Relative abundance") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
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
        fill = "white"),
      panel.background = element_blank()
    ))

ggsave("figures/pollenprimaryisolates_genera.png", height=8, width=12)

