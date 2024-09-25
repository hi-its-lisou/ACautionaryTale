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
ps <- subset_samples(ps, sample_type %in% c("Negative_control", "Food", "Prepupae", "Honey_bee"))

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
                  "Acinetobacter", 
                  "Erwinia", 
                  "Tyzzerella", 
                  "NA_unclassified", 
                  "Staphylococcus",
                  "Enterobacteriaceae_unclassified", 
                  "Aquabacterium", 
                  "Gilliamella", 
                  "Pseudomonas", 
                  "Commensalibacter",
                  "Lactobacillus",
                  "Snodgrassella",
                  "Orbaceae_unclassified",
                  "Bartonella", 
                  "Sphingobium", 
                  "Bifidobacterium",
                  "Enhydrobacter",
                  "Other")



df_Genus$Genus_20 <- factor(df_Genus$Genus_20, levels = custom_order)

df_Genus <- df_Genus %>%
  dplyr::mutate(sample_type2 = dplyr::if_else(sample_type == "Prepupae",
                                              stringr::str_c("ABexp_", AB_treatment),
                                              sample_type), .after = sample_type)



order <- c("Negative_control", "ABexp_control", "ABexp_antibiotic", "Food", "Honey_bee")
df_Genus$sample_type2 <- factor(df_Genus$sample_type2, levels=order)


(absolutePlot <- df_Genus %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    dplyr::distinct(sampleid, Genus_20, .keep_all = TRUE) %>%
    dplyr::mutate(logDNA_weighted = logDNA*(Abundance/100)) %>%
    ggplot(aes(x = sampleid)) +
    geom_col(width = 1, aes(fill = Genus_20, y = logDNA_weighted)) +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2, scales = "free", space = "free") +
    labs(x = "sampleid", y = "Absolute abundance (ng)") +
    theme( 
      axis.text.y = element_text(size=14, face = "bold"),
      axis.title.y = element_text(size=14, face = "bold"),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "none",
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
    ))

(relativeBar <- df_Genus %>%
    dplyr::filter(complete.cases(logDNA)) %>%
    ggplot(aes(x = sampleid, y = Abundance, fill = Genus_20)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = my_palette) +
    facet_nested(~ sample_type2 , scales = "free", space = "free") +
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
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.background = element_blank()
    ))



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
      axis.text.y = element_text(size=12, face =, 'bold'),
      axis.title.y = element_text(size=12, face = 'bold'),
      axis.ticks.y = element_line(linewidth = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
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
    ))

legend <- get_legend(legendPlot)
(abundancesplot_nolegend <- plot_grid(absolutePlot, relativeBar, ncol = 1, align = "hv"))

# Now add the legend to the final combined plot
(Abundances_plot <- cowplot::plot_grid(
  abundancesplot_nolegend,
  legend,  # Add the extracted legend as a separate plot
  ncol = 1,
  rel_heights = c(1.3, 0.5)
))


(Abundances_plot <- cowplot::plot_grid(abundancesplot_nolegend,
                                       cowplot::get_legend(legendPlot),
                                       cols = 1,
                                       rel_heights = c(1.3, 0.2)
))


ggsave("figures/Abundances_plot_nolegend.png", height=8, width=12)
