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

