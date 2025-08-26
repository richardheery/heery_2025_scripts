# Create TMR distribution plots

# Load required packages
library(GenomicRanges)
library(dplyr)
library(ggpubr)
library(doParallel)
library(methodical)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/tmr_plot_functions.R")
source("../auxillary_scripts/granges_functions.R")

# Get low mappability regions 
low_mappability_regions = readRDS("../auxillary_data/low_mappability_regions.rds")

### Make distribution plots for unfiltered TMRs within 5 KB and 50 KB

# Load unfiltered TMR within 5 KB 
cpgea_normal_tmrs_unfiltered = readRDS("tmr_granges/cpgea_normal_tmrs_unfiltered.rds")
cpgea_tumour_tmrs_unfiltered = readRDS("tmr_granges/cpgea_tumour_tmrs_unfiltered.rds")
mcrpc_tmrs_unfiltered = readRDS("tmr_granges/mcrpc_tmrs_unfiltered.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions_unfiltered = bin_relative_tmrs(cpgea_normal_tmrs_unfiltered, width = 4750)
cpgea_tumour_5kb_tmr_distributions_unfiltered = bin_relative_tmrs(cpgea_tumour_tmrs_unfiltered, width = 4750)
mcrpc_tmr_5kb_distributions_unfiltered = bin_relative_tmrs(mcrpc_tmrs_unfiltered, width = 4750)

# Combine TMR distributions
combined_5kb_tmr_distributions_unfiltered = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions_unfiltered, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions_unfiltered,
  mcrpc = mcrpc_tmr_5kb_distributions_unfiltered, .id = "dataset"
  )

# Make a plot of the TMR distributions in the 3 data sets and save
tmrs_5kb_bins_plot_unfiltered = ggplot(combined_5kb_tmr_distributions_unfiltered, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot_unfiltered = customize_ggplot_theme(plot = tmrs_5kb_bins_plot_unfiltered, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_5kb_bins_plot_unfiltered

# # Load unfiltered TMR within 50 KB 
cpgea_normal_tmrs_50kb_unfiltered = readRDS("tmr_granges/cpgea_normal_tmrs_50kb_unfiltered.rds")
cpgea_tumour_tmrs_50kb_unfiltered = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb_unfiltered.rds")
mcrpc_tmrs_50kb_unfiltered = readRDS("tmr_granges/mcrpc_tmrs_50kb_unfiltered.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50kb
cpgea_normal_50kb_tmr_distributions_unfiltered = bin_relative_tmrs(cpgea_normal_tmrs_50kb_unfiltered, width = 49750)
cpgea_tumour_50kb_tmr_distributions_unfiltered = bin_relative_tmrs(cpgea_tumour_tmrs_50kb_unfiltered, width = 49750)
mcrpc_tmr_50kb_distributions_unfiltered = bin_relative_tmrs(mcrpc_tmrs_50kb_unfiltered, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions_unfiltered = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions_unfiltered, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions_unfiltered,
  mcrpc = mcrpc_tmr_50kb_distributions_unfiltered, .id = "dataset"
  )

# Make a plot of the TMR distributions in the 3 data sets and save
tmrs_50kb_bins_plot_unfiltered = ggplot(combined_50kb_tmr_distributions_unfiltered, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_unfiltered = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_unfiltered, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_vline(xintercept = -5000, linetype = "dotted") + geom_vline(xintercept = 5000, linetype = "dotted")
tmrs_50kb_bins_plot_unfiltered

# Combine unfiltered plots for 5 KB and 50 KB
supp_figure3 = ggpubr::ggarrange(tmrs_5kb_bins_plot_unfiltered, tmrs_50kb_bins_plot_unfiltered, labels = c("A", "B"), nrow = 2, common.legend = T, legend = "right")
ggsave(plot = supp_figure3, "../figures/supp_figure3.pdf", width = 16, height = 18)

### Make 50 KB distribution plots for TMRs overlapping poorly mappable regions

# # Load unfiltered TMR within 50 KB 
cpgea_normal_tmrs_50kb_low_mappability = subsetByOverlaps(cpgea_normal_tmrs_50kb_unfiltered, low_mappability_regions)
cpgea_tumour_tmrs_50kb_low_mappability = subsetByOverlaps(cpgea_tumour_tmrs_50kb_unfiltered, low_mappability_regions)
mcrpc_tmrs_50kb_low_mappability = subsetByOverlaps(mcrpc_tmrs_50kb_unfiltered, low_mappability_regions)

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50kb
cpgea_normal_50kb_tmr_distributions_low_mappability = bin_relative_tmrs(cpgea_normal_tmrs_50kb_low_mappability, width = 49750)
cpgea_tumour_50kb_tmr_distributions_low_mappability = bin_relative_tmrs(cpgea_tumour_tmrs_50kb_low_mappability, width = 49750)
mcrpc_tmr_50kb_distributions_low_mappability = bin_relative_tmrs(mcrpc_tmrs_50kb_low_mappability, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions_low_mappability = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions_low_mappability, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions_low_mappability,
  mcrpc = mcrpc_tmr_50kb_distributions_low_mappability, .id = "dataset"
  )

# Make a plot of the TMR distributions in the 3 data sets and save
tmrs_50kb_bins_plot_low_mappability = ggplot(combined_50kb_tmr_distributions_low_mappability, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_low_mappability = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_low_mappability, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_vline(xintercept = -5000, linetype = "dotted") + geom_vline(xintercept = 5000, linetype = "dotted")
tmrs_50kb_bins_plot_low_mappability

# Make 50KB distribution plots after removal of TMRs overlapping poorly mappable regions

# Load filtered TMRs within 50KB
cpgea_normal_tmrs_50kb = readRDS("tmr_granges/cpgea_normal_tmrs_50kb.rds")
cpgea_tumour_tmrs_50kb = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb.rds")
mcrpc_tmrs_50kb = readRDS("tmr_granges/mcrpc_tmrs_50kb.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50kb
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs_50kb, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs_50kb, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs_50kb, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distributions in the 3 data sets and save
tmrs_50kb_bins_plot = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot = customize_ggplot_theme(plot = tmrs_50kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_vline(xintercept = -5000, linetype = "dotted") + geom_vline(xintercept = 5000, linetype = "dotted")
tmrs_50kb_bins_plot

# Combine unfiltered plots for 5 KB and 50 KB
supp_figure4 = ggpubr::ggarrange(tmrs_50kb_bins_plot_low_mappability, tmrs_50kb_bins_plot, labels = c("A", "B"), nrow = 2, common.legend = T, legend = "right")
ggsave(plot = supp_figure4, "../figures/supp_figure4.pdf", width = 16, height = 18)

# Make distribution 5KB plots after removal of TMRs overlapping poorly mappable regions

# Load filtered TMRs
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions = bin_relative_tmrs(tmrs = cpgea_normal_tmrs, width = 4750)
cpgea_tumour_5kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 4750)
mcrpc_tmr_5kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 4750)

# Combine TMR distributions
combined_5kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions,
  mcrpc = mcrpc_tmr_5kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distributions in the 3 data sets and save
tmrs_5kb_bins_plot = ggplot(filter(combined_5kb_tmr_distributions, dataset == "cpgea_normal"), aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot = customize_ggplot_theme(plot = tmrs_5kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
    geom_vline(xintercept = 0, linetype = "dotted") + geom_vline(xintercept = 5000, linetype = "dotted") + geom_vline(xintercept = 5000, linetype = "dotted")
tmrs_5kb_bins_plot
saveRDS(tmrs_5kb_bins_plot, "tmrs_5kb_bins_plot.rds")

### Make introns and exons plots

tmr_list = list(
  "Normal Prostate" = cpgea_normal_tmrs,
  "Prostate Tumours" = cpgea_tumour_tmrs,
  "Prostate Metastases" = mcrpc_tmrs)

# 
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")[names(tss_gr)]

# Separate transcript bodies into 20 equally sized sections
transcripts_gr_sections = tile(transcripts_gr, n = 100)

# Reverse the GRanges for transcripts on the - strand
transcripts_gr_sections[strand(transcripts_gr) == "-"] = GRangesList(lapply(transcripts_gr_sections[strand(transcripts_gr) == "-"], rev))

# Convert transcripts_gr_sections into a flat GRanges
transcripts_gr_sections = unlist(unname(transcripts_gr_sections))

# Add transcript ID and the section number as metadata columns
transcripts_gr_sections$transcript_id = names(transcripts_gr_sections)
transcripts_gr_sections$region = paste("Region", 1:100)

# Convert transcripts_gr_sections back into a GRangesList
tss_grl = split(transcripts_gr_sections, transcripts_gr_sections$transcript_id)

# Load introns and exons for transcripts as a GRanges list
#tss_grl = readRDS("../auxillary_data/pc_transcripts_exons_and_introns_grl.rds")
tss_grl_expanded = expand_transcripts(grl = tss_grl, expand_upstream = seq(500, 5000, 500), expand_downstream = seq(500, 5000, 500))

# Create a vector with the regions to plot
regions = c(paste0("TSS-", rev(seq(500, 5000, 500))), paste("Region", 1:100), paste0("TES+", seq(500, 5000, 500)))
xlabels = rep("", length(regions))
xlabels[c(2, 6, 115, 119)] = c("-4,250", "-2,250", "+2,250", "+4,250")
xlabels[which(regions %in% c("25%", "50%", "75%"))] = c("25%", "50%", "75%")
#regions = c(paste0("TSS-", rev(seq(500, 5000, 500))), exons_introns, paste0("TES+", seq(500, 5000, 500)))

# Create absolute count and normalized counts for CPGEA normal samples
cpgea_normal_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_normal_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = F)
cpgea_normal_tmrs_introns_exons_plot
ggsave(plot = cpgea_normal_tmrs_introns_exons_plot, "cpgea_normal_tmrs_introns_exons_plot.pdf", width = 20.828, height = 7.875698)

# Load TMR example plots and comboine with combined_normal_plots to make figure 5
combined_tmr_correlation_plots = readRDS("combined_tmr_correlation_plots.rds")
figure4 = ggarrange(combined_tmr_correlation_plots, ggpubr::ggarrange(cpgea_normal_tmrs_introns_exons_plot, labels = "C"), heights = c(10.086, 3.9335), nrow = 2)
ggsave(plot = figure4, "../figures/figure4.pdf", width = 20.828, height = 28.07)

# Create absolute count and normalized counts for CPGEA tumour samples
cpgea_tumour_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_tumour_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = F)
cpgea_tumour_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_tumour_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = T)

# Create absolute count and normalized counts for MCRPC samples
mcrpc_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = mcrpc_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = F)
mcrpc_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = mcrpc_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = T)

# Combine tumour and metastases plots
tumour_and_metastases_plots = list(cpgea_tumour_tmrs_introns_exons_plot, mcrpc_tmrs_introns_exons_plot)
combined_tumour_metastases_plots = ggarrange(plotlist = tumour_and_metastases_plots, nrow = 2, ncol = 1, labels = LETTERS[1:2], common.legend = T, legend = "right")
combined_tumour_metastases_plots
ggsave(plot = combined_tumour_metastases_plots, "../figures/supp_figure10.pdf", width = 32, height = 18)

# Make distribution plots just for MANE transcripts
mane_transcripts = readRDS("../auxillary_data/mane_pc_transcript_ids.rds")
cpgea_normal_mane_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_normal_tmrs[cpgea_normal_tmrs$ID %in% mane_transcripts], 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs for MANE Transcripts in Normal Prostate", normalize = F)
cpgea_normal_mane_tmrs_introns_exons_plot
ggsave(plot = cpgea_normal_mane_tmrs_introns_exons_plot, "../figures/supp_figure12.pdf", width = 24, height = 13.5)

# Load Roadmap TMRs
roadmap_tmrs = readRDS("tmr_granges/roadmap_tmrs.rds")

# Bin relative TMRs for roadmap TMRs
roadmap_tmrs_distributions = bin_relative_tmrs(roadmap_tmrs, width = 4750)

# Make distribution plots for Roadmap TMRs
roadmap_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = F)
roadmap_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = T)

# Combine plots and save
combined_roadmap_tmrs_introns_exons_plot = ggpubr::ggarrange(roadmap_tmrs_introns_exons_plot, roadmap_tmrs_introns_exons_plot, 
  align = "hv", nrow = 2, labels = c("A", "B"), common.legend = T, legend = "right")
combined_roadmap_tmrs_introns_exons_plot
ggsave(plot = roadmap_tmrs_introns_exons_plot, "../figures/supp_figure11.pdf", width = 24, height = 13.5)
