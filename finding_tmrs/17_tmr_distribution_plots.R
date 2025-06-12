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

### Make distribution plots for unfiltered TMRs within 5 KB

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
  cpgea_normal = cpgea_normal_5kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions,
  mcrpc = mcrpc_tmr_5kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot_unfiltered = ggplot(combined_5kb_tmr_distributions_unfiltered, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot_unfiltered = customize_ggplot_theme(plot = tmrs_5kb_bins_plot_unfiltered, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_5kb_bins_plot_unfiltered

# Make distribution plots after removal of TMRs overlapping poorly mappable regions

# Load filtered TMRs
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 4750)
cpgea_tumour_5kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 4750)
mcrpc_tmr_5kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 4750)

# Combine TMR distributions
combined_5kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions,
  mcrpc = mcrpc_tmr_5kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot = ggplot(combined_5kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot = customize_ggplot_theme(plot = tmrs_5kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_5kb_bins_plot
saveRDS(tmrs_5kb_bins_plot, "tmrs_5kb_bins_plot.rds")
ggsave(plot = tmrs_5kb_bins_plot , "../figures/supplementary_figure10A.pdf",  width = 27, height = 9)

### Make distribution plots for unfiltered TMRs within 50 KB

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

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_unfiltered = ggplot(combined_50kb_tmr_distributions_unfiltered, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_unfiltered = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_unfiltered, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_50kb_bins_plot_unfiltered

# Make distribution plots after removal of TMRs overlapping poorly mappable regions

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

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot = customize_ggplot_theme(plot = tmrs_50kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_50kb_bins_plot
saveRDS(tmrs_50kb_bins_plot, "tmrs_50kb_bins_plot.rds")




### Make introns and exons plots

tmr_list = list(
  "Normal Prostate" = cpgea_normal_tmrs,
  "Prostate Tumours" = cpgea_tumour_tmrs,
  "Prostate Metastases" = mcrpc_tmrs)

# Load introns and exons for transcripts as a GRanges list
tss_grl = readRDS("../auxillary_data/pc_transcripts_exons_and_introns_grl.rds")
tss_grl_expanded = expand_transcripts(tss_grl, seq(500, 2000, 500), seq(500, 2000, 500))

# Create a vector with the regions to plot
exons = paste0("exon_", 1:10)
introns = paste0("intron_", 1:10)
exons_introns = c(rbind(exons, introns))
exons_introns = exons_introns[-length(exons_introns)]
regions = c(paste0("TSS-", rev(seq(500, 2000, 500))), exons_introns, paste0("TES+", seq(500, 2000, 500)))

# Create absolute count and normalized counts for CPGEA normal samples
cpgea_normal_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_normal_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = F)
cpgea_normal_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_normal_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = T)
combined_normal_plots = ggarrange(plotlist = list(cpgea_normal_tmrs_introns_exons_plot, cpgea_normal_tmrs_introns_exons_plot_normalized), 
  nrow = 2, labels = c("C", "D"), common.legend = T, legend = "right")

# Load TMR example plots and comboine with combined_normal_plots to make figure 5
combined_tmr_correlation_plots = readRDS("combined_tmr_correlation_plots.rds")
figure5 = ggarrange(combined_tmr_correlation_plots, combined_normal_plots, heights = c(10.086, 7.867), nrow = 2)
ggsave(plot = figure5, "../figures/figure5.pdf", width = 20.828, height = 35.945)

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
tumour_and_metastases_plots = list(cpgea_tumour_tmrs_introns_exons_plot, cpgea_tumour_tmrs_introns_exons_plot_normalized,
  mcrpc_tmrs_introns_exons_plot, mcrpc_tmrs_introns_exons_plot_normalized)
combined_tumour_metastases_plots = ggarrange(plotlist = tumour_and_metastases_plots, nrow = 2, ncol = 2, labels = LETTERS[1:4], common.legend = T, legend = "right")
combined_tumour_metastases_plots
ggsave(plot = combined_tumour_metastases_plots, "../figures/supp_figure11.pdf", width = 32, height = 18)

# Make distribution plots just for MANE transcripts
mane_transcripts = readRDS("../auxillary_data/mane_pc_transcript_ids.rds")
cpgea_normal_mane_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_normal_tmrs[cpgea_normal_tmrs$ID %in% mane_transcripts], 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs for MANE Transcripts in Normal Prostate", normalize = T)
cpgea_normal_mane_tmrs_introns_exons_plot_normalized
ggsave(plot = cpgea_normal_mane_tmrs_introns_exons_plot_normalized, "../figures/supp_figure12.pdf", width = 16, height = 9)

# Load Roadmap TMRs without repeats
roadmap_tmrs = readRDS("tmr_granges/roadmap_tmrs.rds")

# Bin relative TMRs for roadmap TMRs
roadmap_tmrs_distributions = bin_relative_tmrs(roadmap_tmrs, width = 4750)

# Make distribution plots for Roadmap TMRs
roadmap_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = F)
roadmap_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = T)

# Combine plots and save
combined_roadmap_tmrs_introns_exons_plot = ggpubr::ggarrange(roadmap_tmrs_introns_exons_plot, roadmap_tmrs_introns_exons_plot_normalized, 
  align = "hv", nrow = 2, labels = c("A", "B"), common.legend = T, legend = "right")
ggsave(plot = combined_roadmap_tmrs_introns_exons_plot, "../figures/supp_figure12.pdf", width = 16, height = 9)
