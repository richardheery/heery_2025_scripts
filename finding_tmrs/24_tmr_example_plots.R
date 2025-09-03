# Create plot of Methodical workflow and show examples of some TMRs

# Load required packages
library(methodical)
library(cowplot)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")

# Load CPGEA normal correlation results
cpgea_normal_correlation_results = readRDS("meth_transcript_cors/cpgea_normal_whole_gene_body_correlations.rds")

# Get GRanges for protein-coding TSS sites
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")

# Gt GRanges for TMRs for normal prostate samples
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs.rds")

# Show TMR workflow using qsox1 (ENST00000367602) as an example

# Get CpG correlation values for qsox1 and remove missing values
qsox1_cpg_correlations = cpgea_normal_correlation_results[["ENST00000367602"]]
qsox1_cpg_correlations = qsox1_cpg_correlations[complete.cases(qsox1_cpg_correlations), ]
row.names(qsox1_cpg_correlations) = NULL

# Get qsox1 TSS site
qsox1_tss = tss_gr[tss_gr$ID == "ENST00000367602"]

# Get TMRs for qsox1
qsox1_tmrs = cpgea_normal_tmrs[cpgea_normal_tmrs$ID == "ENST00000367602"]

# Create qsox1 CpG correlation plot
qsox1_cpg_correlation_plot = plotMethSiteCorCoefs(meth_site_cor_values = qsox1_cpg_correlations, reference_tss = qsox1_tss, 
  value_colours = "set2", xlabel = "Distance to TSS (bp)", ylabel = "DNA Methylation-Transcription Correlation") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "1: Calculate CpG methylation-transcription correlation values")

# Extract the colours for the points
plot_colours = ggplot_build(qsox1_cpg_correlation_plot)$data[[2]][["fill"]]

# Create qsox1 methodical scores plot
qsox1_methodical_scores_plot = plotMethodicalScores(meth_site_values = qsox1_cpg_correlations, reference_tss = qsox1_tss, xlabel = "Distance to TSS (bp)", 
  p_value_threshold = NULL, smooth_scores = F) + 
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "2: Convert correlation values to Methodical scores")  

# Add smoothed curve to qsox1 plot
qsox1_smoothed_methodical_scores_plot = plotMethodicalScores(meth_site_values = qsox1_cpg_correlations, 
  reference_tss = qsox1_tss, p_value_threshold = NULL, smooth_scores = T, smoothed_curve_colour = "hotpink2", curve_alpha = 1, xlabel = "Distance to TSS (bp)") +
geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(legend.position = c(0.9, 0.15)) + guides(fill = "none") +
  labs(title = "3: Smooth Methodical scores using exponential moving average") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) 

# Add TMRs and significance thresholds to plot
qsox1_tmr_plot = plotTMRs(meth_site_plot = qsox1_smoothed_methodical_scores_plot, tmrs_gr = qsox1_tmrs, reference_tss = qsox1_tss) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "4: Identify TMRs where smoothed Methodical curve crosses significance thresholds") +
  geom_hline(yintercept = log10(0.005), linetype = "dashed", colour = "#7B5C90") +
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", colour = "#BFAB25") +
  guides(fill = "none")
qsox1_tmr_plot

# Combine the plots 
qsox1_workflow_list = list(qsox1_cpg_correlation_plot, qsox1_methodical_scores_plot, qsox1_smoothed_methodical_scores_plot, qsox1_tmr_plot)
qsox1_workflow_list = lapply(qsox1_workflow_list, function(x) x + theme(plot.title = element_text(hjust = 0.5, size = 16)))
qsox1_workflow_plot = ggarrange(plotlist = qsox1_workflow_list, nrow = 2, ncol = 2, align = "hv")
qsox1_workflow_plot
ggsave(plot = qsox1_workflow_plot, "../figures/figure3.pdf", width = 20.57, height = 13.5, device = cairo_pdf)

### Create example TMR plots from metastasis samples

# Load CPGEA normal correlation results
mcrpc_correlation_results = readRDS("meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")

# Get GRanges for protein-coding TSS sites
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")

# Get GRanges for TMRs for normal prostate samples
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs.rds")

mcrpc_tmrs[mcrpc_tmrs$gene_name == "FOXD1"]
mcrpc_tmrs[mcrpc_tmrs$gene_name == "PACSIN3"]

# Create FOXD1 CpG correlation plot
foxd1_cpg_correlation_plot = plotMethSiteCorCoefs(meth_site_cor_values = mcrpc_correlation_results[["ENST00000615637"]], reference_tss = tss_gr["ENST00000615637"], 
  value_colours = "set2", xlabel = "Distance to *FOXD1* TSS (bp)", ylabel = "DNA Methylation-Transcription Correlation") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = ggtext::element_markdown(size = 18))

# Add TMRs to FOXD1 plot
foxd1_tmrs_plot = plotTMRs(meth_site_plot = foxd1_cpg_correlation_plot, tmrs_gr = mcrpc_tmrs[mcrpc_tmrs$ID == "ENST00000615637"], 
  reference_tss = tss_gr["ENST00000615637"]) +
  geom_hline(yintercept = 0, linetype = "dashed")

# Create pacsin3 CpG correlation plot
pacsin3_cpg_correlation_plot = plotMethSiteCorCoefs(meth_site_cor_values = mcrpc_correlation_results[["ENST00000298838"]], reference_tss = tss_gr["ENST00000298838"], 
  value_colours = "set2", xlabel = "Distance to *PACSIN3* TSS (bp)", ylabel = "DNA Methylation-Transcription Correlation") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = ggtext::element_markdown(size = 18))

# Add TMRs to pacsin3 plot
pacsin3_tmrs_plot = plotTMRs(meth_site_plot = pacsin3_cpg_correlation_plot, tmrs_gr = mcrpc_tmrs[mcrpc_tmrs$ID == "ENST00000298838"], 
  reference_tss = tss_gr["ENST00000298838"]) +
  geom_hline(yintercept = 0, linetype = "dashed")

# Load TMR correlation results
mcrpc_tmr_correlation_results = data.table::fread("tmr_evaluation_tables/mcrpc_tmrs_correlations_metastases_samples.tsv.gz")

# Get Correlation TMR correlation results for FOXD1
foxd1_tmr_cor_results = filter(mcrpc_tmr_correlation_results, transcript_name == "ENST00000615637")
foxd1_tmr_cor_results$significance = sig_sym(foxd1_tmr_cor_results$q_val)
foxd1_tmr_cor_barplot = ggplot(foxd1_tmr_cor_results, aes(x = genomic_region_name, y = cor, label = significance, fill = factor(sign(cor)))) +
  geom_col(color = "black") + 
  geom_text(nudge_y = ifelse(foxd1_tmr_cor_results$cor >= 0, 0.01, -0.025), size = 8) + theme_classic()
foxd1_tmr_cor_barplot = customize_ggplot_theme(foxd1_tmr_cor_barplot, ylab = "*FOXD1* TMR Correlations",
  x_labels = c("TMR1", "TMR2"), fill_colors = c("#A28CB1", "#D2C465"), show_legend = F) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  geom_hline(yintercept = 0) +
  theme(axis.title.y = ggtext::element_markdown(size = 18))
foxd1_tmr_cor_barplot

# Get Correlation TMR correlation results for pacsin3
pacsin3_tmr_cor_results = filter(mcrpc_tmr_correlation_results, transcript_name == "ENST00000298838")
pacsin3_tmr_cor_results$significance = sig_sym(pacsin3_tmr_cor_results$q_val)
pacsin3_tmr_cor_barplot = ggplot(pacsin3_tmr_cor_results, aes(x = genomic_region_name, y = cor, label = significance, fill = factor(sign(cor)))) +
  geom_col(color = "black") + 
  geom_text(nudge_y = ifelse(pacsin3_tmr_cor_results$cor >= 0, 0.01, -0.025), size = 8) + theme_classic()
pacsin3_tmr_cor_barplot = customize_ggplot_theme(pacsin3_tmr_cor_barplot, ylab = "*PACSIN3* TMR Correlations",
  x_labels = c("TMR1", "TMR2", "TMR3"), fill_colors = c("#A28CB1", "#D2C465"), show_legend = F) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  geom_hline(yintercept = 0) +
  theme(axis.title.y = ggtext::element_markdown(size = 18))
pacsin3_tmr_cor_barplot

# Combine plots 
correlation_plots = list(foxd1_tmrs_plot, foxd1_tmr_cor_barplot, pacsin3_tmrs_plot, pacsin3_tmr_cor_barplot)
combined_tmr_correlation_plots = plot_grid(plotlist = correlation_plots, nrow = 2, ncol = 2, align = "h", 
  rel_widths = c(3.5, 1), labels = c("A", "", "B", "")) 
saveRDS(combined_tmr_correlation_plots, "combined_tmr_correlation_plots.rds")
