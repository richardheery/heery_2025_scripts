# Create plot of Methodical workflow and show examples of some TMRs

# Load required packages
library(methodical)
library(ggpubr)

# Load CPGEA normal correlation results
cpgea_normal_correlation_results = readRDS("meth_transcript_cors/cpgea_normal_whole_gene_body_correlations.rds")

# Get GRanges for protein-coding TSS sites
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")

# Gt GRanges for TMRs for normal prostate samples
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs_with_repeats.rds")

# Show TMR workflow using LILRB1 (ENST00000324602) as an example

# Get CpG correlation values for LILRB1 and remove missing values
LILRB1_cpg_correlations = cpgea_normal_correlation_results[["ENST00000324602"]]
LILRB1_cpg_correlations = LILRB1_cpg_correlations[complete.cases(LILRB1_cpg_correlations), ]
row.names(LILRB1_cpg_correlations) = NULL

# Get LILRB1 TSS site
LILRB1_tss = tss_gr[tss_gr$ID == "ENST00000324602"]

# Get TMRs for LILRB1
LILRB1_tmrs = cpgea_normal_tmrs[cpgea_normal_tmrs$ID == "ENST00000324602"]

# Create LILRB1 CpG correlation plot
LILRB1_cpg_correlation_plot = plotMethSiteCorCoefs(meth_site_cor_values = LILRB1_cpg_correlations, reference_tss = LILRB1_tss, 
  value_colours = "set2", xlabel = "Distance to TSS (bp)", ylabel = "DNA Methylation-Transcription Correlation") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "1: Calculate CpG methylation-transcription correlation values")

# Extract the colours for the points
plot_colours = ggplot_build(LILRB1_cpg_correlation_plot)$data[[2]][["fill"]]

# Create LILRB1 methodical scores plot
LILRB1_methodical_scores_plot = plotMethodicalScores(meth_site_values = LILRB1_cpg_correlations, reference_tss = LILRB1_tss, xlabel = "Distance to TSS (bp)", 
  p_value_threshold = NULL, smooth_scores = F) + 
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "2: Convert correlation values to Methodical scores")  

# Add smoothed curve to LILRB1 plot
LILRB1_smoothed_methodical_scores_plot = plotMethodicalScores(meth_site_values = LILRB1_cpg_correlations, 
  reference_tss = LILRB1_tss, p_value_threshold = NULL, smooth_scores = T, smoothed_curve_colour = "hotpink2", curve_alpha = 1, xlabel = "Distance to TSS (bp)") +
geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(legend.position = c(0.9, 0.15)) + guides(fill = "none") +
  labs(title = "3: Smooth Methodical scores using exponential moving average") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) 

# Add TMRs and significance thresholds to plot
LILRB1_tmr_plot = plotTMRs(meth_site_plot = LILRB1_smoothed_methodical_scores_plot, tmrs_gr = LILRB1_tmrs, reference_tss = LILRB1_tss) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "4: Identify TMRs where smoothed Methodical curve crosses significance thresholds") +
  geom_hline(yintercept = log10(0.005), linetype = "dashed", colour = "#7B5C90") +
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", colour = "#BFAB25") +
  guides(fill = "none")
LILRB1_tmr_plot

# Combine the plots 
LILRB1_workflow_list = list(LILRB1_cpg_correlation_plot, LILRB1_methodical_scores_plot, LILRB1_smoothed_methodical_scores_plot, LILRB1_tmr_plot)
LILRB1_workflow_list = lapply(LILRB1_workflow_list, function(x) x + theme(plot.title = element_text(hjust = 0.5, size = 16)))
LILRB1_workflow_plot = ggarrange(plotlist = LILRB1_workflow_list, nrow = 2, ncol = 2, align = "hv")
LILRB1_workflow_plot
ggsave(plot = LILRB1_workflow_plot, "../figures/figure4.pdf", width = 20.57, height = 13.5)