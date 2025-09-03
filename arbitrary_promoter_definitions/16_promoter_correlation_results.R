# Calculate correlations between promoter methylation and transcript expression for different promoter definitions.

# Load required packages
library(dplyr)
library(plotR)
library(ggpubr)
library(doParallel)
library(cowplot)
library(methodical)
source("../auxillary_scripts/correlate_functions.R")
source("../auxillary_scripts/mcrpc_example_differential_methylation_plot_functions.R")

# Get paths to all CPGEA promoter methylation definition tables
cpgea_promoter_definition_methylation_tables = readRDS("promoter_definition_methylation_tables/cpgea_promoter_definition_methylation_tables.rds")

# Get kallisto output for protein-coding genes subset for normal tumour samples 
cpgea_normalized_counts = data.frame(data.table::fread("../auxillary_data/cpgea_normalized_kallisto_pcg_counts.tsv.gz"), row.names = 1)
cpgea_normalized_counts_normal = dplyr::select(cpgea_normalized_counts, starts_with("N"))
cpgea_normalized_counts_tumour = dplyr::select(cpgea_normalized_counts, starts_with("T"))

# Create a directory for promoter methylation-transcription correlation results
dir.create("promoter_definition_transcript_correlation_tables")

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in normal samples. Took 50 minutes. 
system.time({for(definition in names(cpgea_promoter_definition_methylation_tables)){
  methylation_table = cpgea_promoter_definition_methylation_tables[[definition]]
  feature_matches_df = data.frame(cluster = row.names(methylation_table), row.names(methylation_table))
  correlation_results = cor_tables(table1 = methylation_table, table2 = cpgea_normalized_counts_normal, 
    feature_matches = feature_matches_df, calc_significance = T)
  data.table::fwrite(correlation_results, 
    paste0("promoter_definition_transcript_correlation_tables/", definition, "_definition_normal_sample_correlations.tsv.gz"), 
    sep = "\t", row.names = F, quote = F)
}})

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in tumour samples. Took 50 minutes. 
system.time({for(definition in names(cpgea_promoter_definition_methylation_tables)){
  methylation_table = cpgea_promoter_definition_methylation_tables[[definition]]
  feature_matches_df = data.frame(cluster = row.names(methylation_table), row.names(methylation_table))
  correlation_results = cor_tables(table1 = methylation_table, table2 = cpgea_normalized_counts_tumour, 
    feature_matches = feature_matches_df, calc_significance = T)
  data.table::fwrite(correlation_results, 
    paste0("promoter_definition_transcript_correlation_tables/", definition, "_definition_tumour_sample_correlations.tsv.gz"), 
    sep = "\t", row.names = F, quote = F)
}})

# Get paths to all MCRPC promoter methylation definition tables
mcrpc_promoter_definition_methylation_tables = readRDS("mcrpc_promoter_definition_methylation_tables/mcrpc_promoter_definition_methylation_tables.rds")

# Get kallisto output for protein-coding genes for MCRPC samples
mcrpc_normalized_counts = data.frame(data.table::fread("../auxillary_data/mcrpc_normalized_kallisto_pcg_counts.tsv.gz"), row.names = 1)

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in metastases samples. Took 50 minutes. 
system.time({for(definition in names(mcrpc_promoter_definition_methylation_tables)){
  methylation_table = mcrpc_promoter_definition_methylation_tables[[definition]]
  feature_matches_df = data.frame(cluster = row.names(methylation_table), row.names(methylation_table))
  correlation_results = cor_tables(table1 = methylation_table, table2 = mcrpc_normalized_counts, 
  feature_matches = feature_matches_df, calc_significance = T)
  data.table::fwrite(correlation_results, 
    paste0("promoter_definition_transcript_correlation_tables/", definition, "_definition_metastases_sample_correlations.tsv.gz"), 
    sep = "\t", row.names = F, quote = F)
}})

### Make plots for normal samples

# Get lists of all normal correlation tables
normal_sample_correlation_tables = list.files("promoter_definition_transcript_correlation_tables", pattern = "normal_sample_correlations", full.names = T)
names(normal_sample_correlation_tables) = LETTERS[1:5]

# Create a list with all normal sample correlation tables
normal_sample_correlation_list = lapply(normal_sample_correlation_tables, function(x) 
  setNames(data.table::fread(x, select = 1:5), c("table1_feature", "table2_feature", "cor", "p_value", "q_value")))

# Combine the list into a single table
normal_sample_correlation_tables_combined = bind_rows(normal_sample_correlation_list, .id = "definition")
normal_sample_correlation_tables_combined$definition = 
  factor(normal_sample_correlation_tables_combined$definition, levels = LETTERS[1:5])

# Remove correlations where q_value is NA
normal_sample_correlation_tables_combined = filter(normal_sample_correlation_tables_combined, !is.na(q_value))

# Convert q-value to a significance symbol
normal_sample_correlation_tables_combined$significance = plotR::sig_sym(normal_sample_correlation_tables_combined$q_value, symbol = "\u204E")

# Make violin plots for distributions of correlation values without clusters
normal_sample_all_correlations_violin_plots = 
  ggplot(normal_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28), axis.title = element_text(size = 24), 
    axis.text = element_text(size = 20), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Promoter Methylation vs\nTranscription Correlation Values", 
      title =  NULL)

# Denote whether correlations are positive, negative or uncorrelated
normal_sample_correlation_tables_combined = mutate(normal_sample_correlation_tables_combined, 
  correlation = case_when(
    q_value < 0.05 & cor > 0 ~ "Positive",
    q_value < 0.05 & cor < 0 ~ "Negative",
    q_value > 0.05 ~ "Uncorrelated"
    )
  )

# Convert correlation to a factor
normal_sample_correlation_tables_combined$correlation = factor(normal_sample_correlation_tables_combined$correlation, levels = c("Negative", "Uncorrelated", "Positive"))

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
normal_sample_correlation_tables_combined_summary = mutate(
  summarize(group_by(normal_sample_correlation_tables_combined, definition, correlation), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated unchanged promoters for each definition
normal_correlation_proportions_barplot = ggplot(filter(normal_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")

# Adjust theme of barplot and save
normal_correlation_proportions_barplot = customize_ggplot_theme(normal_correlation_proportions_barplot, 
  title = NULL, 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = c("#7B5C90B4", "#bfab25B4"), 
  plot_title_size = 28, axis_title_size = 24, axis_text_size = 20, legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  scale_y_continuous(limits = c(0, 0.1), exp = c(0, 0)) + theme(legend.position = c(0.85, 0.9))

# Combine normal plots and save
combined_normal_correlation_plots = ggarrange(plotlist = list(normal_sample_all_correlations_violin_plots, normal_correlation_proportions_barplot),
  labels = c("D", "E"))

### Repeat for tumour samples

# Get lists of all tumour correlation tables
tumour_sample_correlation_tables = list.files("promoter_definition_transcript_correlation_tables", pattern = "tumour_sample_correlations", full.names = T)
names(tumour_sample_correlation_tables) = LETTERS[1:5]

# Create a list with all tumour sample correlation tables
tumour_sample_correlation_list = lapply(tumour_sample_correlation_tables, function(x) 
  setNames(data.table::fread(x, select = 1:5), c("table1_feature", "table2_feature", "cor", "p_value", "q_value")))

# Combine the list into a single table
tumour_sample_correlation_tables_combined = bind_rows(tumour_sample_correlation_list, .id = "definition")
tumour_sample_correlation_tables_combined$definition = 
  factor(tumour_sample_correlation_tables_combined$definition, levels = LETTERS[1:5])

# Remove correlations where q_value is NA
tumour_sample_correlation_tables_combined = filter(tumour_sample_correlation_tables_combined, !is.na(q_value))

# Make violin plots for distributions of correlation values without clusters
tumour_sample_all_correlations_violin_plots = 
  ggplot(tumour_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Promoter Methylation vs\nTranscription Correlation Values", 
      title = NULL)

# Denote whether correlations are positive, negative or uncorrelated
tumour_sample_correlation_tables_combined = mutate(tumour_sample_correlation_tables_combined, 
  correlation = case_when(
    q_value < 0.05 & cor > 0 ~ "Positive",
    q_value < 0.05 & cor < 0 ~ "Negative",
    q_value > 0.05 ~ "Uncorrelated"
    )
  )

# Convert correlation to a factor
tumour_sample_correlation_tables_combined$correlation = factor(tumour_sample_correlation_tables_combined$correlation, levels = c("Negative", "Uncorrelated", "Positive"))

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
tumour_sample_correlation_tables_combined_summary = mutate(
  summarize(group_by(tumour_sample_correlation_tables_combined, definition, correlation), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated unchanged promoters for each definition
tumour_correlation_proportions_barplot = ggplot(filter(tumour_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")

# Adjust theme of barplot save
tumour_correlation_proportions_barplot = customize_ggplot_theme(tumour_correlation_proportions_barplot, title = NULL, 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = colour_list$purple_and_gold_light) +
  theme(legend.position = c(0.9, 0.9))

### Repeat for metastases samples

# Get lists of all metastases correlation tables
metastases_sample_correlation_tables = list.files("promoter_definition_transcript_correlation_tables", pattern = "metastases_sample_correlations", full.names = T)
names(metastases_sample_correlation_tables) = LETTERS[1:5]

# Create a list with all metastases sample correlation tables
metastases_sample_correlation_list = lapply(metastases_sample_correlation_tables, function(x) 
  setNames(data.table::fread(x, select = 1:5), c("table1_feature", "table2_feature", "cor", "p_value", "q_value")))

# Combine the list into a single table
metastases_sample_correlation_tables_combined = bind_rows(metastases_sample_correlation_list, .id = "definition")
metastases_sample_correlation_tables_combined$definition = 
  factor(metastases_sample_correlation_tables_combined$definition, levels = LETTERS[1:5])

# Remove correlations where q_value is NA
metastases_sample_correlation_tables_combined = filter(metastases_sample_correlation_tables_combined, !is.na(q_value))

# Make violin plots for distributions of correlation values without clusters
metastases_sample_all_correlations_violin_plots = 
  ggplot(metastases_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Promoter Methylation vs\nTranscription Correlation Values", 
      title = NULL)

# Denote whether correlations are positive, negative or uncorrelated
metastases_sample_correlation_tables_combined = mutate(metastases_sample_correlation_tables_combined, 
  correlation = case_when(
    q_value < 0.05 & cor > 0 ~ "Positive",
    q_value < 0.05 & cor < 0 ~ "Negative",
    q_value > 0.05 ~ "Uncorrelated"
    )
  )

# Convert correlation to a factor
metastases_sample_correlation_tables_combined$correlation = factor(metastases_sample_correlation_tables_combined$correlation, levels = c("Negative", "Uncorrelated", "Positive"))

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
metastases_sample_correlation_tables_combined_summary = mutate(
  summarize(group_by(metastases_sample_correlation_tables_combined, definition, correlation), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated unchanged promoters for each definition
metastases_correlation_proportions_barplot = ggplot(filter(metastases_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")

# Adjust theme of barplot save
metastases_correlation_proportions_barplot = customize_ggplot_theme(metastases_correlation_proportions_barplot, title = NULL, 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = colour_list$purple_and_gold_light) +
  theme(legend.position = c(0.9, 0.9))

# Combine metastases plots and save
combined_tumour_and_metastases_plots = ggarrange(plotlist = list(tumour_sample_all_correlations_violin_plots, metastases_sample_all_correlations_violin_plots, 
  tumour_correlation_proportions_barplot, metastases_correlation_proportions_barplot),
  labels = c("A", "B", "C", "D"))
combined_tumour_and_metastases_plots
ggsave(plot = combined_tumour_and_metastases_plots, filename = "../figures/supp_figure8.pdf", width = 20.57, height = 24.34)

### Make example plots for MCRPC samples

# Add significance symbol
metastases_sample_correlation_tables_combined$significance = sig_sym(metastases_sample_correlation_tables_combined$q_value)

# Load MCRPC correlation results
mcrpc_correlation_results = readRDS("../finding_tmrs/meth_transcript_cors/mcrpc_whole_gene_body_correlations.rds")

# Load location of TSS sites
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")

# Create plots of the promoter correlation values for FOXD1 and PACSIN3 using the ENST00000615637 and ENST00000298838 transcripts
foxd1_promoter_cor_plot = plot_promoter_correlations(transcript = "ENST00000615637", title = "*FOXD1* Promoter Methylation-<br>Transcription Correlation")
pacsin3_promoter_cor_plot = plot_promoter_correlations(transcript = "ENST00000298838", title = "*PACSIN3* Promoter Methylation-<br>Transcription Correlation")

# Create plots of the individual CpG correlation values for FOXD1 and PACSIN3 using the ENST00000615637 and ENST00000298838 transcripts
foxd1_cpg_cor_plot = plot_cpg_cor_values_transcript(transcript = "ENST00000615637", xlabel = "Distance to *FOXD1* TSS (bp)")
pacsin3_cpg_cor_plot = plot_cpg_cor_values_transcript(transcript = "ENST00000298838", xlabel = "Distance to *PACSIN3* TSS (bp)")

# Load promoter coordinates plot
promoter_region_plot = readRDS("promoter_region_plot.rds")

# Combine plots into a single figure along with promoter definition plot and save
promoters_and_examples_plot_list = list(promoter_region_plot, NULL, 
  foxd1_cpg_cor_plot, foxd1_promoter_cor_plot,
  pacsin3_cpg_cor_plot, pacsin3_promoter_cor_plot)
promoters_and_examples_plot = plot_grid(plotlist = promoters_and_examples_plot_list, nrow = 3, ncol = 2, align = "hv", 
  rel_heights = c(1.5, 6, 6), rel_widths = c(3.5, 1), labels = c("A", "", "B", "", "C", ""))

# Add barplot to complete plot
figure2_plot_list = list(promoters_and_examples_plot, combined_normal_correlation_plots)
figure2_plot = plot_grid(plotlist = figure2_plot_list, nrow = 2, ncol = 1, rel_heights = c(13.5, 5.5))
figure2_plot
ggsave(plot = figure2_plot, filename = "../figures/figure2.pdf", 
  width = 20.57, height = 30.07, device = cairo_pdf)
