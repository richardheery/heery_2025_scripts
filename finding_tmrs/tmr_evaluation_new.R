# Evaluate the TMR correlations

# Load required packages
library(methodical)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")

# Load CPGEA methylation RSE and MCRPC RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
mcrpc_meth_rse =  HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")

# Get transcript counts for CPGEA and MCRPC
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_transcript_counts.tsv.gz"), row.names = 1)
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_transcript_counts.tsv.gz"), row.names = 1)

# Get CPGEA normal and tumour samples and MCRPC samples with both WGBS and RNA-seq data
common_cpgea_normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_cpgea_tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Get TMR GRanges
# cpgea_normal_tmrs = rtracklayer::import.bed("tmr_bed_files/normal_prostate_tmrs.bed")
# cpgea_normal_tmrs$ID = gsub("_.*", "", cpgea_normal_tmrs$tmr_name)
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs_with_repeats.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs_with_repeats.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs_with_repeats.rds")

### Evaluate TMR correlations in normal prostate samples

#bpparam = BiocParallel::MulticoreParam(workers = 4)
bpparam = BiocParallel::SerialParam()

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BPPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_normal_samples, "tmr_evaluation_tables_repeats/cpgea_normal_tmrs_correlations_normal_samples.tsv.gz")

# # Filter for correlated and uncorrelated TMRs
# correlated_tmrs = filter(cpgea_normal_tmrs_correlations_normal_samples, q_val < 0.05)
# uncorrelated_tmrs = filter(cpgea_normal_tmrs_correlations_normal_samples, q_val > 0.05)
# boxplot(mean_tmr_coverage[uncorrelated_tmrs$genomic_region_name], mean_tmr_coverage[correlated_tmrs$genomic_region_name])
# 
# # Get the mean coverage for correlated and uncorrelated TMRs
# system.time({tmr_coverage = summarizeRegionMethylation(meth_rse = cpgea_meth_rse, assay = 2, genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, BPPARAM = bpparam)})
# system.time({tmr_meth = summarizeRegionMethylation(meth_rse = cpgea_meth_rse, assay = 1, genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, BPPARAM = bpparam)})
# mean_tmr_coverage = rowMeans(tmr_coverage, na.rm = T)
# boxplot(mean_tmr_coverage[uncorrelated_tmrs$genomic_region_name], mean_tmr_coverage[correlated_tmrs$genomic_region_name])

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_normal_samples, "tmr_evaluation_tables_repeats/cpgea_tumour_tmrs_correlations_normal_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_normal_samples, "tmr_evaluation_tables_repeats/mcrpc_tmrs_correlations_normal_samples.tsv.gz")

### Evaluate TMR correlations in prostate tumour samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_repeats/cpgea_normal_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_repeats/cpgea_tumour_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_tumour_samples, "tmr_evaluation_tables_repeats/mcrpc_tmrs_correlations_tumour_samples.tsv.gz")

### Evaluate TMR correlations in metastases samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_normal_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_repeats/cpgea_normal_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_repeats/cpgea_tumour_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", BBPARAM = bpparam)})
data.table::fwrite(mcrpc_tmrs_correlations_metastases_samples, "tmr_evaluation_tables_repeats/mcrpc_tmrs_correlations_metastases_samples.tsv.gz")

# Get paths to all tmr evaluation tables
tmr_evaluation_tables_paths = list.files("tmr_evaluation_tables/", full.names = T)
names(tmr_evaluation_tables_paths) = gsub(".tsv.gz", "", basename(tmr_evaluation_tables_paths))

# Read in all tables as a list
tmr_evaluation_tables = lapply(tmr_evaluation_tables_paths, function(x) data.table::fread(x))

# Add direction to TMRs
for(evaluation in names(tmr_evaluation_tables)){

  tmr_group = get(gsub("_correlations.*", "", evaluation))
  tmr_evaluation_tables[[evaluation]]$direction = tmr_group$direction[match(
    tmr_evaluation_tables[[evaluation]]$genomic_region_name, tmr_group$tmr_name)]

}

### Create plots of tmr correlations

# Get paths to all tmr evaluation tables
tmr_evaluation_tables_paths = list.files("tmr_evaluation_tables_repeats//", full.names = T)
names(tmr_evaluation_tables_paths) = gsub(".tsv.gz", "", basename(tmr_evaluation_tables_paths))

# Read in all tables as a list
tmr_evaluation_tables = lapply(tmr_evaluation_tables_paths, function(x) data.table::fread(x))

# Add direction to TMRs
for(evaluation in names(tmr_evaluation_tables)){

  tmr_group = get(gsub("_correlations.*", "", evaluation))
  tmr_evaluation_tables[[evaluation]]$direction = tmr_group$direction[match(
    tmr_evaluation_tables[[evaluation]]$genomic_region_name, tmr_group$tmr_name)]

}


# Combine all tables 
tmr_evaluation_tables = bind_rows(tmr_evaluation_tables, .id = "evaluation")

# Update q-values
tmr_evaluation_tables$q_val = p.adjust(tmr_evaluation_tables$p_val, method = "fdr")

# Add TMR group and sample group to table
tmr_evaluation_tables$tmr_group = gsub("_correlations.*", "", tmr_evaluation_tables$evaluation)
tmr_evaluation_tables$samples = gsub(".*_correlations_", "", tmr_evaluation_tables$evaluation)

# Add a column indicating which samples the TMRs are from so that external correlations can be easily identified
tmr_evaluation_tables$tmr_samples = case_when(
  tmr_evaluation_tables$tmr_group == "cpgea_normal_tmrs" ~ "normal_samples",
  tmr_evaluation_tables$tmr_group == "cpgea_tumour_tmrs" ~ "tumour_samples",
  tmr_evaluation_tables$tmr_group == "mcrpc_tmrs" ~ "metastases_samples")

# Convert samples and tmr_group to factor and put levels in correct order
tmr_evaluation_tables$samples = factor(tmr_evaluation_tables$samples, levels = c("normal_samples", "tumour_samples", "metastases_samples"))
tmr_evaluation_tables$tmr_group = factor(tmr_evaluation_tables$tmr_group, 
  levels = c("cpgea_normal_tmrs", "cpgea_tumour_tmrs", "mcrpc_tmrs"))

# Create violin plots showing the distribution TMR correlations
correlation_distribution_violins = ggplot(tmr_evaluation_tables, aes(x = samples, y = cor, fill = direction)) +
  geom_violin() + geom_hline(yintercept = 0, linetype = "dashed")
correlation_distribution_violins = customize_ggplot_theme(correlation_distribution_violins, facet = "tmr_group", facet_nrow = 3, fill_colors = colour_list$purple_and_gold_light, 
  x_labels = c("Normal Prostate Samples", "Prostate Tumour Samples", "Prostate Metastasis Samples"), 
  ylab = "Methylation-Transcription Correlation", fill_title = "TMR Direction",
  fill_labels = c("Negative", "Positive"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumour TMRs", "Prostate Metastasis TMRs")) +
  theme(strip.background = element_blank())
correlation_distribution_violins

# Count number of significant correlations by direction, tmr_group and sample group
significant_correlation_df = summarize(group_by(tmr_evaluation_tables, samples, tmr_group, direction), 
  proportion_sig = sum(q_val < 0.05, na.rm = T)/dplyr::n())

# Create barplots showing the number of significant correlations
significant_tmr_cor_barplots = ggplot(significant_correlation_df, aes(x = samples, y = proportion_sig, fill = direction)) +
  geom_col(position = "dodge", color = "black")
significant_tmr_cor_barplots = customize_ggplot_theme(significant_tmr_cor_barplots, facet = "tmr_group", facet_nrow = 3, fill_colors = colour_list$purple_and_gold_light, 
  x_labels = c("Normal Prostate Samples", "Prostate Tumour Samples", "Prostate Metastasis Samples"), 
  ylab = "Proportion of Correlations that are Significant", fill_title = "TMR Direction",
  fill_labels = c("Negative", "Positive"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumour TMRs", "Prostate Metastasis TMRs")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(strip.background = element_blank())
significant_tmr_cor_barplots

# Combine correlation_distribution_violins and significant_tmr_cor_barplots
combined_violins_and_barplot = correlation_distribution_violins / significant_tmr_cor_barplots

# Load other plots to make supplementary figure 10
tmrs_5kb_bins_plot = readRDS("tmrs_5kb_bins_plot.rds")
tmr_stats_barplot = readRDS("tmr_stats_barplot.rds")
tmr_stats_heatmaps = readRDS("tmr_stat_heatmaps.rds")

plotlist = list(tmrs_5kb_bins_plot, tmr_stats_barplot, tmr_stats_heatmaps, correlation_distribution_violins, significant_tmr_cor_barplots)
supp_figure_10 = ggarrange(plotlist = plotlist, nrow = 5, labels = LETTERS[1:5], widths = 27, heights = c(9, 9, 9, 11.815, 11.815))
ggsave(plot = supp_figure_10, "../figures/supp_figure10.pdf", width = 27, height = 50.63, limitsize = FALSE)

# Count proportion of significant correlations for external TMRs 
tmr_evaluation_tables_external = filter(tmr_evaluation_tables, samples != tmr_samples)
sum(tmr_evaluation_tables_external$q_val < 0.05, na.rm = T)/sum(!is.na(tmr_evaluation_tables_external$q_val))

# Count proportion of significant correlations for external TMRs involving prostate tumours and prostate metastases 
tmr_evaluation_tables_external_cancer = filter(tmr_evaluation_tables_external, samples != "normal_samples" & tmr_samples != "normal_samples")
sum(tmr_evaluation_tables_external_cancer$q_val < 0.05, na.rm = T)/sum(!is.na(tmr_evaluation_tables_external_cancer$q_val))

### Evaluate cpgea_normal_tmrs_1_or_more_cpg. Took 7 minutes
cpgea_normal_tmrs_1_or_more_cpg = readRDS("tmr_granges/cpgea_normal_tmrs_1_or_more_cpg_without_repeats.rds")
bpparam = BiocParallel::MulticoreParam(workers = 2)
system.time({cpgea_normal_tmrs_1_or_more_cpg_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs_1_or_more_cpg, genomic_region_names = cpgea_normal_tmrs_1_or_more_cpg$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_normal_tmrs_1_or_more_cpg$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", BPPARAM = bpparam)})

# Add direction and number of CpGs in TMRs
cpgea_normal_tmrs_1_or_more_cpg_correlations_normal_samples$direction = cpgea_normal_tmrs_1_or_more_cpg $direction
cpgea_normal_tmrs_1_or_more_cpg_correlations_normal_samples$cpg_count = cpgea_normal_tmrs_1_or_more_cpg $meth_site_count

# Get mean correlation for negative and positive TMRs by number of CpGs
mean_cor_by_number_of_cpgs = summarise(group_by(cpgea_normal_tmrs_1_or_more_cpg_correlations_normal_samples, direction, cpg_count), mean_cor = mean(cor, na.rm = T))

# Make a scatter plot of CpG count and mean correlation
cor_cpg_plot = ggplot(mean_cor_by_number_of_cpgs, aes(x = cpg_count, y = mean_cor, color = direction)) +
  geom_point()
cor_cpg_plot = customize_ggplot_theme(cor_cpg_plot, xlab = "Number of CpG Sites", ylab = "Mean TMR Correlation", 
  colors = c("#A28CB1", "#D2C465")) + theme(legend.position = "top")
ggsave(plot = cor_cpg_plot, "../figures/supp_figure2.pdf", width = 9, height = 9)