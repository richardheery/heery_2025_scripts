# Load required packages
library(methodical)
library(dplyr)
library(ggplot2)
library(doParallel)
source("../auxillary_scripts/plotting_functions.R")

# Load GRanges for Gencode 36 MANE transcripts for hg19
gencode_mane_36_hg19 = readRDS("../auxillary_data/gencode_36_mane_transcripts_hg19_gr.rds")

# Get TSS for transcripts and expand
gencode_tss_hg19 = resize(gencode_mane_36_hg19, 1, fix = "start")
tss_flanking_regions = expand_granges(gencode_tss_hg19, upstream = 5000, downstream = 5000)

# Load meth RSE for TCGA hg19
tcga_meth_rse_hg19 = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/tcga_450k_array_hg19/")

# Get paths to DESeq2 normalized count tables for TCGA projects
deseq2_normalized_count_files = list.files("../auxillary_data/rnaseq_data/tcga_rna_seq/deseq_normalized_counts/", full.names = T)
names(deseq2_normalized_count_files) = gsub("_deseq_normalized_counts.tsv.gz", "", basename(deseq2_normalized_count_files))

# Calculate correlation for normal samples for each project
# Took 3 hours with 5 cores
bpparam = BiocParallel::MulticoreParam(workers = 10)


system.time({foreach(project = names(deseq2_normalized_count_files)) %dopar% {

  # Print name of current project
  message(paste("Starting", project))

  # Read in counts table for project
  counts_table = data.frame(data.table::fread(deseq2_normalized_count_files[project]), row.names = 1)

  # Find common samples to tcga_meth_rse_hg19 and counts table
  common_samples = intersect(names(counts_table), colnames(tcga_meth_rse_hg19))

  # Subset for normal samples
  common_normal_samples = common_samples[endsWith(common_samples, "11")]

  # If there are less than 10 normal samples, skip to next iteration
  if(length(common_normal_samples) < 10){
    message("Less than 10 common normal samples so skipping to next project")
    next
  }

  # Calculate correlation results for normal samples
  cor_results = methodical::calculateMethSiteTranscriptCors(
    meth_rse = tcga_meth_rse_hg19,
    transcript_expression_table = counts_table, samples_subset = common_normal_samples,
    tss_gr = gencode_tss_hg19, tss_associated_gr = tss_flanking_regions,
    cor_method = "s", BPPARAM = bpparam)

  # Save correlation results
  saveRDS(cor_results, paste0("tcga_probe_cors/", project, "_normal_sample_cors.rds"))

}})

# Load results for all normal samples
normal_sample_cor_results = list.files("tcga_probe_cors/", pattern = "normal", full.names = T)
names(normal_sample_cor_results) = gsub("_.*", "", basename(normal_sample_cor_results))
normal_sample_cor_results_list = lapply(normal_sample_cor_results, function(x) bind_rows(readRDS(x), .id = "transcript_id"))

# Correct p_values
normal_sample_cor_results_list  = lapply(normal_sample_cor_results_list, function(x) mutate(x, q_val = p.adjust(p_val, method = "fdr")))

# Combine all results into a single table
normal_sample_cor_results_combined = dplyr::bind_rows(normal_sample_cor_results_list, .id = "tissue")

# Bin CpGs by distance to TSS
normal_sample_cor_results_combined$bin = plyr::round_any(normal_sample_cor_results_combined$distance_to_tss, 500)

# Bin correlations and count number of significant correlations and all correlations in bin. 
normal_sample_cor_results_combined_binned = data.frame(summarize(group_by(normal_sample_cor_results_combined, bin), 
  num_sig = sum(q_val < 0.05, na.rm = T), total = sum(!is.na(p_val))))

# Calculate proportion of significant correlations
normal_sample_cor_results_combined_binned$prop_sig = normal_sample_cor_results_combined_binned$num_sig/normal_sample_cor_results_combined_binned$total

# Create barplot of significant correlations per bin and the percentage of significant correlations covered by probes
normal_prop_sig_bins_plot = ggplot(filter(normal_sample_cor_results_combined_binned, abs(bin) < 5000), aes(x = bin, y = prop_sig)) +
  geom_col(position = "identity", fill = "#2a5674") 
normal_prop_sig_bins_plot = customize_ggplot_theme(normal_prop_sig_bins_plot, xlab = "Distance to TSS (bp)", 
  ylab = "Proportion of CpG Sites\nDisplaying Significant Correlations", title = NULL)
normal_prop_sig_bins_plot = normal_prop_sig_bins_plot + scale_x_continuous(expand = c(0, 0), labels = scales::comma)
normal_prop_sig_bins_plot
saveRDS(normal_prop_sig_bins_plot, "all_tcga_normal_prop_sig_bins_plot.rds")

# Load results for all tumour samples
tumour_sample_cor_results = list.files("tcga_probe_cors/", pattern = "tumour", full.names = T)
names(tumour_sample_cor_results) = gsub("_.*", "", basename(tumour_sample_cor_results))
tumour_sample_cor_results_list = lapply(tumour_sample_cor_results, function(x) bind_rows(readRDS(x), .id = "transcript_id"))

# Correct p_values
tumour_sample_cor_results_list  = lapply(tumour_sample_cor_results_list, function(x) mutate(x, q_val = p.adjust(p_val, method = "fdr")))

# Combine all results into a single table
tumour_sample_cor_results_combined = dplyr::bind_rows(tumour_sample_cor_results_list, .id = "tissue")

# Bin CpGs by distance to TSS
tumour_sample_cor_results_combined$bin = plyr::round_any(tumour_sample_cor_results_combined$distance_to_tss, 500)

# Bin correlations and count number of significant correlations and all correlations in bin. 
tumour_sample_cor_results_combined_binned = data.frame(summarize(group_by(tumour_sample_cor_results_combined, bin), 
  num_sig = sum(q_val < 0.05, na.rm = T), total = sum(!is.na(p_val))))

# Calculate proportion of significant correlations
tumour_sample_cor_results_combined_binned$prop_sig = tumour_sample_cor_results_combined_binned$num_sig/tumour_sample_cor_results_combined_binned$total

# Create barplot of significant correlations per bin and the percentage of significant correlations covered by probes
tumour_prop_sig_bins_plot = ggplot(filter(tumour_sample_cor_results_combined_binned, abs(bin) < 5000), aes(x = bin, y = prop_sig)) +
  geom_col(position = "identity", fill = "#2a5674") 
tumour_prop_sig_bins_plot = customize_ggplot_theme(tumour_prop_sig_bins_plot, xlab = "Distance to TSS (bp)", 
  ylab = "Proportion of CpG Sites\nDisplaying Significant Correlations", title = NULL)
tumour_prop_sig_bins_plot = tumour_prop_sig_bins_plot + scale_x_continuous(expand = c(0, 0), labels = scales::comma)
tumour_prop_sig_bins_plot
saveRDS(tumour_prop_sig_bins_plot, "all_tcga_tumour_prop_sig_bins_plot.rds")

# Load plots necessary for completing figure 8
cpgs_per_bin_plot = readRDS("cpgs_per_bin_plot.rds")
combined_probes_per_bin_plot = readRDS("combined_probes_per_bin_plot.rds")

# Make combined plot with CpGs, probes and tumour plots
figure8_plot = cowplot::plot_grid(plotlist = 
    list(cpgs_per_bin_plot, combined_probes_per_bin_plot, normal_prop_sig_bins_plot, tumour_prop_sig_bins_plot), 
  nrow = 2, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"), byrow = T)
figure8_plot
ggsave(plot = figure8_plot, "../figures/figure8.pdf", width = 32, height = 18)

# Calculate correlation for tumour samples for each project
# Took 8 hours with 5 cores
bpparam = BiocParallel::MulticoreParam(workers = 10)
system.time({foreach(project = names(deseq2_normalized_count_files)) %do% {
  
  # Print name of current project
  message(paste("Starting", project))
  
  # Read in counts table for project
  counts_table = data.frame(data.table::fread(deseq2_normalized_count_files[project]), row.names = 1)
  
  # Find common samples to tcga_meth_rse_hg19 and counts table
  common_samples = intersect(names(counts_table), colnames(tcga_meth_rse_hg19))
  
  # Subset for tumour samples. 03 (Peripheral Blood) is used as tumour samples for LAML
  if(project == "LAML"){
    common_tumour_samples = common_samples[endsWith(common_samples, "03")]
  } else {
    common_tumour_samples = common_samples[endsWith(common_samples, "01")]
  }
  
  # If there are less than 10 tumour samples, skip to next iteration
  if(length(common_tumour_samples) < 10){
    message("Less than 10 common tumour samples so skipping to next project")
    next
  }
  
  # Calculate correlation results for tumour samples
  cor_results = methodical::calculateMethSiteTranscriptCors(
    meth_rse = tcga_meth_rse_hg19, 
    transcript_expression_table = counts_table, samples_subset = common_tumour_samples, 
    tss_gr = gencode_tss_hg19, tss_associated_gr = tss_flanking_regions, 
    cor_method = "s", BPPARAM = bpparam)
  
  # Save correlation results
  saveRDS(cor_results, paste0("tcga_probe_cors/", project, "_tumour_sample_cors.rds"))
  
}})

# Create a function which combines correlation results and creates summary plots
plot_cor_results = function(cor_results, title = NULL, breaks = 2500, results_dir){

  print(title)

  # Correct p-values
  cor_results = methodical::correct_correlation_pvalues(cor_results, p_adjust_method = "fdr")

  # Combine the results into a single table
  cor_results_combined = dplyr::bind_rows(cor_results, .id = "gene_id")

  # Add column with 500 bp bin which probe overlaps and filter for results in bins less than 5000
  cor_results_combined$bin = plyr::round_any(cor_results_combined$distance_to_tss, 500)
  cor_results_combined = filter(cor_results_combined, abs(bin) < 5000)

  # Add sign of correlation as a column
  cor_results_combined$sign = sign(cor_results_combined$cor)

  # Remove any correlation values of 0 or NA
  cor_results_combined = filter(cor_results_combined, !is.na(cor) & abs(cor) > 0)

  # Add the mean negative and positive correlations for each bin
  # and the proportion of significant results (defined as p-values < 0.005)
  cor_results_combined_summary = summarize(group_by(cor_results_combined, bin),
    mean_cor = mean(cor, na.rm = T), abs_cor = sum(abs(cor) > 0.3, na.rm = T)/sum(!is.na(cor)),
    prop_sig = sum(p_val < 0.05, na.rm = T)/sum(!is.na(p_val)))

  # Create a barplot with the proportion of significant results in each bin
  # prop_sig_plot = ggplot(cor_results_combined_summary, aes(x = bin, y = prop_sig, fill = factor(sign))) +
  #   geom_col(position = "dodge")
  # prop_sig_plot = customize_ggplot_theme(prop_sig_plot, xlab = "Distance to TSS (bp)",
  #   ylab = "Proportion of Significant Correlations", title = title,
  #   scale_x = scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(-breaks, breaks, breaks/2), labels = scales::comma),
  #   fill_colors = colour_list$purple_and_gold_light, fill_title = "Correlation\nDirection") +
  #   scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma)

  # Create a barplot with the proportion of significant results in each bin
  abs_cor_plot = ggplot(cor_results_combined_summary, aes(x = bin, y = abs_cor)) +
    geom_col(position = "identity", fill = "#d1eeea", color = "black")
  abs_cor_plot = customize_ggplot_theme(abs_cor_plot, xlab = "Distance to TSS (bp)",
    ylab = "Proportion of Illumina 450K Probes\nDisplaying Significant Correlations", title = title,
    scale_x = scale_x_continuous(expand = c(0, 0), labels = scales::comma),
    fill_colors = colour_list$purple_and_gold_light, fill_title = "Correlation\nDirection")

  # # Create a barplot with the mean correlation of each bin
  # mean_cor_plot = ggplot(cor_results_combined_summary, aes(x = bin, y = mean_cor, fill = factor(sign))) +
  #   geom_col()
  # mean_cor_plot = customize_ggplot_theme(mean_cor_plot, xlab = NULL, ylab = "Mean Correlation",
  #   scale_x = scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(-4000, 4000, 2000), labels = scales::comma),
  #   fill_colors = colour_list$purple_and_gold_light, fill_title = "Correlation\nDirection") +
  #   scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), labels = scales::comma)
  #
  # # Combine plots
  # combined_plot = ggpubr::ggarrange(plotlist = list(prop_sig_plot, mean_cor_plot), ncol = 1, common.legend = T, legend = "right")
  # combined_plot = ggpubr::annotate_figure(combined_plot, top = ggpubr::text_grob(title, size = 20),
  #   bottom = ggpubr::text_grob("Distance to TSS (bp)", size = 16))
  ggsave(plot = abs_cor_plot, paste0(results_dir, "/", title, ".pdf"), width = 16, height = 9)
  saveRDS(abs_cor_plot, paste0(results_dir, "/", title, ".rds"))
  return(NULL)

}

# Get paths to all nomral sample and all tumour sample results files
normal_results = list.files("tcga_probe_cors", pattern = "normal", full.names = T)
tumour_results = list.files("tcga_probe_cors", pattern = "tumour", full.names = T)

# Create correlations summary plots for all normal samples
system.time({
  normal_plots = lapply(normal_results, function(x)
    plot_cor_results(
      cor_results = readRDS(x),
      title = basename(gsub("_normal_sample_cors.rds", "", x)),
      results_dir = "normal_plots"
    ))
})

# Create correlations summary plots for all tumour samples
system.time({
  tumour_plots = lapply(tumour_results, function(x)
    plot_cor_results(
      cor_results = readRDS(x),
      title = basename(gsub("_tumour_sample_cors.rds", "", x)),
      results_dir = "tumour_plots"
    ))
})

selected_normal_plots = lapply(list.files("normal_plots/", full.names = T, pattern = ".rds")[c(2, 5, 8, 11)],
  readRDS)
selected_normal_plots =title_ggplots(
  selected_normal_plots,
  c(
    "Normal Breast Samples",
    "Normal Head and Neck Samples",
    "Normal Liver Samples",
    "Normal Endometrial Samples"
  )
)

selected_tumour_plots = lapply(list.files("tumour_plots/", full.names = T, pattern = ".rds")[c(2, 6, 16, 20)],
  readRDS)
selected_tumour_plots =title_ggplots(
  selected_tumour_plots,
  c(
    "Bladder Tumour Samples",
    "Colon Tumour Samples",
    "Lung Adenocarcinoma Tumour Samples",
    "Pancreatic Tumour Samples"
  )
)

combined_selected_plot_list = c(selected_normal_plots, selected_tumour_plots)
combined_selected_plot = cowplot::plot_grid(
  plotlist = combined_selected_plot_list,
  ncol = 2,
  nrow = 4,
  align = "hv",
  byrow = T
)
combined_selected_plot
ggsave(
  plot = combined_selected_plot,
  "../figures//supp_figure18.pdf",
  width = 32,
  height = 36
)




### 500 kb test
# Took 3 hours
system.time({
  for (project in names(deseq2_normalized_count_files)) {
    # Print name of current project
    message(paste("Starting", project))

    # Read in counts table for project
    counts_table = data.frame(data.table::fread(deseq2_normalized_count_files[project]),
      row.names = 1)

    # Find common samples to tcga_meth_rse_hg19 and counts table
    common_samples = intersect(names(counts_table), colnames(tcga_meth_rse_hg19))

    # Subset for normal samples
    common_normal_samples = common_samples[endsWith(common_samples, "11")]

    # If there are less than 10 normal samples, skip to next iteration
    if (length(common_normal_samples) < 10) {
      message("Less than 10 common normal samples so skipping to next project")
      next
    }

    # Calculate correlation results for normal samples
    cor_results = calculateMethSiteTranscriptCors(
      meth_rse = tcga_meth_rse_hg19,
      transcript_expression_table = counts_table,
      samples_subset = common_normal_samples,
      tss_gr = gencode_tss_hg19,
      expand_upstream = 500000,
      expand_downstream = 500000,
      cor_method = "s",
      BPPARAM = BiocParallel::MulticoreParam(workers = 5)
    )

    # Save correlation results
    saveRDS(
      cor_results,
      paste0(
        "tcga_probe_cors_500kb/",
        project,
        "_normal_sample_cors.rds"
      )
    )

  }
})

# Repeat for tumour. Took 5 hours

system.time({
  for (project in names(deseq2_normalized_count_files)) {
    # Print name of current project
    message(paste("Starting", project))

    # Read in counts table for project
    counts_table = data.frame(data.table::fread(deseq2_normalized_count_files[project]),
      row.names = 1)

    # Find common samples to tcga_meth_rse_hg19 and counts table
    common_samples = intersect(names(counts_table), colnames(tcga_meth_rse_hg19))

    # Subset for tumour samples
    common_tumour_samples = common_samples[endsWith(common_samples, "01")]

    # If there are less than 10 tumour samples, skip to next iteration
    if (length(common_tumour_samples) < 10) {
      message("Less than 10 common tumour samples so skipping to next project")
      next
    }

    # Calculate correlation results for tumour samples
    cor_results = calculateMethSiteTranscriptCors(
      meth_rse = tcga_meth_rse_hg19,
      transcript_expression_table = counts_table,
      samples_subset = common_tumour_samples,
      tss_gr = gencode_tss_hg19,
      expand_upstream = 500000,
      expand_downstream = 500000,
      cor_method = "s",
      BPPARAM = BiocParallel::MulticoreParam(workers = 20)
    )

    # Save correlation results
    saveRDS(
      cor_results,
      paste0(
        "tcga_probe_cors_500kb/",
        project,
        "_tumour_sample_cors.rds"
      )
    )

  }
})

# Get paths to all normal sample and all tumour sample results files
normal_results_500kb = list.files("tcga_probe_cors_500kb",
  pattern = "normal",
  full.names = T)
tumour_results_500kb = list.files("tcga_probe_cors_500kb",
  pattern = "tumour",
  full.names = T)

# Create correlations summary plots for all normal samples
system.time({
  normal_plots_500kb = lapply(normal_results_500kb, function(x)
    plot_cor_results(
      cor_results = readRDS(x),
      title = basename(gsub("_normal_sample_cors.rds", "", x)),
      breaks = 400000
    ))
})
pdf_save(plotlist = normal_plots_500kb,
  nrows = 1,
  filename = "tcga_normal_plots_500kb_abs.pdf")

# Create correlations summary plots for all tumour samples
system.time({
  tumour_plots_500kb = lapply(tumour_results_500kb, function(x)
    plot_cor_results(
      cor_results = readRDS(x),
      title = basename(gsub("_tumour_sample_cors.rds", "", x)),
      breaks = 400000
    ))
})
pdf_save(plotlist = tumour_plots_500kb,
  nrows = 1,
  filename = "tcga_tumour_plots_500kb_abs.pdf")


###

# Get paths to normal correlation result files
normal_cor_results_files = list.files("tcga_probe_cors", pattern = "normal", full.names = T)
names(normal_cor_results_files) = gsub("_normal_sample_cors.rds",
  "",
  basename(normal_cor_results_files))

blca_plots = plot_cor_results(readRDS(normal_cor_results_files["BLCA"]))
esca_plots = plot_cor_results(readRDS(normal_cor_results_files["ESCA"]))
coad_plots = plot_cor_results(readRDS(normal_cor_results_files["COAD"]))

###
# Took 22 minutes.
system.time({
  brca_cor_test_normal_samples = calculateMethSiteTranscriptCors(
    meth_rse = tcga_meth_rse_hg19,
    transcript_expression_table = brca_table,
    samples_subset = common_normal_samples,
    tss_gr = gencode_tss,
    expand_upstream = 5000,
    expand_downstream = 5000,
    cor_method = "s",
    BPPARAM = BiocParallel::MulticoreParam(workers = 20)
  )
})
saveRDS(brca_cor_test, "tcga_probe_cors/brca_cor_test.rds")
saveRDS(
  brca_cor_test_normal_samples,
  "tcga_probe_cors/brca_cor_test_normal_samples.rds"
)

brca_cor_test  = readRDS("tcga_probe_cors/brca_cor_test.rds")
brca_cor_test_normal_samples = readRDS("tcga_probe_cors/brca_cor_test_normal_samples.rds")


plot_cor_results(brca_cor_test, "BRCA Probe Correlation Summary")
plot_cor_results(brca_cor_test_normal_samples,
  "BRCA Probe Normal Samples Correlation Summary")


brca_cor_test$bin = plyr::round_any(brca_cor_test$distance_to_tss, 500)
brca_cor_test$sign = sign(brca_cor_test$cor)
brca_cor_test = filter(brca_cor_test, !is.na(cor) & abs(cor) > 0)
brca_cor_test_summary = summarize(
  group_by(brca_cor_test, bin, sign),
  mean_cor = mean(cor, na.rm = T),
  prop_sig = sum(p_val < 0.05, na.rm = T) / sum(!is.na(p_val))
)

ggplot(brca_cor_test_summary,
  aes(
    x = factor(bin),
    y = prop_sig,
    fill = factor(sign)
  )) +
  geom_col(position = "dodge")

ggplot(brca_cor_test_summary,
  aes(
    x = factor(bin),
    y = mean_cor,
    fill = factor(sign)
  )) +
  geom_col()
