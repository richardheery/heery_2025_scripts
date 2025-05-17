### Compare methylation change at TMRs to other genomic features in prostate tumours and TCGA samples

# Load required packages
library(methodical)
library(dplyr)
library(doParallel)
source("../auxillary_scripts/plotting_functions.R")

# Get methylation values from CPGEA 
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")

# Get genome annotation for hg38
genome_annotation = readRDS("../auxillary_data/complete_regulatory_annotation.rds")

# Get repeat ranges for hg38
repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

# Remove any features which overlap repeats
genome_annotation = subsetByOverlaps(genome_annotation, repeat_ranges, invert = T, ignore.strand = T)

# Get a list of TMRs
tmr_list = readRDS("../finding_tmrs/tmr_granges/tmr_list.rds")
tmr_list = lapply(tmr_list, unname)

# Rename TMRs
names(tmr_list) = c("Normal Prostate TMRs -", "Normal Prostate TMRs +", "Prostate Tumour TMRs -", 
  "Prostate Tumour TMRs +", "Prostate Metastasis TMRs -", "Prostate Metastasis TMRs +")

# Combine into a single GRanges
tmrs_gr = unlist(GRangesList(tmr_list))

# Remove metadata columns
mcols(tmrs_gr) = NULL

# Create a column with TMR group
tmrs_gr$region_type = names(tmrs_gr)
names(tmrs_gr) = NULL

# Add TMRs to genome annotation
genome_annotation_with_tmrs = c(genome_annotation, tmrs_gr)
gc()

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(5)

# Get methylation values for genomic features and TMRs. Took 3 minutes with 10 cores. 
system.time({genomic_features_with_tmrs_methylation = 
  summarizeRegionMethylation(meth_rse = cpgea_rse, genomic_regions = genome_annotation_with_tmrs, BPPARAM = bpparam)})

# Save table
data.table::fwrite(genomic_features_with_tmrs_methylation, "genomic_features_with_tmrs_methylation.tsv.gz")

# Add region type from genomic_features_with_tmrs
genomic_features_with_tmrs_methylation$region_type = genome_annotation_with_tmrs$region_type
gc()

# Split by region_type
genomic_features_with_tmrs_methylation = split(
  select(genomic_features_with_tmrs_methylation, -region_type), genomic_features_with_tmrs_methylation$region_type)
gc()

# Get the mean methylation change for each region
genomic_feature_mean_methylation_change = lapply(genomic_features_with_tmrs_methylation, function(x)
  rowMeans(select(x, starts_with("T")) - select(x, starts_with("N")), na.rm = T))

# Convert genomic_feature_mean_methylation_change into a single table and save
genomic_feature_mean_methylation_change = data.frame(
  region_type = rep(names(genomic_feature_mean_methylation_change), times = lengths(genomic_feature_mean_methylation_change)),
  values = unlist(genomic_feature_mean_methylation_change), row.names = NULL)
data.table::fwrite(genomic_feature_mean_methylation_change, "genomic_feature_mean_methylation_change.tsv.gz")

# Create a vector with the selected genomic features to plot
selected_genomic_features = c("CpG Islands", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", 
  "CTCF BS", "Exon", "Intron")

# Filter genomic_feature_mean_methylation_change for selected features
genomic_feature_mean_methylation_change_selected = filter(genomic_feature_mean_methylation_change, 
  region_type %in% c(selected_genomic_features, c("Prostate Tumour TMRs -", "Prostate Tumour TMRs +")))

# Put region_type in desired order
genomic_feature_mean_methylation_change_selected$region_type = 
  factor(genomic_feature_mean_methylation_change_selected$region_type, 
    levels = rev(c("Prostate Tumour TMRs -", "Prostate Tumour TMRs +", selected_genomic_features)))

# Set colours for genomic features
feature_colors = rev(c(rep(colour_list$purple_and_gold_light, each = 1), 
  rev(RColorBrewer::brewer.pal(7, name = "BrBG"))))

# Randomly select 250 regions for each region type for plotting points
set.seed(123)
sample_regions = slice_sample(group_by(genomic_feature_mean_methylation_change_selected, region_type), n = 250)

# Create boxplots of mean methylation change for genomic features and TMRs
set.seed(123)
methylation_change_boxplots = ggplot(genomic_feature_mean_methylation_change_selected, 
  aes(x = region_type, y = values, fill = region_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = sample_regions, alpha = 1, aes(fill = region_type), shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip()
methylation_change_boxplots = customize_ggplot_theme(methylation_change_boxplots, 
  ylab = "Mean Methylation Change", xlab = "TMRs and Regulatory Features", show_legend = F, 
  scale_y = scale_y_continuous(limits = c(-0.5, 0.5)), 
  fill_colors = feature_colors)
methylation_change_boxplots
methylation_change_boxplots = ggpubr::ggarrange(methylation_change_boxplots, labels = "B")
ggsave(plot = methylation_change_boxplots, "../figures/figure6B.pdf", width = 8, height = 9)

### Check methylation change at regions in 8 TCGA patients with WGBS data for matching tumour and normal samples 

# Load meth RSE for TCGA WGBS data
tcga_wgbs_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/tcga_wgbs_hg38/")

# Split genome_annotation_with_tmrs into a list
genome_annotation_with_tmrs_list = split(genome_annotation_with_tmrs, genome_annotation_with_tmrs$region_type)

# Remove TMRs from normal prostate and prostate tumour samples and rename metastasis TMRs
genome_annotation_with_tmrs_list = genome_annotation_with_tmrs_list[grep("Normal|Tumour", names(genome_annotation_with_tmrs_list), invert = T, value = T)]
names(genome_annotation_with_tmrs_list)[c(11, 12)] = c("Positive TMRs", "Negative TMRs")

# Get the mean methylation of CpGs for each genomic region. Took 10 minutes with 1 core locally
system.time({tcga_wgbs_genomic_region_mean_meth = foreach(region = names(genome_annotation_with_tmrs_list), .packages = "methodical") %do% {
  
  message(paste("Starting", region))
  
  region_cpgs = subsetByOverlaps(tcga_wgbs_meth_rse, genome_annotation_with_tmrs_list[[region]])
  system.time({region_methylation = colMeans(assay(region_cpgs), na.rm = T)})
  region_methylation
  
}})

# Add names to results and combine into a single table
names(tcga_wgbs_genomic_region_mean_meth) = names(genome_annotation_with_tmrs_list)
tcga_wgbs_genomic_region_mean_meth = dplyr::bind_rows(tcga_wgbs_genomic_region_mean_meth, .id = "region_type")
data.table::fwrite(tcga_wgbs_genomic_region_mean_meth, "tcga_wgbs_genomic_region_mean_meth_mcrpc_tmrs.tsv.gz")
tcga_wgbs_genomic_region_mean_meth = data.table::fread("tcga_wgbs_genomic_region_mean_meth_mcrpc_tmrs.tsv.gz")
tcga_wgbs_genomic_region_mean_meth = tibble::column_to_rownames(tcga_wgbs_genomic_region_mean_meth, "region_type")

# Identify submitters with matching tumour and normal samples. There is one submitter each for BLCA, BRCA, COAD, LUAD, LUSC, READ, STAD and UCEC
normal_samples = rownames(colData(tcga_wgbs_meth_rse))[which(colData(tcga_wgbs_meth_rse)$sample_type == 11)]
matching_tumour_samples = gsub("_11", "_01", normal_samples)

# Find mean methylation change for all regions 
tcga_wgbs_genomic_region_mean_meth_tumour = select(tcga_wgbs_genomic_region_mean_meth, all_of(matching_tumour_samples))
tcga_wgbs_genomic_region_mean_meth_normal = select(tcga_wgbs_genomic_region_mean_meth, all_of(normal_samples))
tcga_wgbs_genomic_region_mean_meth_change = tcga_wgbs_genomic_region_mean_meth_tumour - tcga_wgbs_genomic_region_mean_meth_normal

# Change row names to a column and convert to long format
tcga_wgbs_genomic_region_mean_meth_change = tibble::rownames_to_column(tcga_wgbs_genomic_region_mean_meth_change, "region_type")
tcga_wgbs_genomic_region_mean_meth_change = tidyr::pivot_longer(tcga_wgbs_genomic_region_mean_meth_change, 
  cols = -region_type, names_to = "sample_name", values_to = "mean_meth_change")

# Add tumour type
tcga_wgbs_genomic_region_mean_meth_change$cancer = colData(tcga_wgbs_meth_rse)[gsub("_1", "_01", tcga_wgbs_genomic_region_mean_meth_change$sample_name), ]$project

# Give TMRs a diamond shape
tcga_wgbs_genomic_region_mean_meth_change$shape = 
  ifelse(grepl("TMR", tcga_wgbs_genomic_region_mean_meth_change$region_type), 23, 21)

# Filter tcga_wgbs_genomic_region_mean_meth_change  for selected regions
tcga_wgbs_genomic_region_mean_meth_change = filter(tcga_wgbs_genomic_region_mean_meth_change, 
  region_type %in% c("Negative TMRs", "Positive TMRs", selected_genomic_features))

# Convert region_type to a factor with levels in correct order
genomic_feature_plot_order = c("Negative TMRs", "Positive TMRs", selected_genomic_features)
tcga_wgbs_genomic_region_mean_meth_change$region_type = 
  factor(tcga_wgbs_genomic_region_mean_meth_change$region_type, levels = genomic_feature_plot_order)

# Reverse order of cancer type so they are in alphabetical order from top to bottom
tcga_wgbs_genomic_region_mean_meth_change$cancer = 
  factor(tcga_wgbs_genomic_region_mean_meth_change$cancer, levels = rev(unique(tcga_wgbs_genomic_region_mean_meth_change$cancer)))
levels(tcga_wgbs_genomic_region_mean_meth_change$cancer) = 
  c("Endometrial", "Stomach", "Rectal", "Lung Squamous", "Lung Adenocarcinoma", "Colon", "Breast", "Bladder")

# Create a plot of the mean methylation change of probes for different genomic features in TCGA
tcga_wgbs_feature_meth_change_plot = ggplot(tcga_wgbs_genomic_region_mean_meth_change,
  aes(y = cancer, x = mean_meth_change, fill = region_type, color = region_type)) +
  geom_point(shape = tcga_wgbs_genomic_region_mean_meth_change$shape,
   size = 6, position = position_jitter(seed = 2, height = 0.3), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed")  +
  geom_hline(yintercept = seq(0.5, length(levels(tcga_wgbs_genomic_region_mean_meth_change$cancer)), by = 1),
    color = "gray", linewidth = .5, alpha = .5)
tcga_wgbs_feature_meth_change_plot = 
  customize_ggplot_theme(tcga_wgbs_feature_meth_change_plot, 
  ylab = "Cancer Type", xlab = "Mean Methylation Change", color_title = "Genomic Feature",
  show_legend = T, scale_x = scale_x_continuous(limits = c(-0.15, 0.3), expand = c(0, 0)), 
  fill_colors = rev(feature_colors), colors = rev(feature_colors), legend_text_size = 16) + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
    panel.grid.major.y = element_blank()) +
    guides(fill = "none")
tcga_wgbs_feature_meth_change_plot
ggsave(plot = tcga_wgbs_feature_meth_change_plot, filename = "../figures/figure6C.pdf", height = 9, width = 16)