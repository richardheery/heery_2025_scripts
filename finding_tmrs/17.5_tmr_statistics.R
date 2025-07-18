# Create plots of TMR numbers and TMR overlaps

# Load required packages
source("../auxillary_scripts/granges_functions.R")
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/enrichment_tests.R")

# Load TMRs
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs.rds")

# Create a list with the different TMR types
tmr_list = list(
  cpgea_normal = cpgea_normal_tmrs,
  cpgea_tumour = cpgea_tumour_tmrs,
  mcrpc = mcrpc_tmrs
)

# Get transcripts and genes associated with each dataset
tmr_transcripts = lapply(tmr_list, function(x) unique(x$ID))
tmr_genes = lapply(tmr_list, function(x) unique(x$gene_name))
lengths(tmr_transcripts)
lengths(tmr_genes)

# Get list of TMRs separated by direction
tmr_list = readRDS("tmr_granges/tmr_list.rds")

# Create a data.frame summarizing the number of TMRs, transcripts and genes associated with each TMR group
tmr_stats = data.frame(
  dataset = gsub("_negative|_positive", "", names(tmr_list)),
  direction = stringr::str_to_title(gsub(".*_", "", names(tmr_list))),
  tmr_count = lengths(tmr_list),
  transcript_count = sapply(tmr_list, function(x) length(unique(x$ID))),
  gene_count = sapply(tmr_list, function(x) length(unique(x$gene_name))),
  row.names = NULL
)

# Put dataset levels in correct order
tmr_stats$dataset = factor(tmr_stats$dataset, unique(tmr_stats$dataset))

# Convert tmr_stats into long format
tmr_stats = tidyr::pivot_longer(tmr_stats, cols = c("tmr_count", "transcript_count", "gene_count"))
tmr_stats$name = factor(tmr_stats$name, levels = c("tmr_count", "transcript_count", "gene_count"))

# Create barplots with the number of TMRs, TMR-associated transcripts and TMR-associated genes for each dataset
tmr_stats_barplot = ggplot(tmr_stats, aes(x = dataset, y = value, fill = direction)) +
  geom_col(position = "dodge", color = "black")
tmr_stats_barplot = customize_ggplot_theme(tmr_stats_barplot, 
  title = NULL, xlab = "Dataset", y = "Count", fill_title = "TMR Direction",
  fill_colors = colour_list$purple_and_gold_light, x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), 
  facet = "name", facet_labels = c("Number of TMRs", "Number of TMR-Associated Transcripts", "Number of TMR-Associated Genes"), 
  facet_scales = "fixed", x_labels_angle = 30) + theme(strip.background = element_blank())
tmr_stats_barplot
saveRDS(tmr_stats_barplot, "tmr_stats_barplot.rds")

# Calculate the proportion overlap between TMR groups
tmr_overlaps = calculate_regions_overlap_list(grl = tmr_list, ignore.strand = T, overlap_threshold = 0.25)

# Get overlap proportions for transcripts and genes from different TMR groups
transcript_overlaps = intersect_lengths_all_pairwise(lapply(tmr_list, function(x) x$ID), proportion = T)
gene_overlaps = intersect_lengths_all_pairwise(lapply(tmr_list, function(x) x$gene_name), proportion = T)

# Create row and column labels for heatmaps 
labels = c("Normal Prostate -", "Normal Prostate +", "Prostate Tumour -", 
  "Prostate Tumour +", "Prostate Metastasis -", "Prostate Metastasis +")

# Create heatmaps of TMR overlaps, for Jaccard indices for transcripts and genes
tmr_overlaps_plot = heatmap_without_clustering(tmr_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Genomic Regions\nCovered by TMRs", title_size = 15, return_ggplot = T)
transcript_overlaps_plot = heatmap_without_clustering(transcript_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Transcripts\nAssociated with TMRs", title_size = 15, return_ggplot = T)
gene_overlaps_plot = heatmap_without_clustering(gene_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Genes\nAssociated with TMRs", title_size = 15, return_ggplot = T)

# Combine heatmaps and save
tmr_heatmap_list = list(tmr_overlaps_plot, transcript_overlaps_plot, gene_overlaps_plot)
combined_tmr_heatmaps = ggpubr::ggarrange(plotlist = tmr_heatmap_list, nrow = 1, ncol = 3)
combined_tmr_heatmaps
saveRDS(combined_tmr_heatmaps, "tmr_stat_heatmaps.rds")

# Test the enrichment of hallmark pathways in genes only associated with prostate tumour and prostgate metastases TMRs

# Get background genes  
background_genes = unique(readRDS("../auxillary_data/cage_supported_gencode_tss.rds")$gene_name)

# Load MSigDB gene sets
msigdb_gene_set_list = readRDS("../auxillary_data/msigdb_complete_gene_set_list.rds")

# Get genes associated with TMRs both in prostate tumours and prostate metastases but not in normal samples
cancer_tmr_genes = setdiff(intersect(cpgea_tumour_tmrs$gene_name, mcrpc_tmrs$gene_name), cpgea_normal_tmrs$gene_name)

# Test the enrichment of hallmark pathways and filter for significant pathways
cancer_tmr_genes_hallmark_pathways_enrichment = fisher_test_apply(test = cancer_tmr_genes, universe = background_genes, query_list = msigdb_gene_set_list$h, return_overlap = F)
cancer_tmr_genes_hallmark_pathways_enrichment = dplyr::filter(cancer_tmr_genes_hallmark_pathways_enrichment, q_value < 0.05)
