# Calculate tmr methylation in CPGEA samples for different tmr definitions and perform differential methylation for different definitions

# Load required packages
library(dplyr)
library(cowplot)
library(GenomicRanges)
library(SummarizedExperiment)
source("../auxillary_scripts/diff_meth_methylsig.R")
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/enrichment_tests.R")

### Test differential methylation of TMRs foind in different datasets

# Load TSS
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")
background_genes = unique(tss_gr$gene_name)

# Get list of TMRs
tmr_list = readRDS("../finding_tmrs/tmr_granges/tmr_list.rds")

# Combine TMRs into a single GRanges
combined_tmr_gr = unlist(GRangesList(tmr_list))
names(combined_tmr_gr) = paste(names(combined_tmr_gr), combined_tmr_gr$tmr_name, sep = "_")

# Load CPGEA meth RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(20)

# Perform differential analysis for different tmr definitions. Takes 15 minutes with 20 cores. 
system.time({tmr_diff_meth_results = diff_meth_methylsig(meth_rse = cpgea_meth_rse, genomic_regions = combined_tmr_gr, meth_reads_assay = "beta", coverage_assay = "cov",
  max_sites_per_chunk = floor(625000000/ncol(cpgea_meth_rse)), group_column = "condition", case = "Tumour", control = "Normal", BPPARAM = bpparam)})
saveRDS(tmr_diff_meth_results, "tmr_diff_meth_results_methylsig.rds")
tmr_diff_meth_results = data.frame(mcols(tmr_diff_meth_results))
tmr_diff_meth_results = readRDS("tmr_diff_meth_results_methylsig.rds")

# Add columns indicating dataset and direction
tmr_diff_meth_results$dataset = gsub("tmrs_.*", "tmrs", row.names(tmr_diff_meth_results))
tmr_diff_meth_results$direction = gsub(".*_tmrs_", "", gsub(".ENST.*", "", row.names(tmr_diff_meth_results)))

# Add the transcript and gene name for results
tmr_diff_meth_results$transcript = gsub("_tmr_.*", "", gsub(".*_ENST", "ENST", row.names(tmr_diff_meth_results)))
tmr_diff_meth_results$gene = tss_gr$gene_name[match(tmr_diff_meth_results$transcript, tss_gr$ID)]

# Denote whether promoters are hypermethylated, hypomethylated or unchanged 
tmr_diff_meth_results = mutate(tmr_diff_meth_results, 
  meth_change = case_when(
    fdr < 0.05 & meth_diff > 0 ~ "Hypermethylated",
    fdr < 0.05 & meth_diff < 0 ~ "Hypomethylated",
    fdr > 0.05 ~ "Unchanged"
    )
  )

# Separate negative and positive hypermethylated and hypomethylated results
tmr_diff_meth_results_grouped = list(
  hypermethylated_negative = filter(tmr_diff_meth_results, meth_change == "Hypermethylated" & direction == "negative"),
  hypomethylated_negative = filter(tmr_diff_meth_results, meth_change == "Hypomethylated" & direction == "negative"),
  hypermethylated_positive = filter(tmr_diff_meth_results, meth_change == "Hypermethylated" & direction == "positive"),
  hypomethylated_positive = filter(tmr_diff_meth_results, meth_change == "Hypomethylated" & direction == "positive")
)

# Get unique genes for each TMR group
tmr_diff_meth_genes = lapply(tmr_diff_meth_results_grouped, function(x)
  split(x$gene, x$dataset))
tmr_diff_meth_genes = lapply(tmr_diff_meth_genes, function(x)
  lapply(x, function(y) unique(y)))

# Convert meth_change to a factor
tmr_diff_meth_results$meth_change = factor(tmr_diff_meth_results$meth_change, levels = c("Hypomethylated", "Unchanged", "Hypermethylated"))

# Remove rows with missing values
tmr_diff_meth_results = tmr_diff_meth_results[complete.cases(tmr_diff_meth_results), ]

# Get proportion of hypermethylated, hypomethylated and unchanged TMRs for each definition
tmr_diff_meth_results_summary = mutate(
  summarize(group_by(tmr_diff_meth_results, dataset, direction, meth_change), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated and unchanged tmrs for each definition
tmr_direction_meth_change_barplot = ggplot(tmr_diff_meth_results_summary, 
  aes(x = dataset, y = count, fill = meth_change, label = paste0(round(freq, 2)*100, "%"))) +
  geom_col(position = "dodge", color = "black") + 
  geom_text(mapping = aes(x = dataset, y = count + 100, group = meth_change), position = position_dodge(width = 0.9), size = 4) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 24), 
  	axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
  	legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"),
    strip.text.x = element_text(size = 20)) +
  scale_x_discrete(labels = c("Normal\nProstate", "Prostate\nTumours", "Prostate\nMetastases")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = c("#4B878BFF", "grey", "#D01C1FFF")) +
  labs(x = "TMR Source", y = "Number of TMRs", fill = NULL,  title = NULL) +
  facet_wrap(as.formula(paste("~", "direction")), nrow = 2, 
    labeller = as_labeller(setNames(c("Negative TMRs", "Positive TMRs"), c("negative", "positive"))))
tmr_direction_meth_change_barplot
saveRDS(tmr_direction_meth_change_barplot, "tmr_direction_meth_change_barplot.rds")

### Create enrichment plots

# Load MSigDB gene sets
msigdb_gene_set_list = readRDS("../auxillary_data/msigdb_complete_gene_set_list.rds")

# Make a function which will take a list of genes and perform pathway enrichment
plot_tmr_gene_enrichment = function(gene_list, pathways, filter_groups = NULL, 
  facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), ylab = "Pathway", title = NULL){
  
  # Enrichment results
  enrichment_results = lapply(gene_list, function(x) 
    fisher_test_apply(test = x, universe = background_genes, query_list = pathways, return_overlap = F))

  # Combine enrichment_results into a single table
  combined_enrichment_results = bind_rows(c(enrichment_results), .id = "group")

  # Filter for significant results
  combined_enrichment_results = filter(combined_enrichment_results, q_value < 0.05)
  combined_enrichment_results = mutate(group_by(combined_enrichment_results, group), 
    rank = rank(-relative_enrichment), ranking_name = paste0(query_name, rank(-relative_enrichment)))

  # Rank the results separately for each dataset
  combined_enrichment_results$group = factor(combined_enrichment_results$group, unique(combined_enrichment_results$group))
  combined_enrichment_results = arrange(combined_enrichment_results, desc(rank), group)
  combined_enrichment_results$ranking = factor(combined_enrichment_results$ranking_name, unique(combined_enrichment_results$ranking_name))
  
  if(nrow(combined_enrichment_results) == 0){
    return(NULL)
  }
  
  if(!is.null(filter_groups)){
    combined_enrichment_results = dplyr::filter(combined_enrichment_results, group %in% filter_groups)
  }

  combined_enrichment_plot = ggplot(combined_enrichment_results, 
    aes(y = ranking, x = relative_enrichment, color = q_value, size = test_overlap_size)) + 
      geom_point() + geom_vline(xintercept = 1, linetype = "dashed")
  combined_enrichment_plot = customize_ggplot_theme(combined_enrichment_plot, title = title, xlab = "Relative Enrichment", ylab = ylab, color_title = "Adjusted\np-value",
    facet = "group", facet_labels = facet_labels, facet_scales = "free_y", strip_text_size = 20) + 
    theme(axis.text.y = element_text(size = 14), strip.background = element_blank()) +
    scale_y_discrete(labels = function(x) gsub("[0-9.]*$", "", gsub("_", " ", x))) +
    labs(size = "Overlap Size") +
    scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T))

  combined_enrichment_plot
  
}

# Make plots of enriched KEGG and Hallmark pathways in prostate tumours and metastases
system.time({kegg_enrichment_plots_cancer = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$`CP:KEGG`, 
    filter_groups = c("cpgea_normal_tmrs", "cpgea_tumour_tmrs"), facet_labels =  c("Normal Prostate TMRs", "Prostate Tumours TMRs"), ylab = "KEGG Pathway"))})
system.time({hallmark_enrichment_plots_cancer = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$h, 
    filter_groups = c("cpgea_normal_tmrs", "cpgea_tumour_tmrs", "mcrpc_tmrs"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumours TMRs", "Prostate Metastases TMRs"), ylab = "MSigDB Hallmark Pathway"))})

# Combine KEGG and Hallmark enrichment plots
combined_kegg_and_hallmark_plots = ggpubr::ggarrange(plotlist = list(
  customize_ggplot_theme(kegg_enrichment_plots_cancer[[1]], title = "Hypermethylated Negative TMRs", xlab = "Relative Enrichment", ylab = "KEGG Pathway") + theme(strip.background = element_blank()), 
  customize_ggplot_theme(hallmark_enrichment_plots_cancer[[2]], title = "Hypomethylated Negative TMRs", xlab = "Relative Enrichment", ylab = "MSigDB Hallmark Pathway") + theme(strip.background = element_blank())),  
    nrow = 2, labels = c("A", "B"))
combined_kegg_and_hallmark_plots
ggsave(plot = combined_kegg_and_hallmark_plots, "../figures/supp_figure15.pdf", width = 27, height = 18)

### Compare methylation of DMRs and TMRs

# Download DMRs
system("wget wget https://s3-us-west-2.amazonaws.com/feng.genomics/WGBS/secondary_data.zip")
system("tar -xvzf secondary_data.zip")

# Read in DMRs and convert to a GRanges and set chromosome style to UCSC
dmrs = data.table::fread("secondary_data/DSS/DSS_mAdeno_benignprostate.tsv")
dmrs_gr = makeGRangesFromDataFrame(dmrs, keep.extra.columns = T)
seqlevels(dmrs_gr)[c(23, 24)] = c("X", "Y")
seqlevelsStyle(dmrs_gr) = "UCSC"

# Check proportion of hypermethylated and hypomethylated DMRs. 96% are hypomethylated
prop.table(table(sign(dmrs_gr$diff.Methy)))

# Check proportion of hypermethylated and hypomethylated DMRs overlapping negative TMRs in metastasis samples. 22% are hypermethylated 
prop.table(table(sign(subsetByOverlaps(dmrs_gr, tmr_list$mcrpc_tmrs_negative)$diff.Methy)))

mcrpc_tmrs = readRDS("../finding_tmrs/tmr_granges/mcrpc_tmrs.rds")

tmr_genes_overlapping_dmrs = unique(subsetByOverlaps(mcrpc_tmrs, dmrs_gr)$gene_name)

filter(fisher_test_apply(test = tmr_genes_overlapping_dmrs, universe = background_genes, query_list = msigdb_gene_set_list$`CP:KEGG`, return_overlap = F), q_value < 0.05)
filter(fisher_test_apply(test = tmr_genes_overlapping_dmrs, universe = background_genes, query_list = msigdb_gene_set_list$h, return_overlap = F), q_value < 0.05)
