# Calcualte number of TMRs found with different numbers of samples

# Load required packages
library(methodical)
library(dplyr)
library(foreach)
library(plotR)

# Get TSS Granges
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")[names(tss_gr)]

# Expand transcripts_gr 5 KB upstream and downstream
transcripts_gr = methodical:::expand_granges(transcripts_gr, 5000, 5000)

# Load CPGEA methylation RSE and transcript counts
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/cpgea_normalized_kallisto_pcg_counts.tsv.gz"), row.names = 1)

# Mask CpGs with less than 10 reads covering them
assay(cpgea_meth_rse, 2)[is.na(assay(cpgea_meth_rse, 2))] = 0
assay(cpgea_meth_rse, 1)[assay(cpgea_meth_rse, 2) < 10] = NA

# Get CPGEA normal and tumour samples
normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)

# Create a bpparam object
bpparam = BiocParallel::SnowParam(workers = 8)

# Create 10 random samples each for sample sizes of 20, 40, 60, 80 and 100
set.seed(123)
random_samples = setNames(lapply(seq(20, 100, 20), function(x) 
  setNames(lapply(1:10, function(y) sample(normal_samples, x)), paste0("_", 1:10))), seq(20, 100, 20))
random_samples = unlist(random_samples, recursive = F)

system.time({random_sample_tmrs = foreach(sample_subset = random_samples, i = names(random_samples), .packages = "methodical") %do% {
  
  print(i)
  
  # Calculate correlations using sample
  sample_subset_cors = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
    transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = sample_subset, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
    cor_method = "spearman", min_number_complete_pairs = 20, BPPARAM = bpparam, add_distance_to_region = T)
  
  # Find TMRs for cors
  sample_subset_tmrs = findTMRs(correlation_list = sample_subset_cors, 
    p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)
  
  sample_subset_tmrs
  
}})
names(random_sample_tmrs) = names(random_samples)
saveRDS(random_sample_tmrs, "random_sample_tmrs.rds")
random_sample_tmrs = readRDS("random_sample_tmrs.rds")

# Make a data.frame with the number of TMRs in each random sample
random_sample_df = data.frame(
  sample_name = gsub("\\..*", "", names(random_sample_tmrs)),
  number_tmrs = lengths(random_sample_tmrs),
  row.names = NULL
)
random_sample_df$sample_name = factor(random_sample_df$sample_name, levels = as.character(seq(20, 200, 20)))

# Make boxplots of the number of TMRs per sample group
sample_boxplots = ggplot(random_sample_df, aes(x = sample_name, y = number_tmrs)) +
  geom_boxplot(fill = "#A28CB1") +
  geom_point()
sample_boxplots = customize_ggplot_theme(sample_boxplots, xlab = "Sample Number", ylab = "TMR Number")
sample_boxplots
ggsave(plot = sample_boxplots, "../figures/supp_figure6.pdf", width = 16, height = 9)
