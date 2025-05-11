# Evalute TMRs

#
library(SummarizedExperiment)
library(methodical)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/tmr_plot_functions.R")

# Load repeats
repeats = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

# Load TMRs
cpgea_normal_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_normal_tmrs.rds")
cpgea_normal_tmrs_50kb = readRDS("../finding_tmrs/tmr_granges/cpgea_normal_tmrs_50kb.rds")
cpgea_tumour_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs.rds")
cpgea_tumour_tmrs_50kb = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs_50kb.rds")
mcrpc_tmrs = readRDS("../finding_tmrs/tmr_granges/mcrpc_tmrs.rds")
mcrpc_tmrs_50kb = readRDS("../finding_tmrs/tmr_granges/mcrpc_tmrs_50kb.rds")

# Subset for TMRs overlapping repeats
cpgea_tumour_tmrs_repeats = subsetByOverlaps()

# Load CPGEA meth RSE and subset for cpg_sites
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")

# Change NA coverage values to 0
assay(cpgea_meth_rse, 2)[is.na(assay(cpgea_meth_rse, 2))] = 0

bpparam = BiocParallel::MulticoreParam(3)

# Get mean coverage of all TMRs. Took 1 minuite
system.time({cpgea_tumour_tmr_coverage = summarizeRegionMethylation(meth_rse = cpgea_meth_rse, genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, assay = 2, BPPARAM = bpparam)})

cpgea_tumour_tmr_coverage_mean = rowMeans(cpgea_tumour_tmr_coverage, na.rm = T)

# Check negaitve and positive TMR coverage
fivenum(cpgea_tumour_tmr_coverage_mean[cpgea_tumour_tmrs$direction == "Negative"])
fivenum(cpgea_tumour_tmr_coverage_mean[cpgea_tumour_tmrs$direction == "Positive"])

# Check repeat overlap coverage
fivenum(cpgea_tumour_tmr_coverage_mean[cpgea_tumour_tmrs %over% repeats])

### Plot distributions of tmrs

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

test_list = list(
  normal_5kb = cpgea_normal_tmrs,
  tumour_5kb = cpgea_tumour_tmrs,
  metastasis_5kb = mcrpc_tmrs,
  normal_50kb = cpgea_normal_tmrs_50kb,
  tumour_50kb = cpgea_tumour_tmrs_50kb,
  metastasis_50kb = mcrpc_tmrs_50kb
)

lapply(test_list, function(x) prop.table(table(subsetByOverlaps(x, repeats)$direction)))

test_list_negative = list(
  normal_5kb = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Negative"],
  tumour_5kb = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Negative"],
  metastasis_5kb = mcrpc_tmrs[mcrpc_tmrs$direction == "Negative"],
  normal_50kb = cpgea_normal_tmrs_50kb[cpgea_normal_tmrs_50kb$direction == "Negative"],
  tumour_50kb = cpgea_tumour_tmrs_50kb[cpgea_tumour_tmrs_50kb$direction == "Negative"],
  metastasis_50kb = mcrpc_tmrs_50kb[mcrpc_tmrs_50kb$direction == "Negative"]
)

test_list_positive = list(
  normal_5kb = cpgea_normal_tmrs[cpgea_normal_tmrs$direction == "Positive"],
  tumour_5kb = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Positive"],
  metastasis_5kb = mcrpc_tmrs[mcrpc_tmrs$direction == "Positive"],
  normal_50kb = cpgea_normal_tmrs_50kb[cpgea_normal_tmrs_50kb$direction == "Positive"],
  tumour_50kb = cpgea_tumour_tmrs_50kb[cpgea_tumour_tmrs_50kb$direction == "Positive"],
  metastasis_50kb = mcrpc_tmrs_50kb[mcrpc_tmrs_50kb$direction == "Positive"]
)
test_list_positive = lapply(test_list_positive, unname)
test_list_positive_gr = unlist(GRangesList(test_list_positive))
