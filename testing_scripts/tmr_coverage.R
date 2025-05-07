# Evalute TMRs

#
library(SummarizedExperiment)
library(methodical)

# Load repeats
repeats = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

# Load TMRs
cpgea_tumour_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs.rds")
cpgea_tumour_tmrs_50kb = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs.rds")

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
