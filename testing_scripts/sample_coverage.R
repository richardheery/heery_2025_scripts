# Calculate mean coverage of CpG sites per sample

# Load required packages
library(HDF5Array)
library(methodical)

# Load CPGEA and MCRPC meth RSEs
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")
mcrpc_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")

# Combine both RSEs
cpgea_rse = sort(cpgea_rse)
colData(cpgea_rse) = NULL
mcrpc_rse = sort(mcrpc_rse)
colData(mcrpc_rse) = NULL
combined_rse = cbind(cpgea_rse, mcrpc_rse)

# Extract 10,000 random CpG sites from chromosome 1 and set the coverage of not covered sites to 0
system.time({random_chr1_sites_rse = sampleMethSites(combined_rse, n_sites = 10000, seqnames_filter = "chr1")})
assay(random_chr1_sites_rse, 2)[is.na(assay(random_chr1_sites_rse, 2))] = 0

# Calculate the mean coverage of random sites
mean_coverage = colMeans(assay(random_chr1_sites_rse, 2))

### Calculate mean coverage for different region classes

# Set the coverage of not covered sites to 0
assay(combined_rse, 2)[is.na(assay(combined_rse, 2))] = 0

# Load genome annotation and filter for regions on chr1 and split into a list
genome_annotation = readRDS("complete_genome_annotation.rds")
genome_annotation_chr1 = genome_annotation[seqnames(genome_annotation) == "chr1"]
genome_annotation_chr1_list = split(genome_annotation_chr1, genome_annotation_chr1$region_type)

# Calculate the mean coverage for each region class in all samples. Took 17 minutes.
system.time({region_mean_coverage = lapply(genome_annotation_chr1_list, function(x)
  colMeans(assay(subsetByOverlaps(combined_rse, x), 2)))})
region_mean_coverage_df = data.frame(t(data.frame(region_mean_coverage)))

# CPGEA samples appear to have very low coverage of CpG islands (1X), predicted promoters (3X), CTCF binding sites (7x) 
# compared to MCRPC samples where coverage is still less than other regions but comparable
# LINEs, SINEs and LTRs are among the most covered regions in both CPGEA ands MCRPC samples 
rowMeans(select(region_mean_coverage_df, starts_with("D")))
rowMeans(select(region_mean_coverage_df, starts_with("N")))
rowMeans(select(region_mean_coverage_df, starts_with("T")))

# Calculate the mean methylation for each region class in all samples. Took 21 minutes.
system.time({region_mean_methylation = lapply(genome_annotation_chr1_list, function(x)
  colMeans(assay(subsetByOverlaps(combined_rse, x), 1), na.rm = T))})
region_mean_methylation_df = data.frame(t(data.frame(region_mean_methylation)))

rowMeans(select(region_mean_methylation_df, starts_with("D")))
rowMeans(select(region_mean_methylation_df, starts_with("N")))
rowMeans(select(region_mean_methylation_df, starts_with("T")))
