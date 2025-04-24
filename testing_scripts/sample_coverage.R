# Calculate mean coverage of CpG sites per sample

# Load required packages
library(HDF5Array)
library(methodical)

# Load CPGEA and MCRPC meth RSEs
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("cpgea_meth_rse/")
mcrpc_rse = HDF5Array::loadHDF5SummarizedExperiment("mcrpc_meth_rse/")

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
