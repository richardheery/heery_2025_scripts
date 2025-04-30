# Calculate promoter methylation in CPGEA and MCRPC samples for different promoter definitions 

# Load required packages
library(methodical)
library(dplyr)

# Get promoter definition list
promoter_definition_list = readRDS("promoter_definition_list.rds")

# Load CPGEA methylation RSE
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(3)

# Make a directory to store result files
dir.create("cpgea_promoter_definition_methylation_tables")

# Get methylation values for all promoter methylation definitions and save to cpgea_promoter_definition_methylation_tables. Took 20 minutes. 
system.time({for(definition in names(promoter_definition_list)){
  message(definition)
  promoter_definition_methylation = methodical::summarizeRegionMethylation(
    meth_rse = cpgea_rse, genomic_regions = promoter_definition_list[[definition]], 
    BPPARAM = bpparam)
  promoter_definition_methylation = tibble::rownames_to_column(promoter_definition_methylation, "transcript_id")
  data.table::fwrite(promoter_definition_methylation, 
    file = sprintf("cpgea_promoter_definition_methylation_tables/%s_definition_methylation_table.tsv.gz", definition), 
    sep = "\t", row.names = F, na = "NA", quote = F)
}})

# Get paths to all promoter methylation definition tables
promoter_definition_methylation_table_files = list.files("cpgea_promoter_definition_methylation_tables", full.names = T, pattern = "_methylation_table.tsv.gz")
names(promoter_definition_methylation_table_files) = gsub("_definition_methylation_table.tsv.gz", "", basename(promoter_definition_methylation_table_files))

# Read in tables as a list and save
cpgea_promoter_definition_methylation_tables = lapply(promoter_definition_methylation_table_files, function(x)
  data.frame(data.table::fread(x), row.names = 1))

# Save list of promoter methylation tables
saveRDS(cpgea_promoter_definition_methylation_tables, "cpgea_promoter_definition_methylation_tables/cpgea_promoter_definition_methylation_tables.rds")

### Repeat for MCRPC samples

# Load mcrpc methylation RSE
mcrpc_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(3)

# Make a directory to store result files
dir.create("mcrpc_promoter_definition_methylation_tables")

# Get methylation values for all promoter methylation definitions and save to mcrpc_promoter_definition_methylation_tables. Took 20 minutes. 
system.time({for(definition in names(promoter_definition_list)){
  message(definition)
  promoter_definition_methylation = methodical::summarizeRegionMethylation(
    meth_rse = mcrpc_rse, genomic_regions = promoter_definition_list[[definition]], 
    BPPARAM = bpparam)
  promoter_definition_methylation = tibble::rownames_to_column(promoter_definition_methylation, "transcript_id")
  data.table::fwrite(promoter_definition_methylation, 
  file = sprintf("mcrpc_promoter_definition_methylation_tables/%s_definition_methylation_table.tsv.gz", definition), 
  sep = "\t", row.names = F, na = "NA", quote = F)
}})

# Get paths to all promoter methylation definition tables
promoter_definition_methylation_table_files = list.files("mcrpc_promoter_definition_methylation_tables", full.names = T, pattern = "_methylation_table.tsv.gz")
names(promoter_definition_methylation_table_files) = gsub("_definition_methylation_table.tsv.gz", "", basename(promoter_definition_methylation_table_files))

# Read in tables as a list and save
mcrpc_promoter_definition_methylation_tables = lapply(promoter_definition_methylation_table_files, function(x)
  data.frame(data.table::fread(x), row.names = 1))

# Save list of promoter methylation tables
saveRDS(mcrpc_promoter_definition_methylation_tables, "mcrpc_promoter_definition_methylation_tables/mcrpc_promoter_definition_methylation_tables.rds")