# Download all STAR count files from TCGA and create TPM, FPKM and FPKM-UQ tables from them

# Load required packages
library(doParallel)

# Make cluster
cl = makeCluster(11)
registerDoParallel(cl)

# Load manifest for TCGA
tcga_star_counts_manifest = read.csv("gdc_manifest_all_tcga_star_counts_2023_11_09.txt", sep = "\t", header = T)

# Add URL by prepending https://api.gdc.cancer.gov/data/ to id
tcga_star_counts_manifest$url = paste0("https://api.gdc.cancer.gov/data/", tcga_star_counts_manifest$id)

# Download all files
# Took 40 minutes with 20 cores
system.time({foreach(file = tcga_star_counts_manifest$url) %dopar% {
  
  download.file(file, destfile = paste("all_star_count_files", basename(file), sep = "/"))
  
}})

# Get paths to all files
star_count_files = list.files("all_star_count_files", full.names = T)

# Extract first columns from each file. Took 11 minutes. 
system.time({first_columns = foreach(file = star_count_files, .combine = data.frame) %dopar% {
  
  data.table::fread(file, nThread = 1)[[1]]
  
}})

# Check that the first column is the same for all files. All columns duplicated
table(duplicated(as.list(first_columns)))

# Get the gene IDs from the first column of one file
gene_ids = data.table::fread(star_count_files[1])[[1]]
gene_ids = grep("ENSG", gene_ids, value = T)

# Remove version ID while keeping PAR information
gene_ids = gsub("\\.[0-9]*", "", gene_ids)

# Read in the JSON file associated with the manifest as a list and turn it into a data.frame
json_file = rjson::fromJSON(file = "gdc_json_all_tcga_star_counts_2023_11_09.json")
json_file = split(unname(unlist(json_file)), names(unlist(json_file)))
json_file_df = data.frame(json_file[setdiff(names(json_file), "annotations.annotation_id")])

# Remove TCGA- from project_id
json_file_df$cases.project.project_id = gsub("TCGA-", "", json_file_df$cases.project.project_id)

# Add the TCGA barcode using the file ID and then add submitter and sample type. 
tcga_star_counts_manifest$barcode = TCGAutils::UUIDtoBarcode(tcga_star_counts_manifest$id, from_type = c("file_id"))$associated_entities.entity_submitter_id
tcga_star_counts_manifest$submitter = substr(tcga_star_counts_manifest$barcode, start = 1, stop = 12)
tcga_star_counts_manifest$sample_type = substr(tcga_star_counts_manifest$barcode, start = 14, stop = 15)

# Add the project ID to the tcga_star_counts_manifest from the JSON file
tcga_star_counts_manifest$project = json_file_df$cases.project.project_id[match(tcga_star_counts_manifest$filename, json_file_df$file_name)]

# Create column combining submitter ID and sample type
tcga_star_counts_manifest$submitter_and_sample_type = paste(tcga_star_counts_manifest$submitter, tcga_star_counts_manifest$sample_type, sep = "_")
tcga_star_counts_manifest$submitter_and_sample_type = gsub("-", "_", tcga_star_counts_manifest$submitter_and_sample_type)

# Sort table by project, submitter and barcode and select the first sample for any duplicated entries for submitter_and_sample_type
tcga_star_counts_manifest = dplyr::arrange(tcga_star_counts_manifest, project, submitter, barcode)
tcga_star_counts_manifest = tcga_star_counts_manifest[!duplicated(tcga_star_counts_manifest$submitter_and_sample_type), ]

# Prepend all_star_count_files/ to manifest$id in order to get file path
tcga_star_counts_manifest$id = paste0("all_star_count_files/", tcga_star_counts_manifest$id)

# Split the table into a list for each project
tcga_star_counts_manifest_list = split(tcga_star_counts_manifest, tcga_star_counts_manifest$project)

# Loop through the manifest table for each project 
# and create a table with the TPM values for each project. Took 80 minutes
system.time({foreach(project = names(tcga_star_counts_manifest_list)) %do% {
  
  # Loop through each file for the project in read the TPM values as a data.frame
  tpm_table = foreach(
    file = tcga_star_counts_manifest_list[[project]]$id, 
    sample = tcga_star_counts_manifest_list[[project]]$submitter_and_sample_type,
    .combine = data.frame) %do% {
      
      # Get TPM values for sample from file
      setNames(data.table::fread(file)[5:60664, "tpm_unstranded"], sample)
      
    }
  
  # Add gene IDs to tpm_table and set as a column
  row.names(tpm_table) = gene_ids
  tpm_table = tibble::rownames_to_column(tpm_table, "gene_id")
  
  # Write table to tpm_tables
  data.table::fwrite(tpm_table, file = paste0("tpm_tables/", project, "_tpm_values.tsv.gz"), 
    sep = "\t", col.names = T, row.names = F, quote = F, na = "NA")
  
}})

# Loop through the manifest table for each project 
# and create a table with the raw counts for each project. Took 10 minutes with 11 cores.
system.time({foreach(project = names(tcga_star_counts_manifest_list), .packages = "foreach") %dopar% {
  
  # Loop through each file for the project in read the counts values as a data.frame
  counts_table = foreach(
    file = tcga_star_counts_manifest_list[[project]]$id, 
    sample = tcga_star_counts_manifest_list[[project]]$submitter_and_sample_type,
    .combine = data.frame) %do% {
      
      # Get counts values for sample from file
      setNames(data.table::fread(file)[5:60664, "unstranded"], sample)
      
    }
  
  # Add gene IDs to counts_table and set as a column
  row.names(counts_table) = gene_ids
  counts_table = tibble::rownames_to_column(counts_table, "gene_id")
  
  # Write table to raw_counts_tables
  data.table::fwrite(counts_table, file = paste0("raw_counts_tables/", project, "_raw_counts.tsv.gz"), 
    sep = "\t", col.names = T, row.names = F, quote = F, na = "NA")
  
}})
  
# Get path to all count tables
count_files = list.files("raw_counts_tables", full.names = T)
names(count_files) = gsub("_raw_counts.tsv.gz", "", basename(count_files))

# Loop through each counts table and normalize counts with DESeq2 and save tables
# Took 8 minutes
library(DESeq2)
library(foreach)
system.time({foreach(project = names(count_files)) %do% {
  
  # Read in counts table for project
  counts_table = data.frame(data.table::fread(count_files[project]), row.names = 1)
  
  # Normalize counts using DESeq2
  counts_table_dds = DESeqDataSetFromMatrix(countData = counts_table, colData = data.frame(sample = names(counts_table)), design = ~ 1)
  counts_table_dds  = estimateSizeFactors(counts_table_dds) 
  deseq2_normalized_counts = data.frame(counts(counts_table_dds, normalized = T))
  
  # Set row.names as a column called gene_id
  deseq2_normalized_counts = tibble::rownames_to_column(deseq2_normalized_counts, "gene_id")
  
  # Save table of normalized counts
  data.table::fwrite(deseq2_normalized_counts, file = paste0("deseq_normalized_counts/", project, "_deseq_normalized_counts.tsv.gz"), 
    sep = "\t", col.names = T, row.names = F, quote = F, na = "NA")
  
}})
