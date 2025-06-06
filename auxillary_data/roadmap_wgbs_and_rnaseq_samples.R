# Download matching Roadmap WGBS and RNA-seq data and create gene counts table and a HDF5-backed SummarizedExperiment for the WGBS data 

# Load required packages
library(methrix)
library(dplyr)
library(doParallel)

# Load experiment_report as a data.frame and adjust column names
experiment_report = data.frame(data.table::fread("experiment_report_2025_6_5_13h_12m.tsv", skip = 1))
names(experiment_report) = tolower(gsub("\\.", "_", names(experiment_report)))

# Remove columns with the same value for each row
experiment_report = experiment_report[, names(which(sapply(experiment_report, function(x) length(unique(x))) > 1))]

# Filter for WGBS and RNA-seq experiments
wgbs_experiments = filter(experiment_report, assay_title == "WGBS")
rna_experiments = filter(experiment_report, assay_title == "polyA plus RNA-seq")

# Get sample accessions for WGBS and RNA-seq files
wgbs_accessions = wgbs_experiments$biosample_accession
rna_accessions = rna_experiments$biosample_accession

# Find sample accessions with both WGBS and RNA-seq data
matching_accessions = intersect(wgbs_accessions, rna_accessions)

# Filter for experiments with both WGBS and RNA-seq data
matching_data_experiments = filter(experiment_report, biosample_accession %in% matching_accessions)

# Remove /files/ and trailing / from file names
matching_data_experiments$files = gsub("/files/", "", matching_data_experiments$files)
matching_data_experiments$files = gsub("/$", "", matching_data_experiments$files)

# Create a list matching files accessions to sample accessions
files_list = strsplit(matching_data_experiments$files, "/,")
names(files_list) = matching_data_experiments$biosample_accession

# Turn files list inside out
files_list = with(reshape2::melt(files_list), split(L1, value))

# Get file accessions for experiments
wgbs_files = unlist(strsplit(filter(matching_data_experiments, assay_title == "WGBS")$files, "/,"))
rnaseq_files = unlist(strsplit(filter(matching_data_experiments, assay_title == "polyA plus RNA-seq")$files,"/,"))

# Read in URLs to download all WGBS and RNA-seq files
file_urls = readLines("roadmap_wgbs_and_rnaseq_urls.txt")[-1]

# Create a data.frame with the URLs and their associated accessiona and file type
url_experiment_report = data.frame(
  file_urls = file_urls,
  file_accession = gsub("\\..*", "", basename(file_urls)),
  file_extension = tools::file_ext(gsub(".gz", "", file_urls))
)

# Filter for files from experiments with both WGBS and RNA-seq data
url_experiment_report = filter(url_experiment_report, file_accession %in% c(wgbs_files, rnaseq_files))

# Add a column indicating if files are for WGBS or RNA-seq data
url_experiment_report$assay = NA
url_experiment_report$assay[url_experiment_report$file_accession %in% wgbs_files] = "WGBS"
url_experiment_report$assay[url_experiment_report$file_accession %in% rnaseq_files] = "RNASeq"

# Get URLs for RNA-Seq TSV files with gene counts 
rnaseq_file_urls = filter(url_experiment_report, assay == "RNASeq", file_extension == "tsv")$file_url

# Get URLs for WGBS BED files
wgbs_file_urls = filter(url_experiment_report, assay == "WGBS", file_extension == "bed")$file_url

# Create directories for WGBS and RNA-seq files
dir.create("roadmap_wgbs_bed_files")
dir.create("roadmap_rnaseq_txt_files")

# Make a cluster
cl = makeCluster(19)
registerDoParallel(cl)

# Download RNA-seq files. Took 90 seconds with 19 cores
download_rnaseq_file_function = function(x){system(paste("wget -P roadmap_rnaseq_txt_files", x))}
system.time({foreach::foreach(url = rnaseq_file_urls) %dopar% {download_rnaseq_file_function(url)}})

# Download WGBS files. Took 6 hours
download_wgbs_file_function = function(x){system(paste("wget -P roadmap_wgbs_bed_files", x))}
system.time({foreach::foreach(url = wgbs_file_urls) %dopar% {download_wgbs_file_function(url)}})

### Create HDF5-backed WGBS array

# Get paths to roadmap BED files and name them with their file accession
roadmap_wgbs_files = list.files("roadmap_wgbs_bed_files", full.names = T)
names(roadmap_wgbs_files) = gsub("\\..*", "", basename(roadmap_wgbs_files))

# Get experiment data for WGBS files and select desired columns
wgbs_data_experiments = filter(matching_data_experiments, assay_title == "WGBS")
wgbs_data_experiments = select(wgbs_data_experiments, sample_accession = "biosample_accession", 
  tissue = "biosample_term_name", age = "biosample_age", description = "biosample_summary")

# Create a data.frame with metadata for WGBS files
wgbs_metadata = data.frame(file_accession = names(roadmap_wgbs_files))
row.names(wgbs_metadata) = unlist(files_list[wgbs_metadata$file_accession])
wgbs_metadata = wgbs_metadata[wgbs_data_experiments$sample_accession, , drop = F]
wgbs_metadata[, c("tissue", "age", "description")] = wgbs_data_experiments[, c("tissue", "age", "description")]

# Add a column indicating sex and remove description
wgbs_metadata$sex = ifelse(grepl("female", wgbs_metadata$description), "female", "male")
wgbs_metadata$description = NULL

# Put roadmap_wgbs_files in same order as wgbs_metadata
roadmap_wgbs_files = roadmap_wgbs_files[wgbs_metadata$file_accession]

# Get reference CpG sites for hg38
hg38_cpgs = methrix::extract_CPGs("BSgenome.Hsapiens.UCSC.hg38")

# Create methrix object from CpG methylation files. Took 1.5 hours. 
system.time({roadmap_methrix = methrix::read_bedgraphs(roadmap_wgbs_files, ref_cpgs = hg38_cpgs, coldata = wgbs_metadata, 
  chr_idx = 1, start_idx = 2, strand_idx = 6, cov_idx = 10, beta_idx = 11, 
  zero_based = T, stranded = T, collapse_strands = T, h5 = T, h5_dir = "roadmap_methrix_h5")})

# Convert methrix into a meth RSE and resave
roadmap_methrix = methrix::load_HDF5_methrix("roadmap_methrix_h5/")
system.time({roadmap_meth_rse = methodical::methrixToRSE(roadmap_methrix)})
system.time({HDF5Array::quickResaveHDF5SummarizedExperiment(roadmap_meth_rse)})
file.rename("roadmap_methrix_h5", "../methylation_data/roadmap_meth_rse_hg38/")

### Create table with gene counts

# Get paths to RNA-seq files
roadmap_rnaseq_files = list.files("roadmap_rnaseq_txt_files", full.names = T)
names(roadmap_rnaseq_files) = gsub("\\..*", "", basename(roadmap_rnaseq_files))

# Read in the expected count (the 5th column) from each file
gene_counts = lapply(roadmap_rnaseq_files, function(x)
  data.frame(data.table::fread(x, select = c(1, 5)), row.names = 1))

# Conver gene_counts into a data.frame
gene_counts = setNames(data.frame(gene_counts), names(roadmap_rnaseq_files))

# Rename columns with sample names
names(gene_counts) = files_list[names(gene_counts)]

# Remove version from gene IDs
row.names(gene_counts) = gsub("\\.[0-9]*", "", row.names(gene_counts))

# Save table
data.table::fwrite(tibble::rownames_to_column(gene_counts, "gene_id"), 
  "rnaseq_data/roadmap_gene_estimated_gene_counts.tsv.gz", sep = "\t", quote = F)

# Create supplementary table with file details

rnaseq_files_df = data.frame(
  file_accession = tools::file_path_sans_ext(basename(roadmap_rnaseq_files)),
  biosample_accession = unlist(files_list[names(roadmap_rnaseq_files)]),
  tissue = wgbs_metadata$tissue[match(unlist(files_list[names(roadmap_rnaseq_files)]), row.names(wgbs_metadata))],
  assay_title = "polyA plus RNA-seq",  
  row.names = NULL
)

wgbs_files_df = wgbs_metadata
wgbs_files_df = tibble::rownames_to_column(wgbs_files_df, "biosample_accession")
wgbs_files_df = select(wgbs_files_df, file_accession, biosample_accession, tissue)
wgbs_files_df$assay_title = "WGBS"

supp_table = rbind(wgbs_files_df, rnaseq_files_df)
write.table(supp_table, "roadmap_supp_table.tsv", sep = "\t", quote = F, row.names = F)
