# Download methylation and RNA-seq data from TumourMethData

# Load required packages
library(TumourMethData)
library(DESeq2)

# Download WGBS data from CPGEA, MCRPC and TCGA along with TCGA array data
dir.create("methylation_data")
download_meth_dataset(dataset = "cpgea_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "mcrpc_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "tcga_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "tcga_450k_array_hg19", "methylation_data/")
download_meth_dataset(dataset = "roadmap_wgbs_hg38", "methylation_data/")

# Download transcript counts for CPGEA and MCRPC 
cpgea_transcript_counts_path = download_rnaseq_dataset(dataset = "cpgea_wgbs_hg38")
mcrpc_transcript_counts_path = download_rnaseq_dataset(dataset = "mcrpc_wgbs_hg38")

# Get names of protein-coding transcripts
pc_transcript_ids = readRDS("../auxillary_data/pc_transcripts_gr.rds")$ID

# Load counts for CPGEA
cpgea_transcript_counts = data.frame(data.table::fread(cpgea_transcript_counts_path), row.names = 1)

# Subset for protein-coding transcripts
cpgea_transcript_counts = cpgea_transcript_counts[pc_transcript_ids, ]

# Create a DESeqDataSet from Kallisto pcg counts. 
cpgea_dds = DESeqDataSetFromMatrix(countData = cpgea_transcript_counts, colData = data.frame(sample = names(cpgea_transcript_counts)), design = ~ 1)
cpgea_dds  = estimateSizeFactors(cpgea_dds) 
cpgea_normalized_kallisto_pcg_counts = data.frame(counts(cpgea_dds, normalized = T))
data.table::fwrite(tibble::rownames_to_column(cpgea_normalized_kallisto_pcg_counts, "transcript_id"), 
  "cpgea_normalized_kallisto_pcg_counts.tsv.gz", sep = "\t", row.names = F, quote = F, na = "NA")

# Load counts for MCRPC
mcrpc_transcript_counts = data.frame(data.table::fread(mcrpc_transcript_counts_path), row.names = 1)

# Subset for protein-coding transcripts
mcrpc_transcript_counts = mcrpc_transcript_counts[pc_transcript_ids, ]

# Create a DESeqDataSet from Kallisto pcg counts. 
mcrpc_dds = DESeqDataSetFromMatrix(countData = mcrpc_transcript_counts, colData = data.frame(sample = names(mcrpc_transcript_counts)), design = ~ 1)
mcrpc_dds  = estimateSizeFactors(mcrpc_dds) 
mcrpc_normalized_kallisto_pcg_counts = data.frame(counts(mcrpc_dds, normalized = T))
data.table::fwrite(tibble::rownames_to_column(mcrpc_normalized_kallisto_pcg_counts, "transcript_id"), 
  "mcrpc_normalized_kallisto_pcg_counts.tsv.gz", sep = "\t", row.names = F, quote = F, na = "NA")
