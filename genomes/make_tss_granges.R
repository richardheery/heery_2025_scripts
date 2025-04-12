# Create GRanges objects with TSS sites from gencode.v38.annotation.gff3

# Load required packages
library(GenomicFeatures)
library(dplyr)
library(rtracklayer)
library(genomeTools)

# Download GFF3 annotation file for Gencode version 38. Downloaded 2/6/2024
system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz")

# Read in gencode.v38.annotation.gff3.gz as a GRanges 
system.time({gencode_gr = rtracklayer::import.gff3("gencode.v38.annotation.gff3.gz")})

# Drop unnecessary metadata columns 
mcols(gencode_gr) = mcols(gencode_gr)[c("ID", "gene_id", "gene_name", "transcript_support_level", "tag", "type", "transcript_type")]

# Convert tag column to a character. Took 4 minutes. 
system.time({gencode_gr$tag = sapply(gencode_gr$tag, function(x) paste(x, collapse = "; "))})

# Remove transcript version but keep PAR_Y designation and add PAR_Y to gene name
gencode_gr$ID = gsub("\\.[0-9]*", "", gencode_gr$ID)
names(gencode_gr) = gencode_gr$ID

# Create GRanges objects for protein-coding and transcripts (excluding mitochondrial transcripts)
pc_transcripts_gr = gencode_gr[gencode_gr$type == "transcript" & 
    gencode_gr$transcript_type == "protein_coding" & seqnames(gencode_gr) != "chrM"]

# Identify MANE select protein-coding transcripts. Only protein-coding genes have MANE select transcripts. 
mane_pc_transcript_ids = pc_transcripts_gr[grepl("MANE_Select", pc_transcripts_gr$tag)]$ID

# Select desired metadata columns for GRanges
mcols(pc_transcripts_gr) = mcols(pc_transcripts_gr)[c("ID", "gene_id", "gene_name")]

# Save protein-coding and transcript GRanges and names of protein-coding transcripts
saveRDS(pc_transcripts_gr, "pc_transcripts_gr.rds")
saveRDS(mane_pc_transcript_ids, "mane_pc_transcript_ids.rds")

# Get TSS for protein-coding transcripts and save
pc_transcripts_tss_gr = resize(pc_transcripts_gr, width = 1, fix = "start")
saveRDS(pc_transcripts_tss_gr, "pc_transcripts_tss_gr.rds")
rm(gencode_gr); gc()

# Download prostate CAGE data from Fantom and read in as a GRanges
system("wget http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.tissue.hCAGE/prostate%252c%2520adult%252c%2520pool1.CNhs10628.10022-101D4.hg19.ctss.bed.gz")

# Read in prostate_adult_pool1_hg19_ctss.bed.gz as a GRanges
prostate_cage_hg19_gr = rtracklayer::import.bed("prostate, adult, pool1.CNhs10628.10022-101D4.hg19.ctss.bed.gz")

# Download liftover chain 
system("wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz")
system("gunzip hg19ToHg38.over.chain.gz")
hg19tohg38_chain =  rtracklayer::import.chain("hg19ToHg38.over.chain")

# Liftover prostate_cage_hg19_gr to hg38
prostate_cage_hg38_gr = liftover_granges(prostate_cage_hg19_gr, chain = hg19tohg38_chain)

# Add CAGE counts to pc_transcripts_tss_gr
# 40% of TSS supported by at least one CAGE tag. This drops to 30% for non-MANE transcripts.
pc_transcripts_tss_gr$cage_tags = 0
overlaps_df = data.frame(findOverlaps(pc_transcripts_tss_gr, prostate_cage_hg38_gr))
pc_transcripts_tss_gr$cage_tags[overlaps_df$queryHits] = prostate_cage_hg38_gr$score[overlaps_df$subjectHits]
cage_supported_gencode_tss = pc_transcripts_tss_gr[pc_transcripts_tss_gr$cage_tags >= 10]
saveRDS(cage_supported_gencode_tss, "cage_supported_gencode_tss.rds")

# Filter pc_transcripts_tss_gr for MANE transcripts. 
# 75% of MANE TSS supported by at least 1 CAGE tag. 52% greater than 10. 
mane_pc_tss = pc_transcripts_tss_gr[mane_pc_transcript_ids]

# Save MANE transcripts supported by more than 10 CAGE tags
cage_supported_mane_tss = mane_pc_tss[mane_pc_tss$cage_tags > 10]
saveRDS(cage_supported_mane_tss$ID, "cage_supported_mane_pc_transcript_ids.rds")