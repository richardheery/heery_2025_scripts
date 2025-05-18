# Create a GRanges with protein-coding MANE transcripts from Gencode 36 for hg19

# Download Gencode 36 annotation for hg19
system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh37_mapping/gencode.v36lift37.annotation.gff3.gz")

# Load basic annotation for Gencode 36
gencode_annotation = rtracklayer::import.gff3("gencode.v36lift37.annotation.gff3.gz")

# Filter for protein-coding transcripts
gencode_transcript_annotation = gencode_annotation[gencode_annotation$type == "transcript" & gencode_annotation$transcript_type == "protein_coding"]

# Convert tag into a character vector
gencode_transcript_annotation$tag = lapply(gencode_transcript_annotation$tag, function(x) paste(x, collapse = ", "))

# Filter for MANE transcripts. There are 17,793 MANE protein-coding transcripts
gencode_transcript_annotation_mane = gencode_transcript_annotation[grepl("MANE_Select", gencode_transcript_annotation$tag)]
gencode_transcript_annotation_mane$ID = gsub("\\.[0-9]*", "", gencode_transcript_annotation_mane$ID)

# Remove version ID from gene ID
gencode_transcript_annotation_mane$gene_id = gsub("\\..*", "", gencode_transcript_annotation_mane$gene_id)

# Add PAR_Y to gene IDs on Y PAR.
gencode_transcript_annotation_mane$gene_id[grepl("PAR", gencode_transcript_annotation_mane$ID)] =
  paste0(gencode_transcript_annotation_mane$gene_id[grepl("PAR", gencode_transcript_annotation_mane$ID)], "_PAR_Y")

# Set names of gencode_transcript_annotation_mane as gene_id
names(gencode_transcript_annotation_mane) = gencode_transcript_annotation_mane$gene_id

# Save gencode_transcript_annotation_mane
saveRDS(gencode_transcript_annotation_mane, "gencode_36_mane_transcripts_hg19_gr.rds")