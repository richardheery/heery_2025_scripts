# Make a GRanges object with location of repeats from UCSC RepeatMasker annotation

# Load required packages
library(dplyr)
library(GenomicRanges)

# Download UCSC repeatmasker annotation.
system("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz")

# Read in downloaded RepeatMasker table
repeatmasker_table = data.table::fread("hg38.fa.out.gz", header = F, fill = T, skip = 3)

# Name the columns based on the description in webrepeatmaskerhelp.html (downloaded from https://www.repeatmasker.org/webrepeatmaskerhelp.html)
names(repeatmasker_table) = c("sw_score", "substitution_pc", "deletion_pc", "insertion_pc",
  "seqnames", "start", "end", "bp_to_end_of_seq", "strand", "name", "taxonomy", "offset_in_repeat", 
  "start_in_repeat_database", "end_in_repeat_database", "ID")

# Extract class and family from taxonomy and assign NA if family is missing
repeatmasker_table$class = gsub("/.*", "", repeatmasker_table$taxonomy)
repeatmasker_table$family = ifelse(grepl("/", repeatmasker_table$taxonomy), 
  gsub(".*/", "", repeatmasker_table$taxonomy), NA)

# Select desired columns from table
repeatmasker_table = select(repeatmasker_table, class, family, name, ID, sw_score, 
  substitution_pc, deletion_pc, insertion_pc, seqnames, start, end)

# Replace "C" with "-" in strand column
repeatmasker_table$strand[repeatmasker_table$strand == "C"] = "-"

# Turn the table into a GRanges object
repeatmasker_granges = makeGRangesFromDataFrame(repeatmasker_table, keep.extra.columns = T)

# Remove ranges not on standard chromosomes. 
seqlevels(repeatmasker_granges, pruning.mode = "coarse") = extractSeqlevels(species = "Homo_sapiens", style = "UCSC")

# Save repeatmasker_granges
saveRDS(repeatmasker_granges, "repeatmasker_granges_ucsc.rds")