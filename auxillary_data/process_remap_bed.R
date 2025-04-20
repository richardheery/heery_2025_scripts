# Create a GRanges object for the binding sites of different TFs in prostate from the ReMap database

# Load required packages
library(GenomicRanges)
library(dplyr)

# Import BED file for non-redundant ranges from ReMap and convert to a GRanges object. Contains 69 million regions. 
system("wget https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz")
system.time({remap_table = data.table::fread("remap2022_nr_macs2_hg38_v1_0.bed.gz", sep = "\t", header = F, select = 1:6)})

# Add names to table
names(remap_table) = c("seqnames", "start", "end", "name", "score", "strand")
 
# # Add name of TFs and name of cell lines as metadata columns. 1,210 different TFs and 1,412,581 different cell lines found.
system.time({remap_table$cell_line = gsub(".*:", "", remap_table$name)})

# Get names of any cell lines related to prostate. 
prostate_cell_lines = readLines("prostate_cell_lines.txt")

# Add "prostate" to prostate cell lines 
prostate_cell_lines = c(prostate_cell_lines, "prostate")

# Create a string with the names of prostate cell lines or prostate to use with grep
# "PC-3" also matches the pancreatic cancer cell line BxPC-3 and so "[^x]PC-3" is used instead
prostate_cell_lines[prostate_cell_lines == "PC-3"] = "[^x]PC-3"
prostate_grep_pattern = paste(prostate_cell_lines, collapse = "|")

# Filter remap_table for prostate cell lines
system.time({prostate_remap_table = filter(remap_table, grepl(prostate_grep_pattern, cell_line))})

# Add name of transcription factor 
prostate_remap_table$tf = gsub(":.*", "", prostate_remap_table$name)

# Convert prostate_remap_table to a GRanges object and save
prostate_remap_table = makeGRangesFromDataFrame(prostate_remap_table, starts.in.df.are.0based = T, keep.extra.columns = T)
saveRDS(prostate_remap_gr, "prostate_remap_gr.rds")