# Create a GRanges object with regions where WGBS reads cannot be mapped with confidence

# Load required packages
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Download Bismap 100 bedGraph and read it in
system("wget https://bismap.hoffmanlab.org/raw/hg38/k100.bismap.bedgraph.gz")
k100_bismap = data.table::fread("k100.bismap.bedgraph.gz", skip = 1, col.names = c("seqnames", "start", "end", "score"))

# Convert k100_bismap to a GRanges and save
system.time({k100_bismap_gr = makeGRangesFromDataFrame(k100_bismap, starts.in.df.are.0based = T, keep.extra.columns = T)})
saveRDS(k100_bismap_gr, "k100_bismap_gr.rds")
k100_bismap_gr = readRDS("k100_bismap_gr.rds")

# Get regions with a mappability score of 1
high_mappability_regions = k100_bismap_gr[k100_bismap_gr$score == 1]

# Create a GRanges for hg38
hg38_gr = GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:24]

# Extract regions which do not overlap high mappable regions and save
low_mappability_regions = setdiff(hg38_gr, high_mappability_regions)
saveRDS(low_mappability_regions, "low_mappability_regions.rds")
