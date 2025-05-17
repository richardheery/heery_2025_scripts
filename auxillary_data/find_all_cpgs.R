# Make BED files and GRanges objects with all the C position of all CpG sites in hg38 and hg19.

# Load required packages
library(BSgenome.Hsapiens.UCSC.hg38) 

# Create BSParams object for hg38
params_hg38 = new("BSParams", X = BSgenome.Hsapiens.UCSC.hg38, FUN = matchPattern, exclude = "_")

# Get all CpG sites in hg38 as a list of XStringsViews for standard chromosomes
all_cpgs_hg38  = bsapply(BSParams = params_hg38, pattern = "CG")[extractSeqlevelsByGroup(species = "Homo_sapiens", style = "UCSC", group = "all")]

# Convert all_cpgs_hg38 into a GRanges object for the C bases of each CpG site
cpg_genome_ranges_hg38 = GRanges(seqnames = rep(names(all_cpgs_hg38),  lengths(all_cpgs_hg38)),
  ranges = IRanges(unlist(lapply(all_cpgs_hg38, start)), unlist(lapply(all_cpgs_hg38, start))))

# Remove objects for hg38
rm(all_cpgs_hg38, params_hg38, BSgenome.Hsapiens.UCSC.hg38)
gc() 

# Save cpg_genome_ranges_hg38 as a BED file and a GRanges object
saveRDS(cpg_genome_ranges_hg38, "all_cpg_sites_hg38.rds")
