
# Load required packages
library(dplyr)
library(GenomicRanges)
library(SummarizedExperiment)

# Read in N178 WGBS file and filter for CpG sites
n178 = data.table::fread("~/wgbs/cpgea/processed_files/N178_WGBS_mC_Identification.stat.gz")
names(n178) = c("seqnames", "start", "strand", "meth_reads", "unmeth_reads", "context", "p-val", "q-val") 
n178 = filter(n178, context == "CG")

# Create a Granges object from n178
n178_gr = makeGRangesFromDataFrame(n178, end.field = "start", keep.extra.columns = T)

# Add a column for coverage and meth_proportion
n178_gr$coverage = n178_gr$meth_reads + n178_gr$unmeth_reads

# Change seqlevels style to UCSC
seqlevelsStyle(n178_gr) = "UCSC"

# Read in a GRanges with all CpG sites for hg38
cpg_sites = readRDS("~/genomes/human/cpg_sites/all_cpg_sites_hg38.rds")

# Add columns for coverage and number of methylated reads
cpg_sites$coverage = 0
cpg_sites$meth_reads = 0

# Create separate GRanges for forward and reverse strand
cpg_sites_forward = cpg_sites
strand(cpg_sites_forward) = "+"
cpg_sites_reverse = shift(cpg_sites_forward, 1)
strand(cpg_sites_reverse) = "-"

# Find hits between forward CpGs and n178_gr
forward_hits = data.frame(findOverlaps(cpg_sites_forward, n178_gr))

# Add coverage and meth_reads to forward CpGs
cpg_sites_forward$coverage[forward_hits$queryHits] = n178_gr$coverage[forward_hits$subjectHits]
cpg_sites_forward$meth_reads[forward_hits$queryHits] = n178_gr$meth_reads[forward_hits$subjectHits]

# Find hits between reverse CpGs and n178_gr
reverse_hits = data.frame(findOverlaps(cpg_sites_reverse, n178_gr))

# Add coverage and meth_reads to reverse CpGs
cpg_sites_reverse$coverage[reverse_hits$queryHits] = n178_gr$coverage[reverse_hits$subjectHits]
cpg_sites_reverse$meth_reads[reverse_hits$queryHits] = n178_gr$meth_reads[reverse_hits$subjectHits]

# Combine the coverage and meth_reads from both strands
cpg_sites$coverage = cpg_sites_forward$coverage + cpg_sites_reverse$coverage
cpg_sites$meth_reads = cpg_sites_forward$meth_reads + cpg_sites_reverse$meth_reads

# Add proportion of methylated reads
cpg_sites$meth_proportion = cpg_sites$meth_reads/cpg_sites$coverage

# Get covered CpG sites
# covered_cpg_sites = cpg_sites[cpg_sites$coverage > 0]

# Load CPGEA meth RSE and subset for cpg_sites
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")[, "N178"]

# Change NA coverage values to 0
assay(cpgea_meth_rse, 2)[is.na(assay(cpgea_meth_rse, 2))] = 0

# Check if calculated coverage and methylated reads are the same as in cpgea_meth_rse. They are
system.time({cpgea_meth_rse_meth_proportion = as.matrix(assay(cpgea_meth_rse, 1))})
system.time({cpgea_meth_rse_meth_coverage = as.matrix(assay(cpgea_meth_rse, 2))})
fivenum(as.numeric(cpgea_meth_rse_meth_proportion) - cpg_sites$meth_proportion)
fivenum(as.numeric(cpgea_meth_rse_meth_coverage) - cpg_sites$coverage)

# Load genome annotation and split into a list
genome_annotation = readRDS("complete_genome_annotation.rds")
genome_annotation_list = split(genome_annotation, genome_annotation$region_type)

# Calculate mean coverage and methylation for different classes of genomic regions in the N178 WGBS file
region_mean_methylation = sapply(genome_annotation_list, function(x) mean(subsetByOverlaps(cpg_sites, x)$meth_proportion, na.rm = T))
region_mean_coverage = sapply(genome_annotation_list, function(x) mean(subsetByOverlaps(cpg_sites, x)$coverage, na.rm = T))
saveRDS(region_mean_methylation, "N178_region_mean_methylation.rds")
saveRDS(region_mean_coverage, "N178_region_mean_coverage.rds")
cpgea_region_mean_coverage = readRDS("N178_region_mean_coverage.rds")

# Calculate mean coverage and methylation for different classes of genomic regions in cpgea_meth_rse
system.time({region_mean_methylation_meth_rse = sapply(genome_annotation_list, function(x) colMeans(assay(subsetByOverlaps(cpgea_meth_rse, x), 1), na.rm = T))})
system.time({region_mean_coverage_meth_rse = sapply(genome_annotation_list, function(x) colMeans(assay(subsetByOverlaps(cpgea_meth_rse, x), 2)))})
saveRDS(region_mean_methylation_meth_rse, "N178_region_mean_methylation_meth_rse.rds")
saveRDS(region_mean_coverage_meth_rse, "N178_region_mean_coverage_meth_rse.rds")
fivenum(region_mean_methylation - region_mean_methylation_meth_rse)
fivenum(region_mean_coverage - region_mean_coverage_meth_rse)

# Get forward strand positions for n178_gr
n178_gr_forward = n178_gr[strand(n178_gr) == "+"]

# Count how many of these CpGs are in n178_gr_forward
region_overlaps_cpgea = sapply(genome_annotation_list, function(x) 
  length(subsetByOverlaps(n178_gr_forward, x))/length(subsetByOverlaps(cpg_sites, x)))
saveRDS(region_overlaps_cpgea, "region_overlaps_cpgea.rds")

# Load 
mcrpc = data.table::fread("NT-DTB-077-Pro.merged.CpG_report.txt.gz")
names(mcrpc) = c("seqnames", "start", "strand", "meth_reads", "unmeth_reads", "context1", "context2") 
mcrpc = dplyr::filter(mcrpc, context1 == "CG")

# Create a Granges object from mcrpc
mcrpc_gr = makeGRangesFromDataFrame(mcrpc, end.field = "start", keep.extra.columns = T)

# Create a table for mapping refseq chromosome names to UCSC
refseq_to_ucsc <- c(
  "NC_000001.11" = "chr1",
  "NC_000002.12" = "chr2",
  "NC_000003.12" = "chr3",
  "NC_000004.12" = "chr4",
  "NC_000005.10" = "chr5",
  "NC_000006.12" = "chr6",
  "NC_000007.14" = "chr7",
  "NC_000008.11" = "chr8",
  "NC_000009.12" = "chr9",
  "NC_000010.11" = "chr10",
  "NC_000011.10" = "chr11",
  "NC_000012.12" = "chr12",
  "NC_000013.11" = "chr13",
  "NC_000014.9"  = "chr14",
  "NC_000015.10" = "chr15",
  "NC_000016.10" = "chr16",
  "NC_000017.11" = "chr17",
  "NC_000018.10" = "chr18",
  "NC_000019.10" = "chr19",
  "NC_000020.11" = "chr20",
  "NC_000021.9"  = "chr21",
  "NC_000022.11" = "chr22",
  "NC_000023.11" = "chrX",
  "NC_000024.10" = "chrY",
  "NC_012920.1"  = "chrM"
)

# Put seqlevels in UCSC format
seqlevels(mcrpc_gr, pruning.mode = "coarse") = names(refseq_to_ucsc)
seqlevels(mcrpc_gr) = refseq_to_ucsc[seqlevels(mcrpc_gr)]

# Add a column for coverage and meth_proportion
mcrpc_gr$coverage = mcrpc_gr$meth_reads + mcrpc_gr$unmeth_reads

# Filter mcrpc_gr for covered CpGs
mcrpc_gr = mcrpc_gr[mcrpc_gr$coverage > 0]

# Get forward strand positions for mcrpc_gr
mcrpc_gr_forward = mcrpc_gr[strand(mcrpc_gr) == "+"]

# Count how many of these CpGs are in mcrpc_gr_forward
region_overlaps_mcrpc = sapply(genome_annotation_list, function(x) 
  length(subsetByOverlaps(mcrpc_gr_forward, x))/length(subsetByOverlaps(cpg_sites, x)))
region_overlaps_mcrpc
saveRDS(region_overlaps_mcrpc, "region_overlaps_mcrpc.rds")
region_overlaps_mcrpc = readRDS("region_overlaps_mcrpc.rds")

cpg_density = sapply(genome_annotation_list, function(x) 
  length(subsetByOverlaps(cpg_sites, x))/sum(width(reduce(x, ignore.strand = T)))*1e6)
cor(cpg_density, region_overlaps_mcrpc)
cor(cpg_density, region_overlaps_cpgea)
cor(cpg_density, cpgea_region_mean_coverage)

