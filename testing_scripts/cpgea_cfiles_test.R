library(dplyr)
library(GenomicRanges)
library(SummarizedExperiment)

n178 = data.table::fread("~/wgbs/cpgea/processed_files/N178_WGBS_mC_Identification.stat.gz")
names(n178) = c("seqnames", "start", "strand", "meth_reads", "unmeth_reads", "context", "p-val", "q-val") 
n178 = filter(n178, context == "CG")

n178_gr = makeGRangesFromDataFrame(n178, end.field = "start", keep.extra.columns = T)
n178_gr$total = n178_gr$meth_reads + n178_gr$unmeth_reads
seqlevelsStyle(n178_gr) = "UCSC"
cpg_sites = readRDS("~/genomes/human/cpg_sites/all_cpg_sites_hg38.rds")
cpg_sites$count = 0
cpg_sites_forward = cpg_sites
strand(cpg_sites_forward) = "+"
cpg_sites_reverse = shift(cpg_sites_forward, 1)
strand(cpg_sites_reverse) = "-"

forward_hits = data.frame(findOverlaps(cpg_sites_forward, n178_gr))
cpg_sites_forward$count[forward_hits$queryHits] = n178_gr$total[forward_hits$subjectHits]

reverse_hits = data.frame(findOverlaps(cpg_sites_reverse, n178_gr))
cpg_sites_reverse$count[reverse_hits$queryHits] = n178_gr$total[reverse_hits$subjectHits]

cpg_sites$count = cpg_sites_forward$count + cpg_sites_reverse$count

covered_cpg_sites = cpg_sites[cpg_sites$count > 0]

n178_gr = subsetByOverlaps(n178_gr, cpg_sites)

cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")[, "N178"]
cpgea_meth_rse = subsetByOverlaps(cpgea_meth_rse, covered_cpg_sites)
head(assay(cpgea_meth_rse, 2))
