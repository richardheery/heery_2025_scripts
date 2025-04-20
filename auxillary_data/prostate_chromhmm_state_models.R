# Create GRanges objects for the chromHMM models for prostate cancer (downloaded from https://ngdc.cncb.ac.cn/omix/release/OMIX237)
# from the paper "Integrative Epigenome Map of the Normal Human Prostate Provides Insights Into Prostate Cancer Predisposition"

# Description of the states can be found here: https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html
# 15 state is based on 5 core histone marks (H3K27me3, H3K4me3, H3K4me1, H3K9me3 H3K36me3). 
# 18 state is based on the 5 core histone marks and H3K27ac. 

# Load requireed packages
library(rtracklayer)
library(genomeTools)

# Download data
system("wget https://download.cncb.ac.cn/OMIX/OMIX237/OMIX237-64-02.zip")
system("unzip OMIX237-64-02.zip")

# Get liftover chain
hg19_to_hg38_chain = rtracklayer::import.chain("../genomes/hg19ToHg38.over.chain")

# Create GRanges objects for the 18 model, keep only standard chromosomes and sort. Has 623,664 ranges. No ranges overlap. 
prostate_18_states = import.bed("prostate_18states_joint_model_dense.bed")
seqlevels(prostate_18_states, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)[1:25]
prostate_18_states = sort(prostate_18_states)

# Convert into a GRanges object and remove score, itemRgb and thick columns
prostate_18_states = GRanges(prostate_18_states)
mcols(prostate_18_states)[c("score", "itemRgb", "thick")] = NULL

# Liftover from hg19 to hg38, 620,543 ranges which do not overlap. Median width 1,000 bp.
system.time({prostate_18_states_hg38 = genomeTools::liftover_granges(genomic_regions = prostate_18_states, chain = hg19_to_hg38_chain)})
prostate_18_states_hg38$score = NULL

# Convert chromatin state name to a factor and put in the desired order
prostate_18_states_hg38$name = factor(prostate_18_states_hg38$name, 
  levels = gtools::mixedsort(unique(prostate_18_states_hg38$name)))

# Get chromatin state description 
chrom_state_18_desc = read.table("chromatin_states_18_description.txt", sep = "\t", header = T)

# Add full chromatin state description to prostate_18_states_hg38
prostate_18_states_hg38$description = chrom_state_18_desc$DESCRIPTION[as.numeric(prostate_18_states_hg38$name)]
prostate_18_states_hg38$description = factor(prostate_18_states_hg38$description, chrom_state_18_desc$DESCRIPTION)

# Save prostate_18_states_hg38
saveRDS(prostate_18_states_hg38, "prostate_18_states_hg38_gr.rds")
