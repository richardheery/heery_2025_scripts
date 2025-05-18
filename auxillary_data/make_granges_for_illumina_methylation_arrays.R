# Create GRanges objects with the hg38 genomic coordinates of the Illumina array probes

# Load required packages
library(dplyr)
library(GenomicRanges)
library(BSgenome)

### Create a GRanges for the HumanMethylation450 array

# Download manifest file for EPIC methylation array from Illumina
system("wget https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv")

# Load HumanMethylation450 probe information table
infinium_450k_probe_data_hg19 = data.table::fread("humanmethylation450_15017482_v1-2.csv", header = T, skip = 7, sep = ",", nrows = 485585, fill=TRUE)

# Remove any probes that are not for measuring cytosine methylation (SNP probes and controls)
infinium_450k_probe_data_hg19 = filter(infinium_450k_probe_data_hg19, grepl("c", IlmnID))

# Select the infiniumID, chromosome and genomic coordinate (MAPINFO) columns
infinium_450k_probe_data_hg19 = transmute(infinium_450k_probe_data_hg19, IlmnID, 
  CHR = paste0("chr", CHR), MAPINFO, refgene_id = UCSC_RefGene_Accession, refgene_name = UCSC_RefGene_Name)

# Convert into a GRanges object
infinium_450k_probe_granges_hg19 = makeGRangesFromDataFrame(infinium_450k_probe_data_hg19, 
  seqnames.field = "CHR", start.field = "MAPINFO", end.field = "MAPINFO", keep.extra.columns = T)

# Add probe name to infinium_450k_probe_granges_hg19
infinium_450k_probe_granges_hg19$name = infinium_450k_probe_granges_hg19$IlmnID
infinium_450k_probe_granges_hg19$IlmnID = NULL

# Sort infinium_450k_probe_granges_hg19
infinium_450k_probe_granges_hg19 = sort(infinium_450k_probe_granges_hg19)

# Remove probes which do not overlap a CpG site in hg19
probe_seqs_hg19 = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, resize(infinium_450k_probe_granges_hg19, fix = "start", width = 2)))
infinium_450k_probe_granges_hg19 = infinium_450k_probe_granges_hg19[probe_seqs_hg19 == "CG"]

# Load liftover chain
hg19_to_hg38_liftover_chain = rtracklayer::import.chain("hg19ToHg38.over.chain")

# Liftover infinium_450k_probe_grange to hg38
infinium_450k_probe_granges_hg38 = genomicTools::liftover_granges(infinium_450k_probe_granges_hg19, chain = hg19_to_hg38_liftover_chain)

# Remove probes which do not overlap a CpG site in hg38. Leaves 480,975 probes
probe_seqs_hg38 = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, resize(infinium_450k_probe_granges_hg38, fix = "start", width = 2)))
infinium_450k_probe_granges_hg38 = infinium_450k_probe_granges_hg38[probe_seqs_hg38 == "CG"]

# Put seqlevels in correct order
seqlevels(infinium_450k_probe_granges_hg38) = seqlevels(infinium_450k_probe_granges_hg19)

# Save probe ranges
saveRDS(infinium_450k_probe_granges_hg38, "infinium_450k_probe_granges_hg38.rds")

### Create a GRanges for the EPIC array

# Download manifest file for EPIC methylation array from Illumina
system(" wget https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")
system("unzip infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip")

# Read in the probe annotation file for EPIC array. 
# NB that MAPINFO refers to the position of the C of the CpG on the forward/+ strand.
# Start_hg38 and End_hg38 are 0-based half-open coordinates (BED format) for the CpG 
system.time({epic_probe_table = read.csv("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)})

# Add a column with the type of probe by extracting the first two letters of the probe name
epic_probe_table$probe_type = substr(epic_probe_table$IlmnID, 1, 2)

# Filter table to remove probes which are for detecting SNPs (rs) or bind to non-CpG sites
epic_probe_table = filter(epic_probe_table, !probe_type %in% c("ch", "rs"))

# Select desired columns 
epic_probe_table = select(epic_probe_table, Name, CHR, Strand, MAPINFO, CHR_hg38, Strand_hg38, Start_hg38, End_hg38, GencodeCompV12_Accession, Regulatory_Feature_Group)

# Rename chromsomes so that they are in UCSC format
epic_probe_table$CHR = paste0("chr", epic_probe_table$CHR)

# Turn probe coordinate table into GRanges object for hg38 and save
epic_probe_gr_hg38 = with(filter(epic_probe_table, !is.na(Start_hg38) & !is.na(CHR_hg38)), sort(GRanges(seqnames = CHR_hg38, ranges = IRanges(start = Start_hg38 + 1),
  ensembl_transcript_id = GencodeCompV12_Accession, regulatory_feature = Regulatory_Feature_Group, seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38), names = Name), ignore.strand = T))
saveRDS(epic_probe_gr_hg38, "epic_probe_gr_hg38.rds")

