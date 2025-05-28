# Create a GRanges with annotation of regulatory features and CpG islands

# Load required packages
library(GenomicRanges)
library(biomaRt)
library(doParallel)
library(annotatr)

### Download regulatory feature coordinates from biomart

# Get regulatory feature mart for hg38 for Ensembl version 104
regulatory_mart_hg38 = useEnsembl("ENSEMBL_MART_FUNCGEN", dataset = "hsapiens_regulatory_feature")

# Get names of all regulatory features, except Epigenetically Modified Accessible Regions (EMARs)
regulatory_feature_names_hg38 = listFilterOptions(regulatory_mart_hg38, "regulatory_feature_type_name")
regulatory_feature_names_hg38 = setdiff(regulatory_feature_names_hg38, "EMAR")

# Get Genomic cooridnates for each type of regulatory feature
biomart_regulatory_feature_coordinates_list_hg38 = setNames(foreach(feature = regulatory_feature_names_hg38, .packages = "biomaRt") %do% {
  getBM(mart = regulatory_mart_hg38, filters = c("chromosome_name", "regulatory_feature_type_name"), values = list(1:22, feature), 
        attributes = c("chromosome_start", "chromosome_end", "chromosome_name", "feature_type_description"))
}, nm = regulatory_feature_names_hg38)

# Add "chr" preceding each chromosome name
for (feature in names(biomart_regulatory_feature_coordinates_list_hg38)){
  biomart_regulatory_feature_coordinates_list_hg38[[feature]]$chromosome_name = paste0("chr", biomart_regulatory_feature_coordinates_list_hg38[[feature]]$chromosome_name)
}

# Get sorted GRanges for each feature type 
biomart_regulatory_feature_coordinates_granges_list_hg38 = setNames(foreach(feature = regulatory_feature_names_hg38, .packages = "GenomicRanges") %do% {
  sort(GRanges(biomart_regulatory_feature_coordinates_list_hg38[[feature]]))
}, nm = regulatory_feature_names_hg38)

# Convert biomart_regulatory_feature_coordinates_granges_list_hg38 into a GRanges
all_biomart_features_granges_hg38 = unlist(GRangesList(biomart_regulatory_feature_coordinates_granges_list_hg38))

# Rename feature_type_description as region_type
mcols(all_biomart_features_granges_hg38) = data.frame(region_type = all_biomart_features_granges_hg38$feature_type_description)

# Capitalize regulatory feature names and remove word region from them and rename CTCF Transcription Factor CTCF BS
all_biomart_features_granges_hg38$region_type = stringr::str_to_title(gsub(" region", "", all_biomart_features_granges_hg38$region_type))
all_biomart_features_granges_hg38$region_type[all_biomart_features_granges_hg38$region_type == "Ctcf Transcription Factor"] = "CTCF BS"

### Download CpG islands from UCSC
cpg_island_annotation = build_annotations(genome = "hg38", annotations = "hg38_cpgs")

# Remove hg38 from type and change factor names
cpg_island_annotation$type = factor(gsub("hg38_", "", cpg_island_annotation$type))
cpg_island_annotation$type = dplyr::recode_factor(cpg_island_annotation$type, 
  "cpg_inter" = "Open Sea",
  "cpg_islands" = "CpG Islands",
  "cpg_shelves" = "CpG Shelves",
  "cpg_shores" = "CpG Shores")

# Rename type as region type and remove other metadata columns
mcols(cpg_island_annotation) = data.frame(region_type = cpg_island_annotation$type)

# Load list with locations of exons and introns and convert into a GRanges
pc_transcripts_exons_and_introns_grl = readRDS("pc_transcripts_exons_and_introns_grl.rds")
pc_transcripts_exons_and_introns_gr = unlist(GRangesList(pc_transcripts_exons_and_introns_grl))

# Rename type as region type and remove other metadata columns
mcols(pc_transcripts_exons_and_introns_gr) = data.frame(region_type = mcols(pc_transcripts_exons_and_introns_gr)$region)
pc_transcripts_exons_and_introns_gr$region_type = stringr::str_to_title(pc_transcripts_exons_and_introns_gr$region_type)

# Combine all_biomart_features_granges_hg38 and cpg_island_annotation and save
complete_regulatory_annotation = unname(sort(c(pc_transcripts_exons_and_introns_gr, cpg_island_annotation, all_biomart_features_granges_hg38)))
saveRDS(complete_regulatory_annotation, "complete_regulatory_annotation.rds")