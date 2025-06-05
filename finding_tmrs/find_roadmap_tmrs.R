# Find TMRs in Roadmap samples

# Load required packages
library(methodical)
library(plotR)
library(patchwork)

# Get CAGE-supported MANE transcripts
mane_transcripts = readRDS("~/mounts/local_mount/genomes/genome_annotation_new/transcripts/gencode_v38/mane_pc_transcript_ids.rds")
cage_supported_transcripts = readRDS("~/mounts/local_mount/genomes/fantom/cage_supported_gencode_tss.rds")
cage_supported_mane_transcripts = intersect(mane_transcripts, cage_supported_transcripts)

# Get TSS Granges for TSS and transcripts and subset for cage_supported_mane_transcripts
tss_gr = readRDS("~/mounts/local_mount/genomes/genome_annotation_new/transcripts/gencode_v38/pc_transcripts_tss_gr.rds")[cage_supported_mane_transcripts]
transcripts_gr = readRDS("~/mounts/local_mount/genomes/genome_annotation_new/transcripts/gencode_v38/pc_transcripts_gr.rds")[cage_supported_mane_transcripts]

# Rename using gene ID
names(tss_gr) = gsub("\\.[0-9]*", "", tss_gr$gene_id)
names(transcripts_gr) = gsub("\\.[0-9]*", "", transcripts_gr$gene_id)

# Expand transcripts_gr
transcripts_gr = methodical:::expand_granges(transcripts_gr, 5000, 5000)

# Load Roadmap meth RSE and gene counts
roadmap_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("roadmap_meth_rse_hg38/")
gene_counts = data.frame(data.table::fread("roadmap_gene_estimated_gene_counts.tsv.gz"), row.names = 1)

# Find common samples
common_samples = intersect(names(gene_counts), colnames(roadmap_meth_rse))

# Create a bpparm object
bpparam = BiocParallel::MulticoreParam(workers = 3) 

# Calculate methylation-transcription correlations. Took 22 minutes with 10 cores.  
system.time({transcript_meth_cors_roadmap = calculateMethSiteTranscriptCors(meth_rse = roadmap_meth_rse , 
  transcript_expression_table = gene_counts, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", samples_subset = common_samples, BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_roadmap, "transcript_meth_cors_roadmap.rds")
transcript_meth_cors_roadmap = readRDS("transcript_meth_cors_roadmap.rds")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(20)

# Find TMRs for correlations. Took 1 minute. 
system.time({roadmap_tmrs = findTMRs(correlation_list = transcript_meth_cors_roadmap, 
  p_adjust_method = "fdr", p_value_threshold = 0.1, BPPARAM = bpparam)})
saveRDS(roadmap_tmrs, "roadmap_tmrs.rds")
roadmap_tmrs = readRDS("roadmap_tmrs.rds")

# Get repeat ranges for hg38 and subset TMRs for those not overlapping repeats
# There are 478 TMRs (289 negative and 189 positive)
repeat_ranges = readRDS("~/genomes/repetitive_sequences/repeatmasker/repeatmasker_granges_ucsc.rds")
roadmap_tmrs_no_repeats = subsetByOverlaps(roadmap_tmrs, repeat_ranges, invert = T)
saveRDS(roadmap_tmrs_no_repeats, "roadmap_tmrs_no_repeats.rds")

### Create TMR distribution plots

# Load Roadmap TMRs without repeats
roadmap_tmrs = readRDS("roadmap_tmrs_no_repeats.rds")

# Define a function for ploting overlap of TMRs with transcript regions
plot_tmr_regions = function(tmrs, transcript_regions_gr, regions_filter, normalize = T, title = "Distribution of TMRs"){
  
  # Get title size of each class of region
  region_sizes = sapply(split(transcript_regions_gr, transcript_regions_gr$region), function(x)
    sum(width(reduce(x, ignore.strand = T))))
  
  # Find overlaps between tmrs and transcript_regions_gr
  overlaps_df = data.frame(findOverlaps(tmrs, transcript_regions_gr))
  
  # Add names of transcripts, direction and region to overlap results
  overlaps_df$transcript1 = tmrs$ID[overlaps_df$queryHits]
  overlaps_df$direction = tmrs$direction[overlaps_df$queryHits]
  overlaps_df$transcript2 = transcript_regions_gr$transcript_id[overlaps_df$subjectHits]
  overlaps_df$region = transcript_regions_gr$region[overlaps_df$subjectHits]
  
  # Filter for overlaps for TMRs and their associated transcript
  overlaps_df = dplyr::filter(overlaps_df, transcript1 == transcript2)
  
  # Count number of TMRs per region
  overlaps_summary = summarize(group_by(overlaps_df, direction, region), count = n())
  overlaps_summary = tidyr::complete(tibble(overlaps_summary), direction, region, fill = list(count = 0))
  
  # Normalize count by coverage of region
  overlaps_summary$normalized_count = overlaps_summary$count/region_sizes[overlaps_summary$region]*1e6
  overlaps_summary = data.frame(overlaps_summary)
  
  # Filter for desired regions
  tmr_overlaps_count_filtered = dplyr::filter(overlaps_summary, region %in% regions_filter)
  regions_filter = gsub("intron_", "Intron ", regions_filter)
  regions_filter = gsub("exon_", "Exon ", regions_filter)
  regions_filter = gsub("500", "250", regions_filter)
  regions_filter = gsub("1000", "750", regions_filter)
  regions_filter = gsub("1500", "1250", regions_filter)
  regions_filter = gsub("2000", "1750", regions_filter)
  
  # Convert region into a factor
  tmr_overlaps_count_filtered$region = gsub("intron_", "Intron ", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("exon_", "Exon ", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("500", "250", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("1000", "750", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("1500", "1250", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = gsub("2000", "1750", tmr_overlaps_count_filtered$region)
  tmr_overlaps_count_filtered$region = factor(tmr_overlaps_count_filtered$region, levels = regions_filter)
  
  # Update count to normalized count if normalize is TRUE
  ylab = "TMR Count"
  if(normalize){
    tmr_overlaps_count_filtered$count = tmr_overlaps_count_filtered$normalized_count
    ylab = "Number of TMRs per MB"
  }
  
  # Add a column indicating if regions are upstream, in the transcribed region or downstream
  tmr_overlaps_count_filtered$region_class = "Transcribed Region"
  tmr_overlaps_count_filtered$region_class[grep("TSS", tmr_overlaps_count_filtered$region)] = "Upstream of TSS"
  tmr_overlaps_count_filtered$region_class[grep("TES", tmr_overlaps_count_filtered$region)] = "Downstream of TES"
  tmr_overlaps_count_filtered$region_class = factor(tmr_overlaps_count_filtered$region_class, 
    levels = c("Upstream of TSS", "Transcribed Region", "Downstream of TES"))
  
  # Create plot of number of TMRs per region
  regions_plot_normalized = ggplot(tmr_overlaps_count_filtered, aes(x = region, y = count, fill = direction)) +
    geom_col(position = "dodge") + scale_fill_manual(values = c(plotR::colour_list$purple_and_gold_light), drop = F)
  regions_plot_normalized = plotR::customize_ggplot_theme(regions_plot_normalized, 
    title = NULL, 
    xlab = NULL, ylab = ylab,
    fill_colors = c(plotR::colour_list$purple_and_gold_light), x_labels_angle = 45) +
    scale_fill_manual(values = c(plotR::colour_list$purple_and_gold_light), drop = F)
  regions_plot_normalized + facet_grid(~region_class, scales = "free_x", space = "free_x") + 
    theme(panel.spacing = unit(0,'lines'), strip.background = element_blank(),
      strip.text = element_text(size = 16))
  
}

source("~/promoter_project_gene_body/final_tmrs/expand_transcripts.R")

# Load introns and exons for transcripts as a GRanges list
tss_grl = readRDS("~/genomes/genome_annotation_new/transcripts/gencode_v38/pc_transcripts_exons_and_introns_grl.rds")
cage_transcripts = readRDS("~/genomes/fantom/cage_supported_gencode_tss.rds")
mane_transcripts = readRDS("~/genomes/genome_annotation_new/transcripts/gencode_v38/mane_pc_transcript_ids.rds")
tss_grl = tss_grl[cage_transcripts]
tss_grl_expanded = expand_transcripts(tss_grl, seq(500, 2000, 500), seq(500, 2000, 500))

exons = paste0("exon_", 1:10)
introns = paste0("intron_", 1:10)
exons_introns = c(rbind(exons, introns))
exons_introns = exons_introns[-length(exons_introns)]
regions = c(paste0("TSS-", rev(seq(500, 2000, 500))), exons_introns, paste0("TES+", seq(500, 2000, 500)))

roadmap_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = F)
roadmap_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = roadmap_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = T)
combined_roadmap_tmrs_introns_exons_plot = ggpubr::ggarrange(roadmap_tmrs_introns_exons_plot, roadmap_tmrs_introns_exons_plot_normalized, 
  align = "hv", nrow = 2, labels = c("A", "B"), common.legend = T, legend = "right")
ggsave(plot = combined_roadmap_tmrs_introns_exons_plot, "~/promoter_project_gene_body/final_plots/supp_figure_roadmap_tmr_distribution.pdf", width = 16, height = 9)

# Make regulatory feature and chromatin state plot

# Get a GRanges for transcripts and expand
cage_tss = readRDS("~/genomes/fantom/cage_supported_gencode_tss.rds")
transcripts_gr = readRDS("~/genomes/genome_annotation_new/transcripts/gencode_v38/pc_transcripts_gr.rds")
transcripts_gr = genomicTools::adjust_gr(transcripts_gr, 5000, 5000)
transcripts_gr = transcripts_gr[cage_tss]

# Get CpG islands from UCSC
cpg_island_annotation <- annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs")
cpg_island_annotation$type = factor(gsub("hg38_", "", cpg_island_annotation$type))
cpg_island_annotation$type = recode_factor(cpg_island_annotation$type, 
  "cpg_inter" = "Open Sea",
  "cpg_islands" = "CpG Islands",
  "cpg_shelves" = "CpG Shelves",
  "cpg_shores" = "CpG Shores")
cpg_island_annotation$region_type = cpg_island_annotation$type
cpg_island_annotation$type = NULL

# Load hg38 genome annotation
genome_annotation_hg38 = methodicalFinal::genome_annotation_hg38

# Remove CpG Islands
genome_annotation_hg38 = filter_gr(genome_annotation_hg38, "region_type", values = "CpG Island", invert = T)

# Shorten names of exons and introns
genome_annotation_hg38$region_type[genome_annotation_hg38$region_type == "exon"] = "Exons"
genome_annotation_hg38$region_type[genome_annotation_hg38$region_type == "intron"] = "Introns"

# Remove " Region" from cluster_annotation$region_type
genome_annotation_hg38$region_type = gsub(" Region", "", genome_annotation_hg38$region_type)

# Remove repeats from annotation
genome_annotation_hg38 = plyranges::filter(genome_annotation_hg38, !region_type %in% 
    c("Alu", "CR1", "DNA Transposon", "L1", "L2", "LTR", "Low Complexity", "MIR", "SVA", "Satellite", "Simple Repeat", "Protein-Coding"))

# Remove lncRNA regions
genome_annotation_hg38 = plyranges::filter(genome_annotation_hg38, !grepl("lncRNA", region_type))

# Combine genome_annotation_hg38 and cpg_island_annotation
genome_annotation_hg38 = c(genome_annotation_hg38, cpg_island_annotation)

# Add transcripts_gr as background TMR search space
background_gr = transcripts_gr
mcols(background_gr) = NULL
background_gr$region_type = "Background"
genome_annotation_hg38 = c(genome_annotation_hg38, background_gr)

# Convert genome_annotation_hg38 to a list
genome_annotation_hg38_list = GRangesList(split(genome_annotation_hg38, genome_annotation_hg38$region_type))

# Calculate the size of the overlaps of transcripts_gr with each chromatin state
annotation_widths = sapply(genome_annotation_hg38_list, function(x) calculate_regions_intersections(gr1 = transcripts_gr, gr2 = x))

# Find overlaps between TMRs and chromatin states
tmr_annotation_overlaps = data.frame(findOverlaps(roadmap_tmrs, genome_annotation_hg38))
tmr_annotation_overlaps$tmr_group = "Roadmap"
tmr_annotation_overlaps$direction = roadmap_tmrs$direction[tmr_annotation_overlaps$queryHits]
tmr_annotation_overlaps$annotation = genome_annotation_hg38$region_type[tmr_annotation_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_annotation_overlaps_summary = dplyr::summarise(group_by(tmr_annotation_overlaps, tmr_group, direction, annotation), count = n())

# Normalize count by MB covered by the chromatin state
tmr_annotation_overlaps_summary$normalized_count = tmr_annotation_overlaps_summary$count/annotation_widths[tmr_annotation_overlaps_summary$annotation]*1e6

# Convert region_type to a factor and give specified order
tmr_annotation_overlaps_summary = filter(tmr_annotation_overlaps_summary, annotation %in% 
    c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))
tmr_annotation_overlaps_summary$annotation = factor(tmr_annotation_overlaps_summary$annotation, 
  c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))
tmr_annotation_overlaps_summary = data.frame(complete(tibble(tmr_annotation_overlaps_summary), tmr_group, annotation, direction, fill = list(normalized_count = 0, count = 0)))

# Create a plot annotating TMRs and save
cluster_genomic_feature_annotation_plot = ggplot(tmr_annotation_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black") 
cluster_genomic_feature_annotation_plot = customize_ggplot_theme(cluster_genomic_feature_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "annotation", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 14, 
  legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm"), axis.text.x = element_blank()) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
cluster_genomic_feature_annotation_plot

### Annotate chromatin states

# Get a GRanges with chromatin states for prostate 
prostate_18_states_hg38_gr = readRDS("~/roadmap/prostate_roadmap/prostate_18_states_hg38_gr.rds")

# Shorten "PolyComb" in description to "PC" 
levels(prostate_18_states_hg38_gr$description)[3] = "Flanking TSS 5'"
levels(prostate_18_states_hg38_gr$description)[4] = "Flanking TSS 3'"
levels(prostate_18_states_hg38_gr$description)[5] = "Strong Transcription"
levels(prostate_18_states_hg38_gr$description)[6] = "Weak Transcription"
levels(prostate_18_states_hg38_gr$description)[7] = "Genic Enhancer 1"
levels(prostate_18_states_hg38_gr$description)[8] = "Genic Enhancer 2"
levels(prostate_18_states_hg38_gr$description)[12] = "ZNF Genes/Repeats"
levels(prostate_18_states_hg38_gr$description)[16] = "Repressed PC"
levels(prostate_18_states_hg38_gr$description)[17] = "Weak Repressed PC"

# Split prostate_18_states_hg38_gr into a list
prostate_18_states_hg38_gr$region_type = prostate_18_states_hg38_gr$description
background_gr2 = background_gr
background_gr2$region_type = "Background2"
prostate_18_states_hg38_gr = c(background_gr, prostate_18_states_hg38_gr, background_gr2)
prostate_18_states_hg38_gr$region_type = factor(prostate_18_states_hg38_gr$region_type, 
  levels = c(levels(prostate_18_states_hg38_gr$region_type)[1:10], levels(prostate_18_states_hg38_gr$region_type)[20], 
    levels(prostate_18_states_hg38_gr$region_type)[11:19]))
prostate_18_states_hg38_gr_list = split(prostate_18_states_hg38_gr, prostate_18_states_hg38_gr$region_type)

# Calculate the size of the ovelaps of transcripts_gr with each chromatin state
chromatin_state_widths = sapply(prostate_18_states_hg38_gr_list, function(x) calculate_regions_intersections(gr1 = transcripts_gr, gr2 = x))

# Find overlaps between TMRs and chromatin states
tmr_state_overlaps = data.frame(findOverlaps(roadmap_tmrs, prostate_18_states_hg38_gr))
tmr_state_overlaps$tmr_group = "Roadmap"
tmr_state_overlaps$direction = roadmap_tmrs$direction[tmr_state_overlaps$queryHits]
tmr_state_overlaps$chromatin_state = prostate_18_states_hg38_gr$region_type[tmr_state_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_state_overlaps_summary = summarise(group_by(tmr_state_overlaps, tmr_group, direction, chromatin_state), count = n())

# Normalize count by MB covered by the chromatin state
tmr_state_overlaps_summary$normalized_count = tmr_state_overlaps_summary$count/chromatin_state_widths[tmr_state_overlaps_summary$chromatin_state]*1e6

# Create a plot annotating TMRs and save
cluster_chromatin_state_annotation_plot = ggplot(tmr_state_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black")
cluster_chromatin_state_annotation_plot = customize_ggplot_theme(cluster_chromatin_state_annotation_plot, title = NULL, 
  xlab = NULL, ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "chromatin_state", facet_nrow = 2, facet_scales = "free_x", strip_text_size = 17, axis_text_size = 16, axis_title_size = 24, legend_key_size = 1.5, 
  legend_title_size = 24, legend_text_size = 20) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm"), axis.text.x = element_blank()) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50)))) +
  facet_wrap(~chromatin_state, drop = F, nrow = 2)
cluster_chromatin_state_annotation_plot

### Combine The TMR distribution plots and TMR annotation plots 
patchwork_plot = cluster_genomic_feature_annotation_plot / cluster_chromatin_state_annotation_plot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = patchwork_plot, filename = "tmr_plots/tmr_annotation_plots.pdf", width = 27, height = 27)

### Check overlap of roadmap TMRs with other TMRs

# Get list of TMRs
roadmap_tmrs_negative = roadmap_tmrs[roadmap_tmrs$direction == "Negative"]
roadmap_tmrs_positive = roadmap_tmrs[roadmap_tmrs$direction == "Positive"]
tmr_list = readRDS("~/promoter_project_gene_body/final_tmrs/gene_body_tmrs/cage_tmr_list_no_repeats.rds")
tmr_list_commbined = unlist(GRangesList(tmr_list))

sapply(tmr_list, function(x) length(subsetByOverlaps(roadmap_tmrs_negative, x)))
sapply(tmr_list, function(x) length(subsetByOverlaps(roadmap_tmrs_positive, x)))
sapply(tmr_list, function(x) length(subsetByOverlaps(roadmap_tmrs_negative, x)))/length(roadmap_tmrs_negative)
sapply(tmr_list, function(x) length(subsetByOverlaps(roadmap_tmrs_positive, x)))/length(roadmap_tmrs_positive)
length(subsetByOverlaps(roadmap_tmrs, tmr_list_commbined))/length(roadmap_tmrs)

sort(unique(subsetByOverlaps(roadmap_tmrs, tmr_list_commbined)$gene_name))
