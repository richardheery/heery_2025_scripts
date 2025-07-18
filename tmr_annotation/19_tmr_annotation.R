# Annotate TMRs using genomic features and chromatin states

# Load required packages
library(GenomicRanges)
library(dplyr)
library(patchwork)
source("../auxillary_scripts/granges_functions.R")
source("../auxillary_scripts/plotting_functions.R")

# Load TMRs filtered for CAGE supported TSS not overlapping repeats
cpgea_normal_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("../finding_tmrs/tmr_granges/mcrpc_tmrs.rds")

# Get list of TMRs
tmr_list = list(
  "Normal Prostate" = cpgea_normal_tmrs,
  "Prostate Tumours" = cpgea_tumour_tmrs,
  "Prostate Metastases" = mcrpc_tmrs)

# Combine tmr_list into a single GRanges
tmrs_combined = unlist(GRangesList(lapply(tmr_list, unname)))
tmrs_combined$group = names(tmrs_combined)

# Get a GRanges for transcripts and expand
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")
transcripts_gr = methodical::expand_granges(transcripts_gr, 5000, 5000)

# Load hg38 genome annotation
genome_annotation_hg38 = readRDS("../auxillary_data/complete_regulatory_annotation.rds")

# Remove introns and exons
#genome_annotation_hg38 = genome_annotation_hg38[!genome_annotation_hg38$region_type %in% c("Exon", "Intron")]

# Add transcripts_gr as background TMR search space
background_gr = transcripts_gr
mcols(background_gr) = NULL
background_gr$region_type = "Background"
#genome_annotation_hg38 = c(genome_annotation_hg38, background_gr)

# Convert genome_annotation_hg38 to a list
genome_annotation_hg38_list = GRangesList(split(genome_annotation_hg38, genome_annotation_hg38$region_type))
genome_annotation_hg38_list = lapply(genome_annotation_hg38_list, function(x) GenomicRanges::intersect(x, transcripts_gr, ignore.strand = T))

######
tmr_list = readRDS("../finding_tmrs/tmr_granges/tmr_list.rds")
true = lapply(tmr_list, function(x) 
    sapply(genome_annotation_hg38_list, function(y)
      sum(width(GenomicRanges::intersect(x, y, ignore.strand = T)))/sum(width(reduce(y, ignore.strand = T)))*1000000))


cl = makeCluster(10)
registerDoParallel(cl, 10)
system.time({permutation_results = foreach(i = 1:100, .combine = rbind, .packages = "GenomicRanges") %dopar% {
  
  print(i)
  
  # 
  genome_annotation_hg38_copy = genome_annotation_hg38
  genome_annotation_hg38_copy$region_type = sample(genome_annotation_hg38_copy$region_type)
  genome_annotation_hg38_copy_list = split(genome_annotation_hg38_copy, genome_annotation_hg38_copy$region_type)
  genome_annotation_hg38_copy_list = lapply(genome_annotation_hg38_copy_list, function(x) GenomicRanges::intersect(x, transcripts_gr, ignore.strand = T))
  
  res = data.frame(lapply(genome_annotation_hg38_copy_list, function(x) 
    sapply(tmr_list, function(y)
      sum(width(GenomicRanges::intersect(x, y, ignore.strand = T)))/sum(width(reduce(x, ignore.strand = T)))*1000000)))
  
  res = tibble::rownames_to_column(res, "tmr_group")
  res
  
}})

permutation_results = split(dplyr::select(permutation_results, -tmr_group), permutation_results$tmr_group)
true$cpgea_normal_tmrs_negative
permutation_results$cpgea_normal_tmrs_negative

(rowSums(apply(permutation_results$cpgea_normal_tmrs_negative, 1, function(x) x > true$cpgea_normal_tmrs_negative))+1)/101

#####

# Calculate the size of the overlaps of transcripts_gr with each chromatin state
annotation_widths = sapply(genome_annotation_hg38_list, function(x) calculate_regions_intersections(gr1 = transcripts_gr, gr2 = x))

# Find overlaps between TMRs and chromatin states
tmr_annotation_overlaps = data.frame(findOverlaps(tmrs_combined, genome_annotation_hg38))
tmr_annotation_overlaps$tmr_group = gsub("_tmrs_.*", "_tmrs", tmrs_combined$group[tmr_annotation_overlaps$queryHits])
tmr_annotation_overlaps$direction = tmrs_combined$direction[tmr_annotation_overlaps$queryHits]
tmr_annotation_overlaps$annotation = genome_annotation_hg38$region_type[tmr_annotation_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_annotation_overlaps_summary = dplyr::summarise(group_by(tmr_annotation_overlaps, tmr_group, direction, annotation), count = n())

# Normalize count by MB covered by the chromatin state
tmr_annotation_overlaps_summary$normalized_count = tmr_annotation_overlaps_summary$count/annotation_widths[tmr_annotation_overlaps_summary$annotation]*1e6

# Convert region_type to a factor and give specified order
tmr_annotation_overlaps_summary$annotation = factor(tmr_annotation_overlaps_summary$annotation, 
  c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))

# Put tmr_groups in right order
tmr_annotation_overlaps_summary$tmr_group = factor(tmr_annotation_overlaps_summary$tmr_group, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Create a plot annotating TMRs and save
tmr_genomic_feature_annotation_plot = ggplot(tmr_annotation_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black") 
tmr_genomic_feature_annotation_plot = customize_ggplot_theme(tmr_genomic_feature_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "annotation", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 14, 
  legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
tmr_genomic_feature_annotation_plot

### Annotate chromatin states

# Get a GRanges with chromatin states for prostate 
prostate_18_states_hg38_gr = readRDS("../auxillary_data/prostate_18_states_hg38_gr.rds")

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
tmr_state_overlaps = data.frame(findOverlaps(tmrs_combined, prostate_18_states_hg38_gr))
tmr_state_overlaps$tmr_group = gsub("_tmrs_.*", "_tmrs", tmrs_combined$group[tmr_state_overlaps$queryHits])
tmr_state_overlaps$direction = tmrs_combined$direction[tmr_state_overlaps$queryHits]
tmr_state_overlaps$chromatin_state = prostate_18_states_hg38_gr$region_type[tmr_state_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_state_overlaps_summary = summarise(group_by(tmr_state_overlaps, tmr_group, direction, chromatin_state), count = n())

# Normalize count by MB covered by the chromatin state
tmr_state_overlaps_summary$normalized_count = tmr_state_overlaps_summary$count/chromatin_state_widths[tmr_state_overlaps_summary$chromatin_state]*1e6

# Put tmr_group in right order
tmr_state_overlaps_summary$tmr_group = factor(tmr_state_overlaps_summary$tmr_group, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Create a plot annotating TMRs and save
tmr_chromatin_state_annotation_plot = ggplot(tmr_state_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black")
tmr_chromatin_state_annotation_plot = customize_ggplot_theme(tmr_chromatin_state_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "chromatin_state", facet_nrow = 2, facet_scales = "free_x", strip_text_size = 17, axis_text_size = 16, axis_title_size = 24, legend_key_size = 1.5, 
  legend_title_size = 24, legend_text_size = 20) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
tmr_chromatin_state_annotation_plot

### Combine The TMR distribution plots and TMR annotation plots 
combined_annotation_plot = tmr_genomic_feature_annotation_plot / tmr_chromatin_state_annotation_plot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
saveRDS(combined_annotation_plot, "combined_annotation_plot.rds")

### Make plots for Roadmap TMRs

# Find overlaps between Roadmap TMRs and regulatory features
roadmap_tmrs = readRDS("../finding_tmrs/tmr_granges/roadmap_tmrs.rds")
roadmap_tmr_annotation_overlaps = data.frame(findOverlaps(roadmap_tmrs, genome_annotation_hg38))
roadmap_tmr_annotation_overlaps$roadmap_tmr_group = "Roadmap"
roadmap_tmr_annotation_overlaps$direction = roadmap_tmrs$direction[roadmap_tmr_annotation_overlaps$queryHits]
roadmap_tmr_annotation_overlaps$annotation = genome_annotation_hg38$region_type[roadmap_tmr_annotation_overlaps$subjectHits]

# Make a table summarizing the overlaps
roadmap_tmr_annotation_overlaps_summary = dplyr::summarise(group_by(roadmap_tmr_annotation_overlaps, roadmap_tmr_group, direction, annotation), count = n())

# Normalize count by MB covered by the regulatory features
roadmap_tmr_annotation_overlaps_summary$normalized_count = roadmap_tmr_annotation_overlaps_summary$count/annotation_widths[roadmap_tmr_annotation_overlaps_summary$annotation]*1e6

# Convert region_type to a factor and give specified order
roadmap_tmr_annotation_overlaps_summary = filter(roadmap_tmr_annotation_overlaps_summary, annotation %in% 
    c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))
roadmap_tmr_annotation_overlaps_summary$annotation = factor(roadmap_tmr_annotation_overlaps_summary$annotation, 
  c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))
roadmap_tmr_annotation_overlaps_summary = data.frame(tidyr::complete(tibble(roadmap_tmr_annotation_overlaps_summary), roadmap_tmr_group, annotation, direction, fill = list(normalized_count = 0, count = 0)))

# Create a plot annotating TMRs and save
roadmap_tmr_genomic_feature_annotation_plot = ggplot(roadmap_tmr_annotation_overlaps_summary, aes(x = roadmap_tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black") 
roadmap_tmr_genomic_feature_annotation_plot = customize_ggplot_theme(roadmap_tmr_genomic_feature_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "annotation", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 14, 
  legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm"), axis.text.x = element_blank()) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
roadmap_tmr_genomic_feature_annotation_plot

# Find overlaps between Roadmap TMRs and chromatin states
roadmap_tmr_state_overlaps = data.frame(findOverlaps(roadmap_tmrs, prostate_18_states_hg38_gr))
roadmap_tmr_state_overlaps$roadmap_tmr_group = "Roadmap"
roadmap_tmr_state_overlaps$direction = roadmap_tmrs$direction[roadmap_tmr_state_overlaps$queryHits]
roadmap_tmr_state_overlaps$chromatin_state = prostate_18_states_hg38_gr$region_type[roadmap_tmr_state_overlaps$subjectHits]

# Make a table summarizing the overlaps
roadmap_tmr_state_overlaps_summary = summarise(group_by(roadmap_tmr_state_overlaps, roadmap_tmr_group, direction, chromatin_state), count = n())

# Normalize count by MB covered by the chromatin state
roadmap_tmr_state_overlaps_summary$normalized_count = roadmap_tmr_state_overlaps_summary$count/chromatin_state_widths[roadmap_tmr_state_overlaps_summary$chromatin_state]*1e6

# Create a plot annotating TMRs and save
roadmap_tmr_chromatin_state_annotation_plot = ggplot(roadmap_tmr_state_overlaps_summary, aes(x = roadmap_tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black")
roadmap_tmr_chromatin_state_annotation_plot = customize_ggplot_theme(roadmap_tmr_chromatin_state_annotation_plot, title = NULL, 
  xlab = NULL, ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "chromatin_state", facet_nrow = 2, facet_scales = "free_x", strip_text_size = 17, axis_text_size = 16, axis_title_size = 24, legend_key_size = 1.5, 
  legend_title_size = 24, legend_text_size = 20) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm"), axis.text.x = element_blank()) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50)))) +
  facet_wrap(~chromatin_state, drop = F, nrow = 2)
roadmap_tmr_chromatin_state_annotation_plot

### Combine The TMR distribution plots and TMR annotation plots 
combined_roadmap_annotation_plot = roadmap_tmr_genomic_feature_annotation_plot / roadmap_tmr_chromatin_state_annotation_plot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = combined_roadmap_annotation_plot, filename = "../figures/supp_figure13.pdf", width = 27, height = 27)
