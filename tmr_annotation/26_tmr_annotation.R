
# Load required packages
library(dplyr)
library(GenomicRanges)
library(regioneR)
library(patchwork)
source("../auxillary_scripts/plotting_functions.R")

# Get a GRanges for transcripts and expand
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")
transcripts_gr = methodical::expand_granges(transcripts_gr, 5000, 5000)

tmr_list = readRDS("../finding_tmrs/tmr_granges/tmr_list.rds")

genome = GRanges(seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38))[1:24]
non_transcript_regions = GenomicRanges::setdiff(genome, transcripts_gr, ignore.strand = T)

# Load hg38 genome annotation
genome_annotation_hg38 = readRDS("../auxillary_data/complete_regulatory_annotation.rds")
genome_annotation_hg38$region_type = gsub(" ", "_", genome_annotation_hg38$region_type)

# Convert genome_annotation_hg38 to a list
genome_annotation_hg38_list = GRangesList(split(genome_annotation_hg38, genome_annotation_hg38$region_type))

tmr_overlap_permutation_test = function(tmrs, regions_grl, n){
  
  # Calculate size of the overlaps between TMRs and annotation regions
  tmr_overlaps = sapply(regions_grl, function(x) sum(width(GenomicRanges::intersect(x, tmrs, ignore.strand = T))))
  
  # Flatten regions_grl
  regions_gr = unlist(GRangesList(regions_grl))
  
  # Create a GRangeLsit with randomized TMRs
  system.time({random_tmrs = GRangesList(parallel::mclapply(seq.int(n), function(x)
    randomizeRegions(A = tmrs, genome = genome, mask = non_transcript_regions, allow.overlaps = T, per.chromosome = F), mc.cores = 18))})
  
  # Name GRanges with permuation number and convert to a flat GRanges
  names(random_tmrs) = paste0("permutation_", seq_along(random_tmrs))
  random_tmrs = unlist(random_tmrs)
  names(random_tmrs) = gsub("\\..*", "", names(random_tmrs))
  
  # Find overlaps between random_tmrs and regions_gr
  random_overlaps = data.frame(findOverlaps(random_tmrs, regions_gr, ignore.strand = T))
  
  # Add permutation number and region type
  random_overlaps$permutation = names(random_tmrs)[random_overlaps$queryHits]
  random_overlaps$region_type = names(regions_gr)[random_overlaps$subjectHits]
  
  # Find size of the intersections
  random_overlaps$intersection = width(pintersect(random_tmrs[random_overlaps$queryHits], regions_gr[random_overlaps$subjectHits]))
  
  # Set region_type as a factor
  random_overlaps$region_type = factor(random_overlaps$region_type, levels = names(regions_grl))
  
  # Summarize the intersections
  random_overlaps_summary = data.frame(dplyr::summarise(
    dplyr::group_by(random_overlaps, permutation, region_type, .drop = FALSE), 
      intersection = sum(intersection)))
  
  # Convert into a wide data.frame and ensure columns are in the same order as regions_grl
  random_overlaps_summary = tidyr::pivot_wider(random_overlaps_summary, names_from = "region_type", values_from = "intersection")
  random_overlaps_summary = data.frame(tibble::column_to_rownames(random_overlaps_summary, "permutation"))
  random_overlaps_summary = random_overlaps_summary[, names(regions_grl)]
  
  # Calculate p-values
  p_value = (rowSums(apply(random_overlaps_summary, 1, function(x) x >= tmr_overlaps)) + 1)/(n+1)
  
  # Create a data.frame with the results and return
  results_df = data.frame(
    region_type = names(tmr_overlaps),
    tmr_overlaps = tmr_overlaps,
    random_overlaps = colMeans(random_overlaps_summary, na.rm =T),
    enrichment = tmr_overlaps/colMeans(random_overlaps_summary, na.rm = T),
    p_value = p_value, row.names = NULL
  )
  return(results_df)
  
}

# 
set.seed(123)
system.time({regulatory_regions_randomization_results = lapply(tmr_list, function(x) tmr_overlap_permutation_test(x, regions_grl = genome_annotation_hg38_list, n = 1000))})
saveRDS(regulatory_regions_randomization_results, "regulatory_regions_randomization_results.rds")
regulatory_regions_randomization_results = readRDS("regulatory_regions_randomization_results.rds")

# Combine results into a single table
regulatory_regions_randomization_results = dplyr::bind_rows(regulatory_regions_randomization_results, .id = "tmr_group")

# Remove introns and exons
regulatory_regions_randomization_results = dplyr::filter(regulatory_regions_randomization_results, !region_type %in% c("Intron", "Exon", "Open_Sea"))

# Correct p-values and add significance symbols
regulatory_regions_randomization_results$q_value = p.adjust(regulatory_regions_randomization_results$p_value, method = "fdr")
regulatory_regions_randomization_results$significance = sig_sym(regulatory_regions_randomization_results$q_value)

# Add name of dataset and TMR direction
regulatory_regions_randomization_results$dataset = gsub("_tmrs_.*", "", regulatory_regions_randomization_results$tmr_group)
regulatory_regions_randomization_results$direction = gsub(".*_tmrs_", "", regulatory_regions_randomization_results$tmr_group)

# Set levels of region_type
regulatory_regions_randomization_results$region_type = gsub("_", " ", regulatory_regions_randomization_results$region_type)
regulatory_regions_randomization_results$region_type = factor(regulatory_regions_randomization_results$region_type, levels = sort(unique(regulatory_regions_randomization_results$region_type)))

# Create a barplot showing the enrichment of region classes among TMRs
regulatory_regions_enrichment_barplot = ggplot(regulatory_regions_randomization_results, aes(x = dataset, y = enrichment, fill = direction)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(mapping = aes(label = significance, y = enrichment + 0.05, x = dataset, group = direction), position = position_dodge(width = 0.9), size = 10)
regulatory_regions_enrichment_barplot = customize_ggplot_theme(regulatory_regions_enrichment_barplot, xlab = "Dataset", ylab = "Enrichment",
  x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), x_labels_angle = 55, fill_labels = c("Negative", "Positive"), fill_colors = c("#A28CB1", "#D2C465"),
  facet = "region_type", facet_nrow = 1, facet_scales = "fixed") + 
  theme(strip.background = element_blank())
regulatory_regions_enrichment_barplot

###

# Get a GRanges with chromatin states for prostate 
prostate_18_states_hg38_gr = readRDS("../auxillary_data/prostate_18_states_hg38_gr.rds")

# Split prostate_18_states_hg38_gr into a list
prostate_18_states_hg38_gr$region_type = prostate_18_states_hg38_gr$description
prostate_18_states_hg38_grl= split(prostate_18_states_hg38_gr, prostate_18_states_hg38_gr$region_type)
names(prostate_18_states_hg38_grl) = gsub(" ", "_", names(prostate_18_states_hg38_grl))
names(prostate_18_states_hg38_grl) = gsub("/", "_", names(prostate_18_states_hg38_grl))
names(prostate_18_states_hg38_grl) = gsub("'", "", names(prostate_18_states_hg38_grl))
names(prostate_18_states_hg38_grl) = gsub("&", "_", names(prostate_18_states_hg38_grl))

set.seed(123)
system.time({chromatin_states_randomization_results = lapply(tmr_list, function(x) tmr_overlap_permutation_test(tmrs = x, regions_grl = prostate_18_states_hg38_grl, n = 1000))})
saveRDS(chromatin_states_randomization_results , "chromatin_states_randomization_results.rds")

chromatin_states_randomization_results = readRDS("chromatin_states_randomization_results.rds")

# Combine results into a single table
chromatin_states_randomization_results = dplyr::bind_rows(chromatin_states_randomization_results, .id = "tmr_group")

# Remove introns and exons
chromatin_states_randomization_results = dplyr::filter(chromatin_states_randomization_results, !region_type %in% c("Intron", "Exon"))

# Correct p-values and add significance symbols
chromatin_states_randomization_results$q_value = p.adjust(chromatin_states_randomization_results$p_value, method = "fdr")
chromatin_states_randomization_results$significance = sig_sym(chromatin_states_randomization_results$q_value)

# Add name of dataset and TMR direction
chromatin_states_randomization_results$dataset = gsub("_tmrs_.*", "", chromatin_states_randomization_results$tmr_group)
chromatin_states_randomization_results$direction = gsub(".*_tmrs_", "", chromatin_states_randomization_results$tmr_group)

# Set levels of region_type
chromatin_states_randomization_results$region_type = gsub("_", " ", chromatin_states_randomization_results$region_type)
chromatin_states_randomization_results$region_type = factor(chromatin_states_randomization_results$region_type, levels = chromatin_states_randomization_results$region_type[1:18])

# Shorten "PolyComb" in description to "PC" 
levels(chromatin_states_randomization_results$region_type)[3] = "Flanking TSS 5'"
levels(chromatin_states_randomization_results$region_type)[4] = "Flanking TSS 3'"
levels(chromatin_states_randomization_results$region_type)[5] = "Strong Transcription"
levels(chromatin_states_randomization_results$region_type)[6] = "Weak Transcription"
levels(chromatin_states_randomization_results$region_type)[7] = "Genic Enhancer 1"
levels(chromatin_states_randomization_results$region_type)[8] = "Genic Enhancer 2"
levels(chromatin_states_randomization_results$region_type)[12] = "ZNF Genes/Repeats"
levels(chromatin_states_randomization_results$region_type)[14] = "Bivalent/Poised TSS"
levels(chromatin_states_randomization_results$region_type)[16] = "Repressed PC"
levels(chromatin_states_randomization_results$region_type)[17] = "Weak Repressed PC"

# Create a barplot showing the enrichment of region classes among TMRs
chromatin_states_enrichment_barplot = ggplot(chromatin_states_randomization_results, aes(x = dataset, y = enrichment, fill = direction)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(mapping = aes(label = significance, y = enrichment + 0.05, x = dataset, group = direction), position = position_dodge(width = 0.9), size = 10)
chromatin_states_enrichment_barplot = customize_ggplot_theme(chromatin_states_enrichment_barplot, xlab = "Dataset", ylab = "Enrichment",
  x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), x_labels_angle = 55, fill_labels = c("Negative", "Positive"), fill_colors = c("#A28CB1", "#D2C465"),
  facet = "region_type", facet_nrow = 2, facet_scales = "fixed") + 
  theme(strip.background = element_blank())
chromatin_states_enrichment_barplot

### Combine The TMR distribution plots and TMR annotation plots 
combined_annotation_plot = regulatory_regions_enrichment_barplot / chromatin_states_enrichment_barplot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
saveRDS(combined_annotation_plot, "combined_annotation_plot.rds")

### Make plots for Roadmap TMRs

# Find overlaps between Roadmap TMRs and roadmap_regulatory features
roadmap_tmrs = readRDS("../finding_tmrs/tmr_granges/roadmap_tmrs.rds")
roadmap_tmrs_list = split(roadmap_tmrs, roadmap_tmrs$direction)

set.seed(123)
system.time({roadmap_tmrs_regulatory_regions_randomization_results = lapply(roadmap_tmrs_list, function(x)
  tmr_overlap_permutation_test(x, regions_grl = genome_annotation_hg38_list, n = 1000))})
saveRDS(roadmap_tmrs_regulatory_regions_randomization_results, "roadmap_tmrs_regulatory_regions_randomization_results.rds")

# Load regulatory permutation test results
roadmap_regulatory_regions_randomization_results = readRDS("roadmap_tmrs_regulatory_regions_randomization_results.rds")

# Combine results into a single table
roadmap_regulatory_regions_randomization_results = dplyr::bind_rows(roadmap_regulatory_regions_randomization_results, .id = "tmr_group")

# Remove introns and exons
roadmap_regulatory_regions_randomization_results = dplyr::filter(roadmap_regulatory_regions_randomization_results, !region_type %in% c("Intron", "Exon", "Open_Sea"))

# Correct p-values and add significance symbols
roadmap_regulatory_regions_randomization_results$q_value = p.adjust(roadmap_regulatory_regions_randomization_results$p_value, method = "fdr")
roadmap_regulatory_regions_randomization_results$significance = sig_sym(roadmap_regulatory_regions_randomization_results$q_value)

# Add name of dataset and TMR direction
roadmap_regulatory_regions_randomization_results$dataset = gsub("_tmrs_.*", "", roadmap_regulatory_regions_randomization_results$tmr_group)
roadmap_regulatory_regions_randomization_results$direction = gsub(".*_tmrs_", "", roadmap_regulatory_regions_randomization_results$tmr_group)

# Set levels of region_type
roadmap_regulatory_regions_randomization_results$region_type = gsub("_", " ", roadmap_regulatory_regions_randomization_results$region_type)
roadmap_regulatory_regions_randomization_results$region_type = factor(roadmap_regulatory_regions_randomization_results$region_type, levels = sort(unique(roadmap_regulatory_regions_randomization_results$region_type)))

# Create a barplot showing the enrichment of region classes among TMRs
roadmap_regulatory_regions_enrichment_barplot = ggplot(roadmap_regulatory_regions_randomization_results, aes(x = dataset, y = enrichment, fill = direction)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(mapping = aes(label = significance, y = enrichment + 0.05, x = dataset, group = direction), position = position_dodge(width = 0.9), size = 10)
roadmap_regulatory_regions_enrichment_barplot = customize_ggplot_theme(roadmap_regulatory_regions_enrichment_barplot, xlab = "Dataset", ylab = "Enrichment",
  x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), x_labels_angle = 55, fill_labels = c("Negative", "Positive"), fill_colors = c("#A28CB1", "#D2C465"),
  facet = "region_type", facet_nrow = 1, facet_scales = "fixed") + 
  theme(strip.background = element_blank())
roadmap_regulatory_regions_enrichment_barplot

set.seed(123)
system.time({roadmap_tmrs_chromatin_states_randomization_results = lapply(roadmap_tmrs_list, function(x)
  tmr_overlap_permutation_test(x, regions_grl = prostate_18_states_hg38_grl, n = 1000))})
saveRDS(roadmap_tmrs_chromatin_states_randomization_results, "roadmap_tmrs_chromatin_states_randomization_results.rds")

#
roadmap_chromatin_states_randomization_results = readRDS("roadmap_tmrs_chromatin_states_randomization_results.rds")

# Combine results into a single table
roadmap_chromatin_states_randomization_results = dplyr::bind_rows(roadmap_chromatin_states_randomization_results, .id = "tmr_group")

# Correct p-values and add significance symbols
roadmap_chromatin_states_randomization_results$q_value = p.adjust(roadmap_chromatin_states_randomization_results$p_value, method = "fdr")
roadmap_chromatin_states_randomization_results$significance = sig_sym(roadmap_chromatin_states_randomization_results$q_value)

# Add name of dataset and TMR direction
roadmap_chromatin_states_randomization_results$direction = roadmap_chromatin_states_randomization_results$tmr_group

# Set levels of region_type
roadmap_chromatin_states_randomization_results$region_type = gsub("_", " ", roadmap_chromatin_states_randomization_results$region_type)
roadmap_chromatin_states_randomization_results$region_type = factor(roadmap_chromatin_states_randomization_results$region_type, 
  levels = roadmap_chromatin_states_randomization_results$region_type[1:18])

# Shorten "PolyComb" in description to "PC" 
levels(roadmap_chromatin_states_randomization_results$region_type)[3] = "Flanking TSS 5'"
levels(roadmap_chromatin_states_randomization_results$region_type)[4] = "Flanking TSS 3'"
levels(roadmap_chromatin_states_randomization_results$region_type)[5] = "Strong Transcription"
levels(roadmap_chromatin_states_randomization_results$region_type)[6] = "Weak Transcription"
levels(roadmap_chromatin_states_randomization_results$region_type)[7] = "Genic Enhancer 1"
levels(roadmap_chromatin_states_randomization_results$region_type)[8] = "Genic Enhancer 2"
levels(roadmap_chromatin_states_randomization_results$region_type)[12] = "ZNF Genes/Repeats"
levels(roadmap_chromatin_states_randomization_results$region_type)[14] = "Bivalent/Poised TSS"
levels(roadmap_chromatin_states_randomization_results$region_type)[16] = "Repressed PC"
levels(roadmap_chromatin_states_randomization_results$region_type)[17] = "Weak Repressed PC"

# Create a barplot showing the enrichment of region classes among TMRs
roadmap_chromatin_states_enrichment_barplot = ggplot(roadmap_chromatin_states_randomization_results, aes(x = direction, y = enrichment, fill = direction)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(mapping = aes(label = significance, y = enrichment + 0.05, x = direction, group = direction), position = position_dodge(width = 0.9), size = 10)
roadmap_chromatin_states_enrichment_barplot = customize_ggplot_theme(roadmap_chromatin_states_enrichment_barplot, xlab = "Dataset", ylab = "Enrichment",
  x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), x_labels_angle = 55, fill_labels = c("Negative", "Positive"), fill_colors = c("#A28CB1", "#D2C465"),
  facet = "region_type", facet_nrow = 2, facet_scales = "fixed") + 
  theme(strip.background = element_blank())
roadmap_chromatin_states_enrichment_barplot

### Combine The TMR distribution plots and TMR annotation plots 
combined_roadmap_annotation_plot = roadmap_regulatory_regions_enrichment_barplot / roadmap_chromatin_states_enrichment_barplot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = combined_roadmap_annotation_plot, filename = "../figures/supp_figure13.pdf", width = 27, height = 27)

