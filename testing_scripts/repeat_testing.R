# Load required packages
library(GenomicRanges)
library(genomicTools)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/tmr_plot_functions.R")
source("../auxillary_scripts/granges_functions.R")

# Get Granges for transcripts
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr.rds")[readRDS("../auxillary_data/cage_supported_gencode_tss.rds")$ID]
transcripts_gr = methodical::expand_granges(transcripts_gr, 50000, 50000)

# Get repeats 
repeats = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")
names(repeats) = paste0("region_", seq_along(repeats))

# Get the mapability score for repeats
repeat_mappability_scores = genomicTools::bigwig_summarize_over_regions(bw_filepaths = "~/genomes/human/bismap_hg38/k100.bismap.bw",
  gr = repeats, statistic = "mean0", ncores = 1, bed_filepath = NULL, column_names = NULL)

repeats$mappability_score = repeat_mappability_scores$k100.bismap


repeats = subsetByOverlaps(repeats, transcripts_gr)

# Make lists of repeat classes, families and subfamilies
repeat_class_grl = split(repeats, repeats$class)
repeat_families_grl = split(repeats, repeats$family)

repeat_subfamilies_grl = split(repeats[repeats$class != "Simple_repeat"], repeats[repeats$class != "Simple_repeat"]$name)

# Get TMRs
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs.rds")
cpgea_normal_tmrs_repeats = readRDS("tmr_granges/cpgea_normal_tmrs_with_repeats.rds")
cpgea_normal_tmrs_repeats_50kb = readRDS("tmr_granges/cpgea_normal_tmrs_50kb.rds")
cpgea_tumour_tmrs_repeats_50kb  = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb.rds")
names(cpgea_tumour_tmrs_repeats_50kb) = cpgea_tumour_tmrs_repeats_50kb$tmr_name

# Get the maping for TMRs
tmr_mappability_scores = genomicTools::bigwig_summarize_over_regions(bw_filepaths = "~/genomes/human/bismap_hg38/k100.bismap.bw",
  gr = cpgea_tumour_tmrs_repeats_50kb, statistic = "mean0", ncores = 1, bed_filepath = NULL, column_names = NULL)
cpgea_tumour_tmrs_repeats_50kb$mapping = tmr_mappability_scores$k100.bismap

length(subsetByOverlaps(cpgea_normal_tmrs_repeats , repeats, invert = T))

query_gr = cpgea_normal_tmrs_repeats
names(query_gr) = query_gr$tmr_name
subject_gr = repeats

individual_proportion_overlaps = function(query_gr, subject_gr, ignore.strand = T){
  
  # Give an error if any names of query_gr are duplicated
  if(anyDuplicated(names(query_gr))){stop("Names of query_gr cannot be duplicated")}
  
  # Set names for query_gr if they are missing
  if(is.null(names(query_gr))){
    names(query_gr) = paste0("region_", seq_along(query_gr))
  }
  
  # Merge overlapping regions of subject_gr
  subject_gr = reduce(subject_gr, ignore.strand = ignore.strand)
  
  # Find proportion of regions in query_gr involved in each overlap with subject_gr
  overlaps_df = with(data.frame(findOverlaps(query_gr, subject_gr, ignore.strand = ignore.strand)), 
    data.frame(
      name = names(query_gr[queryHits]), 
      overlap = width(pintersect(query_gr[queryHits], subject_gr[subjectHits]))/width(query_gr[queryHits]))
    )
  
  # Get the sum of the proportions for each region from query_gr
  overlaps_df = data.frame(dplyr::summarise(dplyr::group_by(overlaps_df, name), 
    overlap = sum(overlap)))
  
  # Add back in non-overlapping regions
  overlaps_df$name = factor(overlaps_df$name, levels = names(query_gr))
  overlaps_df = tidyr::complete(overlaps_df, name, fill = list(overlap = 0))
  
  # Put in the same order as query_gr and return
  overlaps_df = overlaps_df[match(names(query_gr), overlaps_df$name), ]
  
}

plot_tmrs = function(tmrs){
  
  binned_tmrs =  bin_relative_tmrs(tmrs, 49750)
  
  tmrs_5kb_bins_plot = ggplot(binned_tmrs, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
  tmrs_5kb_bins_plot = customize_ggplot_theme(plot = tmrs_5kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma)) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
  tmrs_5kb_bins_plot
  
}

plot_tmrs(cpgea_tumour_tmrs_repeats_50kb)
plot_tmrs(subsetByOverlaps(cpgea_tumour_tmrs_repeats_50kb, bad, invert = T))
plot_tmrs(subsetByOverlaps(cpgea_tumour_tmrs_repeats_50kb, repeats, invert = T))

plot_tmrs(tmrs = cpgea_tumour_tmrs_repeats_50kb[cpgea_tumour_tmrs_repeats_50kb$mapping < 1])
plot_tmrs(tmrs = cpgea_tumour_tmrs_repeats_50kb[cpgea_tumour_tmrs_repeats_50kb$mapping == 1])

class_overlaps = sapply(repeat_class_grl, function(x) genomicTools::calculate_regions_intersections(gr1 = cpgea_tumour_tmrs_repeats_50kb, x, overlap_measure = "jaccard"))
options(scipen = 30)
sort(class_overlaps)
barplot(sort(class_overlaps))

family_overlaps = sapply(repeat_families_grl, function(x) genomicTools::calculate_regions_intersections(gr1 = cpgea_tumour_tmrs_repeats_50kb, x, overlap_measure = "jaccard"))
sort(family_overlaps)
barplot(sort(family_overlaps))

subfamily_overlaps = sapply(repeat_subfamilies_grl, function(x) genomicTools::calculate_regions_intersections(gr1 = cpgea_tumour_tmrs_repeats_50kb, x, overlap_measure = "jaccard"))
barplot(sort(subfamily_overlaps))


names(cpgea_tumour_tmrs_repeats_50kb) = cpgea_tumour_tmrs_repeats_50kb$tmr_name
cpgea_tumour_tmrs_repeats_50kb_props = individual_proportion_overlaps(cpgea_tumour_tmrs_repeats_50kb, repeats)

x = seq(0, 1, 0.1)
plots = lapply(x, function(y) plot_tmrs(tmrs = cpgea_tumour_tmrs_repeats_50kb[cpgea_tumour_tmrs_repeats_50kb_props$overlap <= y]))
plotR::pdf_save(plots, filename = "~/temp/50kb_plots.pdf")
plot_tmrs(tmrs = cpgea_tumour_tmrs_repeats_50kb[cpgea_tumour_tmrs_repeats_50kb_props$overlap < 0.01])

# Read in bsimap bedGraph
bismap_bed = data.table::fread("k100.bismap.bedGraph")
names(bismap_bed) = c("seqnames", "start", "end", "score")
bismap_gr = makeGRangesFromDataFrame(bismap_bed, starts.in.df.are.0based = T, keep.extra.columns = T)
saveRDS(bismap_gr, "~/genomes/human/bismap_hg38/bismap_gr.rds")
bismap_gr = readRDS("~/genomes/human/bismap_hg38/bismap_gr.rds")
bad = bismap_gr[bismap_gr$score < 1]

# Expand non_high_mappability_regions 
unmappable_regions = methodical::expand_granges(unmappable_regions, 150, 150)

cpgea_normal_tmrs_50kb = readRDS("tmr_granges/cpgea_normal_tmrs_50kb_unfiltered.rds")
cpgea_tumour_tmrs_50kb = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb_unfiltered.rds")
mcrpc_tmrs_50kb = readRDS("tmr_granges/mcrpc_tmrs_50kb_unfiltered.rds")

cpgea_noarmal_tmrs_50kb = subsetByOverlaps(cpgea_normal_tmrs_50kb, unmappable_regions, invert = T)
cpgea_tumour_tmrs_50kb = subsetByOverlaps(cpgea_tumour_tmrs_50kb, unmappable_regions, invert = T)
mcrpc_tmrs_50kb = subsetByOverlaps(mcrpc_tmrs_50kb, unmappable_regions, invert = T)
