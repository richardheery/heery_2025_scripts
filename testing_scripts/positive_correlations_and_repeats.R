# Load TMRs after filtering for coverage
cpgea_normal_tmrs = readRDS("new_tmr_granges/cpgea_normal_tmrs_50kb.rds")
cpgea_tumour_tmrs = readRDS("new_tmr_granges/cpgea_tumour_tmrs_50kb.rds")
mcrpc_tmrs = readRDS("new_tmr_granges/mcrpc_tmrs_50kb.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_coverage_filtered = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_coverage_filtered = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_coverage_filtered, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.50, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_50kb_bins_plot_coverage_filtered

### Make plot without repeats

repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

cpgea_normal_tmrs = subsetByOverlaps(cpgea_normal_tmrs, repeat_ranges, invert = T)
cpgea_tumour_tmrs = subsetByOverlaps(cpgea_tumour_tmrs, repeat_ranges, invert = T)
mcrpc_tmrs = subsetByOverlaps(mcrpc_tmrs, repeat_ranges, invert = T) 

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_coverage_filtered_no_lines = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_coverage_filtered_no_lines = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_coverage_filtered_no_lines, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.50, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
x11(); tmrs_50kb_bins_plot_coverage_filtered_no_lines

### Filter just for LINEs

line_ranges = repeat_ranges[repeat_ranges$class == "SINE"]

cpgea_normal_tmrs = subsetByOverlaps(cpgea_normal_tmrs, line_ranges, invert = T)
cpgea_tumour_tmrs = subsetByOverlaps(cpgea_tumour_tmrs, line_ranges, invert = T)
mcrpc_tmrs = subsetByOverlaps(mcrpc_tmrs, line_ranges, invert = T) 

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_coverage_filtered_no_lines = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_coverage_filtered_no_lines = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_coverage_filtered_no_lines, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.50, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_50kb_bins_plot_coverage_filtered_no_lines

##### New TMRs

# Load TMRs after filtering for coverage
cpgea_normal_tmrs = readRDS("tmr_granges/cpgea_normal_tmrs_50kb.rds")
cpgea_tumour_tmrs = readRDS("tmr_granges/cpgea_tumour_tmrs_50kb.rds")
mcrpc_tmrs = readRDS("tmr_granges/mcrpc_tmrs_50kb.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_standard = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_standard = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_standard, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.50, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_50kb_bins_plot_standard

### Make plot without repeats

repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

cpgea_normal_tmrs = subsetByOverlaps(cpgea_normal_tmrs, repeat_ranges, invert = T)
cpgea_tumour_tmrs = subsetByOverlaps(cpgea_tumour_tmrs, repeat_ranges, invert = T)
mcrpc_tmrs = subsetByOverlaps(mcrpc_tmrs, repeat_ranges, invert = T) 

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 49750)
cpgea_tumour_50kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 49750)
mcrpc_tmr_50kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 49750)

# Combine TMR distributions
combined_50kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions,
  mcrpc = mcrpc_tmr_50kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_standard_no_repeats = ggplot(combined_50kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_50kb_bins_plot_standard_no_repeats = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_standard_no_repeats, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.50, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
x11(); tmrs_50kb_bins_plot_standard_no_repeats

map_track = data.table::fread("~/genomes/human/bismap_hg38/k100.bismap.bed.gz", sep = "\t", skip = 1)
names(map_track) = c("seqnames", "start", "end", "name", "score", "strand")
map_track_gr = makeGRangesFromDataFrame(map_track, starts.in.df.are.0based = T, keep.extra.columns = T)

hg38 = as(seqinfo(BSgenome.Hsapiens.UCSC.hg38), "GRanges")[1:24]
non_covered = setdiff(hg38, map_track_gr, ignore.strand = T)

### Repeat and non-covered overlaps

cpgea_tumour_tmrs = readRDS("new_tmr_granges/cpgea_tumour_tmrs_50kb.rds")
cpgea_tumour_tmrs_positive = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Positive"]
cpgea_tumour_tmrs_negative = cpgea_tumour_tmrs[cpgea_tumour_tmrs$direction == "Negative"]

length(subsetByOverlaps(cpgea_tumour_tmrs_negative, repeat_ranges))/length(cpgea_tumour_tmrs_negative)
#
length(subsetByOverlaps(cpgea_tumour_tmrs_positive, repeat_ranges))/length(cpgea_tumour_tmrs_positive)

repeat_ranges = readRDS("../auxillary_data/repeatmasker_granges_ucsc.rds")

# Filter repeats for thosw within 50KB of TSS
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_tss_gr.rds")
repeats_near_transcripts = subsetByOverlaps(repeat_ranges, expand_granges(transcripts_gr, 50000, 50000))
subfamilies_1000 = names(which(table(repeats_near_transcripts$name) > 1000))

repeats_near_transcripts = repeats_near_transcripts[repeats_near_transcripts$name %in% subfamilies_1000]
repeats_near_transcripts_list = split(repeats_near_transcripts, repeats_near_transcripts$name)

system.time({positive_repeat_jaccards = sapply(repeats_near_transcripts_list, function(x)
  genomicTools::calculate_regions_intersections(x, cpgea_tumour_tmrs_positive, overlap_measure = "jaccard"))})

system.time({positive_repeat_proportions = sapply(repeats_near_transcripts_list, function(x)
  genomicTools::calculate_regions_intersections(x, cpgea_tumour_tmrs_positive, overlap_measure = "proportion"))})

system.time({negative_repeat_proportions = sapply(repeats_near_transcripts_list, function(x)
  genomicTools::calculate_regions_intersections(x, cpgea_tumour_tmrs_negative, overlap_measure = "proportion"))})

positive_tmr_repeat_overlaps = sapply(repeats_near_transcripts_list, function(x) 
  length(subsetByOverlaps(cpgea_tumour_tmrs_positive, x, ignore.strand = T))/length(cpgea_tumour_tmrs_positive))
negative_tmr_repeat_overlaps = sapply(repeats_near_transcripts_list, function(x) 
  length(subsetByOverlaps(cpgea_tumour_tmrs_negative, x, ignore.strand = T))/length(cpgea_tumour_tmrs_negative))

