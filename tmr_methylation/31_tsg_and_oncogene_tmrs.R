# Find TMRs overlapping TSGs and oncogenes

# Load required packages
library(GenomicRanges)
library(methodical)
library(dplyr)
library(ggpubr)

# Get GRanges with TSS
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")

# Get Illumina probes ranges for hg38
illumina_probes_hg38 = readRDS("../auxillary_data/infinium_450k_probe_granges_hg38.rds")
illumina_probes_hg38$region_type = "probe"
illumina_probes_hg38 = promoters(illumina_probes_hg38, 5, 6)

# Load meth RSEs for CPGEA and MCRPC
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse/")
mcrpc_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse/")

# Combine into a single RSE
colData(cpgea_meth_rse) = NULL; colData(mcrpc_meth_rse) = NULL; metadata(mcrpc_meth_rse) = list()
mcrpc_meth_rse = sort(mcrpc_meth_rse)
combined_meth_rse = cbind(cpgea_meth_rse, mcrpc_meth_rse)

# Get TSGs and oncogenes from COSMIC
tsgs = readRDS("../auxillary_data/cosmic_all_tsgs.rds")
oncogenes = readRDS("..//auxillary_data/cosmic_all_oncogenes.rds")

# Get TMRs for CPGEA and MCRPC
cpgea_tumour_tmrs = readRDS("../finding_tmrs/tmr_granges/cpgea_tumour_tmrs.rds")
cpgea_tumour_tmrs = cpgea_tumour_tmrs[abs(cpgea_tumour_tmrs$distance_to_tss) < 5000]
mcrpc_tmrs = readRDS("../finding_tmrs/tmr_granges/mcrpc_tmrs.rds")
mcrpc_tmrs = mcrpc_tmrs[abs(mcrpc_tmrs$distance_to_tss) < 5000]

# Find TMRs overlapping TSGs and oncogenes
tsg_cpgea_tumour_tmrs = cpgea_tumour_tmrs[cpgea_tumour_tmrs$gene_name %in% tsgs]
oncogene_cpgea_tumour_tmrs = cpgea_tumour_tmrs[cpgea_tumour_tmrs$gene_name %in% oncogenes]
tsg_mcrpc_tmrs = mcrpc_tmrs[mcrpc_tmrs$gene_name %in% tsgs]
oncogene_mcrpc_tmrs = mcrpc_tmrs[mcrpc_tmrs$gene_name %in% oncogenes]

# Create a function which shows methylation change at CpG around TMRs associated with TSGs and oncogenes
plot_tmr_methylation = function(transcript, samples){
  
  if(samples == "cpgea"){
    tumour_pattern = "T"
    tmrs_gr = cpgea_tumour_tmrs
    ylabel = "Prostate Tumour\nMethylation Change"
  } else if(samples == "mcrpc"){
    tumour_pattern = "DTB"
    tmrs_gr = mcrpc_tmrs
    ylabel = "Prostate Metastases\nMethylation Change"
  }
  
  # 
  tss = tss_gr[tss_gr$ID == transcript]
  tss_expanded = promoters(tss, upstream = 5000, downstream = 5001)
  
  # Get methylation values around TSS
  tss_methylation_values = extractGRangesMethSiteValues(meth_rse = combined_meth_rse, genomic_regions = tss_expanded)
  tss_methylation_values = rowMeans(select(tss_methylation_values, starts_with(tumour_pattern)), na.rm = T) - 
    rowMeans(select(tss_methylation_values, starts_with("N")), na.rm = T)
  tss_methylation_values = setNames(data.frame(tss_methylation_values), "meth_change")
  
  meth_change_plot = plotMethylationValues(meth_site_values = tss_methylation_values, sample_name = "meth_change", 
    reference_tss = tss, title = NULL, ylabel = ylabel, 
    xlabel = sprintf("Distance to *%s* TSS (bp)", tss$gene_name)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(hjust = 0.5, size = 20), )
  plotTMRs(meth_site_plot = meth_change_plot, tmrs_gr = tmrs_gr[tmrs_gr$ID == transcript], reference_tss = tss) +
  theme(legend.position = "right") + guides(fill = "none")
  
}
c("ENST00000398606", "ENST00000359013", "ENST00000327367", "ENST00000411767")
# Create plots of methylation change for TSGs: GSTP1 and TGFBR2 in CPGEA samples and PTCH1 and PTPN13 in MCRPC samples
system.time({tsg_tmr_plots_cpgea = lapply(c("ENST00000398606", "ENST00000359013", "ENST00000327367", "ENST00000411767"), function(x) 
  {print(x); plot_tmr_methylation(x, samples = "cpgea")})})
plotR::pdf_save(tsg_tmr_plots_cpgea, filename = "tsg_tmr_plots_cpgea.pdf")

sort(unique(tsg_cpgea_tumour_tmrs$ID))
system.time({tsg_tmr_plots_mcrpc = lapply(sort(unique(tsg_mcrpc_tmrs$ID)), function(x) 
  {print(x); plot_tmr_methylation(x, samples = "mcrpc")})})
plotR::pdf_save(tsg_tmr_plots_mcrpc, filename = "tsg_tmr_plots_mcrpc.pdf")


# Combine CPGEA and MCRPC plots
combined_tsg_plot = ggarrange(plotlist = tsg_tmr_plots_cpgea, nrow = 2, ncol = 2)
ggsave(plot = combined_tsg_plot, "combined_tsg_plot.pdf",device = cairo_pdf, width = 24, height = 13.5)

ENST00000267101, ENST00000354725
# Create plots of methylation change for TSGs: FOXA1 and MYC in CPGEA samples and PTCH1 and PTPN13 in MCRPC samples
system.time({oncogene_tmr_plots_cpgea = lapply(c("ENST00000267101", "ENST00000354725"), function(x) 
  {print(x); plot_tmr_methylation(x, samples = "cpgea")})})
plotR::pdf_save(oncogene_tmr_plots_cpgea, filename = "oncogene_tmr_plots_cpgea.pdf")


system.time({oncogene_tmr_plots_mcrpc = lapply(c("ENST00000536559", "ENST00000683015"), function(x) 
  {print(x); plot_tmr_methylation(x, samples = "mcrpc")})})

# Combine CPGEA and MCRPC plots
combined_oncogene_plot_list = c(oncogene_tmr_plots_cpgea, oncogene_tmr_plots_mcrpc)
combined_oncogene_plot = ggarrange(plotlist = oncogene_tmr_plots_cpgea, nrow = 2, ncol = 1)
ggsave(plot = combined_oncogene_plot, "combined_oncogene_plot.pdf",device = cairo_pdf, width = 24, height = 13.5)

# 
complete_plot_list = c(tsg_tmr_plots_cpgea, oncogene_tmr_plots_cpgea)
complete_plot_list = lapply(complete_plot_list, function(x) x + 
    theme(legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.key.size = unit(2, "cm")))
complete_plot = ggarrange(plotlist = complete_plot_list, nrow = 3, ncol = 2, common.legend = T, legend = "top")
ggsave(plot = complete_plot, "../figures/supp_figure16.pdf",device = cairo_pdf, width = 24, height = 20.25)

presentation_plots = ggarrange(plotlist = complete_plot_list[c(1, 3, 5, 7)], nrow = 2, ncol = 2, common.legend = T, legend = "top")
ggsave(plot = presentation_plots, "~/promoter_project/presentation_figures/tsg_and_oncogene_plots.pdf", device = cairo_pdf, width = 24, height = 13.5)

####
complete_plot_list = lapply(complete_plot_list, function(x) 
  annotateMethSitePlot(meth_site_plot = x, annotation_gr = illumina_probes_hg38, reference_tss = T))


# BCOR 113, PATZ1 (102), CEBPA (94), AXIN2 (88), NF1 (87), CDX2 (76), SH2B3 (74), ETV6 (65), FAS (56), # WNK2 (45), PPARG (9), 
# BCOR, CEBPA, NF1, CDX2, FAS and GSTP1, PPARG
# FAS, NF1 and CEBPA

# MCRPC: PPARG (8), FBLN2 (9), IGFBP2 (10), SLC34A2 (11), PTPN13 (12), IKZF1 (19), PCTH1 (31), 
# PTCH1, PTPN13, 
system.time({tsg_tmr_plots = lapply(seq_along(tsg_cpgea_tumour_tmrs), function(x) 
  {print(x); plot_tmr_methylation(tsg_cpgea_tumour_tmrs$transcript_id[x], samples = "cpgea")})})
plotR::pdf_save(tsg_tmr_plots, filename = "tsg_and_oncogene_plots/tsg_tmr_meth_change_plots.pdf")

# AFF3 (8), FOXP1 (23), CTNND2 (37), AFDN (47), ERBB3 (87), SIX1 (96)
# CTNND2
# system.time({oncogene_tmr_plots = lapply(seq_along(oncogene_cpgea_tumour_tmrs), function(x) 
#   {print(x); plot_tmr_methylation(oncogene_cpgea_tumour_tmrs$transcript_id[x])})})
# plotR::pdf_save(oncogene_tmr_plots, filename = "tsg_and_oncogene_plots/oncogene_tmr_meth_change_plots.pdf")

# Check MCRPC samples
system.time({tsg_tmr_plots_mcrpc = lapply(seq_along(tsg_mcrpc_tmrs), function(x) 
  {print(x); plot_tmr_methylation(tsg_mcrpc_tmrs$transcript_id[x], samples = "mcrpc")})})
plotR::pdf_save(tsg_tmr_plots_mcrpc, filename = "tsg_and_oncogene_plots/tsg_tmr_meth_change_plots_mcrpc.pdf")

system.time({oncogene_tmr_plots_mcrpc = lapply(seq_along(oncogene_mcrpc_tmrs), function(x) 
  {print(x); plot_tmr_methylation(oncogene_mcrpc_tmrs$transcript_id[x], samples = "mcrpc")})})
plotR::pdf_save(oncogene_tmr_plots_mcrpc, filename = "tsg_and_oncogene_plots/oncogene_tmr_meth_change_plots_mcrpc.pdf")
"finsihed"