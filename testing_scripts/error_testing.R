meth_rse = cpgea_meth_rse; 
  transcript_expression_table = cpgea_kallisto_deseq2_counts; samples_subset = normal_samples; tss_gr = tss_gr; tss_associated_gr = transcripts_gr; 
  cor_method = "spearman"; min_number_complete_pairs = 30; BPPARAM = bpparam; add_distance_to_region = T

meth_rse; assay_number = 1; transcript_expression_table; 
  samples_subset = NULL; tss_gr; tss_associated_gr; cor_method = "pearson"; 
  add_distance_to_region = TRUE; max_sites_per_chunk = NULL; BPPARAM = BiocParallel::bpparam(); results_dir = NULL

meth_rse = cpgea_meth_rse; 
  transcript_expression_table = cpgea_kallisto_deseq2_counts; samples_subset = normal_samples; tss_gr = tss_gr; tss_associated_gr = transcripts_gr; 
  cor_method = "spearman"; min_number_complete_pairs = 30; BPPARAM = bpparam; add_distance_to_region = T

    # Create an iterator function for TSS sites
    tss_iter <- methodical:::.tss_iterator(meth_values_chunk; tss_region_indices_list; transcript_values_list; 
      tss_gr_chunk_list; cor_method; add_distance_to_region; results_dir)
    
table1 = meth_table; table2 = transcript_table; 
        table1_name = "meth_site"; table2_name = "transcript_name"; 
        cor_method = cor_method; p_adjust_method = "none"

cor_method = "pearson"; table1_name = "table1"; table2_name = "table2"; p_adjust_method = "BH"; n_covariates = 0

table1 = meth_table; table2 = transcript_table; 
        table1_name = "meth_site"; table2_name = "transcript_name"; 
        cor_method = cor_method; p_adjust_method = "none"

correlation_objects = readRDS("correlation_objects_test.rds")
attach(correlation_objects)
meth_table <- t(meth_table)
transcript_table <- setNames(data.frame(t(transcript_table)), transcript_name)
table1 = meth_table; table2 = transcript_table; cor_method = "pearson"; table1_name = "table1"; table2_name = "table2"; p_adjust_method = "BH"; n_covariates = 0

###
library(EPIC)
library(methodical)
library

# Get transcript counts for CPGEA samples
cpgea_transcript_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_transcript_counts.tsv.gz", sep = "\t"), row.names = 1)
mcrpc_transcript_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_transcript_counts.tsv.gz", sep = "\t"), row.names = 1)
normal = grep("N", names(cpgea_transcript_counts))
tumour = grep("T", names(cpgea_transcript_counts))

# Creates genes to transcript list
tss_gr = readRDS("../auxillary_data/cage_supported_gencode_tss.rds")
gene_to_transcript_list = split(tss_gr$ID, tss_gr$gene_name)

gene_expression = methodical::sumTranscriptValuesForGenes(cpgea_transcript_counts, gene_to_transcript_list)
system.time({test = EPIC(gene_expression)})
fivenum(test$cellFractions[normal, "otherCells"])
fivenum(test$cellFractions[tumour, "otherCells"])

gene_expression_mcrpc = methodical::sumTranscriptValuesForGenes(mcrpc_transcript_counts, gene_to_transcript_list)
system.time({test2 = EPIC(gene_expression_mcrpc)})
