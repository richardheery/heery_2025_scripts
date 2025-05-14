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