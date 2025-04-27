#' Rapidly calculate the correlation and the significance of columns in a data.frame with a vector
#'
#' @param df A data.frame
#' @param vec A vector. Should have length equal to the number of rows in df.
#' @param cor_method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor(). Default is "pearson".
#' @param p_adjust_method Method used to adjust p-values. Same as the methods from p.adjust.methods. Default is Benjamini-Hochberg.
#' @return A data.frame with the correlation, p-value and adjusted p-value for each column in df with vec. row.names are the colnames of df
#' @export
rapid_cor_test = function(df, vec, cor_method = "p", p_adjust_method = "BH"){
  
  # Make sure vec is a vector, for example if a column of a data.frame is provided
  vec = as.numeric(vec)
  
  # Check that the length of vec equals the number of rows of df
  if(nrow(df) != length(vec)){
    stop("length of vec doesn't equal the number of rows of df")
  }
  
  # Calculate the number of complete pairs of observations for each column in df with vec
  n = colSums((!is.na(df)) & !is.na(vec))
  
  # Calculate specified correlation values
  cors = as.numeric(cor(df, vec, use = "p", method = cor_method))
  
  # Calculate t-statistics from correlations
  t_stat = cors * sqrt(n-2)/sqrt((1-cors^2))
  
  # Calculate p-values from t-statistics
  p_val =  2*(pt(-abs(t_stat), df = n - 2))
  
  # Calculate q-values from p-values using specified method
  q_val = p.adjust(p = p_val, method = p_adjust_method)
  
  # Return a data.frame with correlations, p-values and q-values
  return(data.frame(correlation = cors, p_value = p_val, q_value = q_val))
  
}


#' Correlation of two vectors using common names
#'
#' Calculates the specified correlation coefficient of two named numeric vectors, matching elements by common names and ignoring elements whose name is not present in both vectors.
#' Can optionally uses a prefix of specified length from the names if only prefixes of names match. Can also calculate the significance of the correlation. 
#'
#' @param x A named numeric vector
#' @param y A named numeric vector
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x} and \code{y}
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor()
#' @param calc_significance A logical value indicating whether to calculate the significance of the correlation using \code{cor.test}. Default is F.
#' @return The specified correlation coefficient of matched elements in \code{x} and \code{y}
#' @export
cor_by_names = function(x, y, prefix = NULL, method = "p", calc_significance = F){
  if(!is.null(prefix)){
    names(x) = substr(names(x), start = 1, stop = prefix)
    names(y) = substr(names(y), start = 1, stop = prefix)}
  common_names = intersect(names(x), names(y))
  
  cor_df = data.frame(x = x[common_names], y = y[common_names])
  cor_df = cor_df[complete.cases(cor_df), ]
  
  if(!calc_significance){
    result = cor(cor_df$x, cor_df$y, method = method, use = "c")
  } else {result = cor.test(cor_df$x, cor_df$y, method = method, use = "c")}
  return(result)
}

#' Calculate partial correlation between two vectors matching elements by names, controlling for another feature which shares at least some of the same names
#'
#' Calculates the partial correlation coefficient of two named numeric vectors while controlling for a third vector,
#' matching elements by common names and ignoring elements whose name is not present in all three vectors.
#' Can optionally uses a prefix of specified length from the names if only prefixes of names match.
#'
#' @param x A numeric vector (optionally named) used to calculate the partial correlation of \code{x} and \code{y} controlling \code{z}
#' @param y A numeric vector (optionally named) used to calculate the partial correlation of \code{x} and \code{y} controlling \code{z}
#' @param z A vector of values (optionally named) used to calculate the partial correlation with \code{x} and \code{y} controlling \code{z}
#' @param use_names A logical value indicating whether elements in vectors should be matched by names 
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x}, \code{y} and \code{z}
#' @param method The correlation method to use. Options are "pearson" and "spearman". Default is "pearson". 
#' @param calc_significance A logical value indicating whether to calculate the significance of the correlation using \code{cor.test}. Default is F.
#' @return The specified partial correlation coefficient of matched elements in \code{x} and \code{y}, controlling for \code{z}
#' @export
partial_cor = function(x, y, z, use_names = F, prefix = NULL, method = "pearson", calc_significance = F){
  
  method = match.arg(arg = method, choices = c("pearson", "spearman"))
  
  if(use_names){if(!is.null(prefix)){
      names(x) = substr(names(x), start = 1, stop = prefix)
      names(y) = substr(names(y), start = 1, stop = prefix)
      names(z) = substr(names(z), start = 1, stop = prefix)}
    common_names = intersect(intersect(names(x), names(y)), names(z))
    x = x[common_names]; y = y[common_names]; z = z[common_names]}
  
  cor_df = data.frame(x = x, y = y, z = z)
  cor_df = cor_df[complete.cases(cor_df), ]
  
  if(method == "spearman"){
    cor_df = lapply(cor_df, rank)
  }
  
  x_without_z = residuals(lm(x ~ z, data = cor_df))
  y_without_z = residuals(lm(y ~ z, data = cor_df))
  if(!calc_significance){result = cor(x_without_z, y_without_z)}
  if(calc_significance){result = cor.test(x_without_z, y_without_z)}
  return(result)
}

#' Correlation values for pairs of features from two different tables
#'
#' Calculate the correlation values for pairs of features from two different tables and return results in a data.frame.
#' Features should be rows and samples should be columns.
#' All pair-wise combinations of features across the two tables are tested by default, though subsets of features to test can be specified from both tables.
#' Samples are matched between tables using column names, with samples not present in both tables ignored.
#' Samples can optionally be matched using a prefix of the sample names of a specified length.
#' Pearson's correlation is calculated by default, though Spearman's correlation can also be used.
#' The significance of the correlation value can also be returned, though this will increase the running time by a factor of ~10.
#' Results can be sorted by absolute correlation value or by significance value (if calc_significance = T) using sort_by_correlation and sort_significance respectively.
#'
#' @param table1 A matrix or data.frame where rows are features and columns are samples. Can also be a named vector.
#' @param table2 A matrix or data.frame where rows are features and columns are samples. Can also be a named vector. 
#' @param table1_features A character vector with a subset of feature names to test from table1
#' @param table2_features A character vector with a subset of feature names to test from table2
#' @param feature_matches An optional data.frame with two columns, the first for features in table1 and the second for features in table2.
#' If provided, only correlations involving pairs of features found in the rows of the data.frame are calculated.
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x} and \code{y}
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor()
#' @param sort_by_correlation A character string indicating whether to sort results by correlation values in decreasing order ("decreasing"), in increasing order ("increasing") or by absolute correlation ("absolute"). Default is not to sort results. 
#' @param calc_significance A logical value indicating whether a p-value should be calculated for the correlation using cor.test(). A q-value will also be returned using FRD correction
#' @return A data.frame with the correlation values for all pairs of features (or just those indicated) from the two input tables
#' @export
cor_tables = function(table1, table2, table1_features = NULL, table2_features = NULL, feature_matches = NULL, prefix = NULL,
  method = "pearson", sort_by_correlation = NULL, calc_significance = F){

    # Check that allowed values for correlation method and sort_by_correlation are supplied
    match.arg(arg = method, choices = c("pearson", "p", "spearman", "s"), several.ok = F)
    match.arg(arg = sort_by_correlation, choices = c("increasing", "decreasing", "absolute"), several.ok = F)

    # If either feature tables are actually vectors, convert them into a data.frame with one row
    if(is.vector(table1)){table1 = data.frame(t(table1), row.names = deparse(substitute(table1)))}
    if(is.vector(table2)){table2 = data.frame(t(table2), row.names = deparse(substitute(table2)))}

    # If features are not specified for table1 or table2, all features are used by default
    if(is.null(table1_features)){table1_features = row.names(table1)}
    if(is.null(table2_features)){table2_features = row.names(table2)}

    # If it is specified to use a prefix of sample names rather than full names, column names are set to the prefixes of specified length
    if(!is.null(prefix)){
      colnames(table1) = substr(colnames(table1), 1, prefix)
      colnames(table2) = substr(colnames(table2), 1, prefix)}

    # Find sample names common to both tables and subset tables so that they contain these samples
    common_names = intersect(colnames(table1), colnames(table2))
    if(length(common_names) == 0){stop("There are no samples in common between table1 and table2")}
    table1 = table1[, common_names]
    table2 = table2[, common_names]
    
    # If feature_matches not provided, calculate all pairwise correlations of features from table1 and table2
    if(is.null(feature_matches)){

      # Subset tables so that they only contain the specified features (if provided) and transpose tables so that features are now columns
      table1 = t(table1[table1_features, , drop = F])
      table2 = t(table2[table2_features, , drop = F])
      
      # Calculate correlation of features and significance using the specified method
      cor_matrix = cor(table1, table2, method = method, use = "p")
      cor_df = data.frame(table1_feature = rep(row.names(cor_matrix), times = ncol(cor_matrix)), 
        table2_feature = rep(colnames(cor_matrix), each = nrow(cor_matrix)), cor = c(cor_matrix))
      if(calc_significance){
        cor_df$p_value = c(apply(table2, 2, function(x)
          apply(table1, 2, function(y) tryCatch(cor.test(x, y, method = method, use = "c")$p.value, error = function(err) NA))))
        cor_df$q_value = p.adjust(cor_df$p_value, method = "fdr")
      }
    } else {

      # If feature_matches provided, calculate correlations for just the indicated pairs of features
      cor_results = apply(feature_matches, 1, function(match)
        tryCatch(cor(unlist(table1[match[[1]], ]), unlist(table2[match[[2]], ]), method = method, use = "p"), error = function(err) NA))
      cor_df = data.frame(table1_feature = feature_matches[,1], table2_feature = feature_matches[,2],
        cor = cor_results)
      if(calc_significance){
        cor_df$p_value = apply(feature_matches, 1, function(match)
          tryCatch(cor.test(unlist(table1[match[[1]], ]), unlist(table2[match[[2]], ]), method = method, use = "p")$p.value, error = function(err) NA))
        cor_df$q_value = p.adjust(cor_df$p_value, method = "fdr")
      }
    }

    # Sort by the absolute correlation value if specified then return results
    if(!is.null(sort_by_correlation)){
      if(sort_by_correlation == "decreasing"){cor_df = dplyr::arrange(cor_df, dplyr::desc(cor))
      } else if(sort_by_correlation == "increasing"){cor_df = dplyr::arrange(cor_df, cor)
      } else if(sort_by_correlation == "absolute"){cor_df = dplyr::arrange(cor_df, dplyr::desc(abs(cor)))}
    }
    return(cor_df)

}

###

#' Partial correlation values for pairs of features from two different tables
#'
#' Calculate the partial correlation values for pairs of features from two different tables while controliing for a specified confounding variable and return results in a data.frame.
#' Features should be rows and samples should be columns.
#' All pair-wise combinations of features across the two tables are tested by default, though subsets of features to test can be specified from both tables.
#' Samples are matched between tables and the confounding variable vector using column names of the tables and names of the confounding vector, 
#' with samples not present in both tables and the confounding vector ignored.
#' Samples can optionally be matched using a prefix of the sample names of a specified length.
#' The significance of the correlation value can also be returned if specified.
#' Results can be sorted by absolute correlation value or by significance value (if calc_significance = T) using sort_by_correlation and sort_significance respectively.


#' @param table1 A matrix or data.frame where rows are features and columns are samples.
#' @param table2 A matrix or data.frame where rows are features and columns are samples.
#' @param confounder_table A matrix or data.frame where rows are features and columns are samples.
#' @param feature_matches A data.frame with three columns, the first for features in table1, the second for features in table2 and the 3rd for confounding features in confounder_table.
#' Partial correlations are calculated using each row of feature_matches between the table1 and table2 features controlling for the confounder_table features. 
#' @param method The type of partial correlation to compute. Options are "pearson" and "spearman". Default is "pearson". 
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x} and \code{y}
#' @param sort_by_correlation A character string indicating whether to sort results by correlation values in decreasing order ("decreasing"), in increasing order ("increasing") or by absolute correlation ("absolute"). Default is not to sort results. 
#' @return A data.frame with the correlation values for all pairs of features (or just those indicated) from the two input tables
#' @export
partial_cor_tables = function(table1, table2, confounder_table, feature_matches, method = "pearson", prefix = NULL,
 sort_by_correlation = NULL){
  
  method = match.arg(arg = method, choices = c("pearson", "spearman"))

  # Check that an allowed value for sort_by_correlation is supplied
  match.arg(arg = sort_by_correlation, choices = c("increasing", "decreasing", "absolute"), several.ok = F)

  # If either feature tables are actually vectors, convert them into a data.frame with one row
  if(is.vector(table1)){table1 = data.frame(t(table1), row.names = deparse(substitute(table1)))}
  if(is.vector(table2)){table2 = data.frame(t(table2), row.names = deparse(substitute(table2)))}
  if(is.vector(confounder_table)){table2 = data.frame(t(confounder_table), row.names = deparse(substitute(confounder_table)))}

  # If it is specified to use a prefix of sample names rather than full names, column names are set to the prefixes of specified length
  if(!is.null(prefix)){
    colnames(table1) = substr(colnames(table1), 1, prefix)
    colnames(table2) = substr(colnames(table2), 1, prefix)
    colnames(confounder_table) = substr(colnames(confounder_table), 1, prefix)}

  # Find sample names common to both tables and subset tables so that they contain these samples
  common_names = intersect(intersect(colnames(table1), colnames(table2)), names(confounder_table))
  table1 = table1[, common_names]
  table2 = table2[, common_names]
  confounder_table = confounder_table[common_names]
  
  # If feature_matches provided, calculate correlations for just the indicated pairs of features
  partial_cor_results = apply(feature_matches, 1, function(match)
    tryCatch({correlateR::partial_cor(x = unlist(table1[match[1], ]), y = unlist(table2[match[2], ]), z = unlist(confounder_table[match[3], ]), use_names = F,
      method = method, calc_significance = calc_significance)}, error = function(err) NA))
  cor_df = data.frame(table1_feature = feature_matches[,1], table2_feature = feature_matches[,2], 
    confounder = feature_matches[,3], cor = partial_cor_results)


  # Sort by the absolute correlation value if specified then return results
  if(!is.null(sort_by_correlation)){
    if(sort_by_correlation == "decreasing"){cor_df = dplyr::arrange(cor_df, dplyr::desc(cor))
    } else if(sort_by_correlation == "increasing"){cor_df = dplyr::arrange(cor_df, cor)
    } else if(sort_by_correlation == "increasing"){cor_df = dplyr::arrange(cor_df, dplyr::desc(abs(cor)))}
  }
  return(cor_df)

}

#' Create a scatter plot using two variables
#'
#' Calculates the specified correlation coefficient of two named numeric vectors, matching elements by common names and ignoring elements whose name is not present in both vectors.
#' Can optionally uses a prefix of specified length from the names if only prefixes of names match.
#'
#' @param x A named numeric vector
#' @param y A named numeric vector
#' @param use_names A logical value indicating whether elements in \code{x} and \code{y} should be matched by their names. Default is F.
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x} and \code{y}
#' @param regression_line A logical value indicating whether to add a regression line to the plot
#' @param correlation_label A character string to precede the correlation value added to the plot. Default is NULL i.e. add no correlation label
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor()
#' @param position Where to place the correlation label. Options are "bottomright", "topright", "bottomleft" or "topleft". Default is "bottom.right"
#' @param label_size Size of the correlation_label. Default is 1
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param title Title of the plot
#' @param alpha Value specifying the alpha
#' @return The created ggplot object
#' @export
plot_features = function(x, y, use_names = F, prefix = NULL, regression_line = F, correlation_label = NULL, method = "p", position = "bottomright",
  plot_title_size = 24, axis_title_size = 20, axis_text_size = 18, label_size = 5, xlab = "feature1", ylab = "feature2", title = "Feature1 Vs Feature2", alpha = 1){
  
  # Check that an allowed value is used for position
  match.arg(arg = position, choices = c("bottomright", "topright", "bottomleft", "topleft"), several.ok = F)
  
# Create a scatterplot of two features, matching them by names (or optionally using a prefix of names)
  if(use_names){
    if(!is.null(prefix)){
      names(x) = substr(names(x), start = 1, stop = prefix)
      names(y) = substr(names(y), start = 1, stop = prefix)
      }
    common_names = intersect(names(x), names(y))
    x = x[common_names]
    y = y[common_names]
  }
  feature_df = data.frame(x = x, y = y)
  plot = ggplot2::ggplot(feature_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(alpha = alpha) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
      ggplot2::scale_x_continuous(expand = c(0, 0.2)) +
      ggplot2::theme_classic() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = plot_title_size), 
        axis.title = element_text(size = axis_title_size), 
        axis.text = element_text(size = axis_text_size)) +
      ggplot2::labs(x = xlab, y = ylab, title = title)
  if(regression_line){plot = plot + ggplot2::geom_smooth(method = "lm", se = F)}
  if(!is.null(correlation_label)){
    plot = plot + ggplot2::annotate("text", label = paste0(correlation_label, round(cor(x, y, method = method, use = "c"), 2)),
      x = ifelse(grepl("right", position), max(x, na.rm = T), min(x, na.rm = T)),
      y = ifelse(grepl("top", position), max(y, na.rm = T), min(y, na.rm = T)),
      hjust = ifelse(grepl("right", position), 1, 0), size = label_size)}
  return(plot)
}

#' Plot distribution of correlation values
#'
#' @param correlation_table A data.frame with a column named cor giving correlation values between different variables as produced by cor_tables, for example. 
#' @param binwidth Width of bins for correlation values in the histogram. Default is 0.1.
#' @param x_axis_breaks A single value indicating how often to label the bins on the X-axis as they move away from 0. Default is to label bins after every change in correlation value of 0.2.
#' @param xlab Label for the x-axis. Default is "Correlation".
#' @param ylab Label for the y-axis. Default is "Proportion"
#' @param title Title of the plot. Default is "Histogram of Correlations"
#' @param annotate_significant_correlations A logical value indicating whether to indicate the number of significant correlations (using q-values) to the plot. Default is FALSE. 
#' If TRUE, correlation_table should contain a column called q_value. 
#' @return The created ggplot object
#' @export
cor_histogram = function(correlation_table, binwidth = 0.1, x_axis_breaks = 0.2, xlab = "Correlation", ylab = "Proportion", title = "Histogram of Correlations", annotate_significant_correlations = F){
  
  # Check that correlation_table has a column called cor
  if(!"cor" %in% names(correlation_table)){stop("correlation_table must have a coilumn called \"cor\"")}
  
  # If annotate_significant_correlations is TRUE, check that there is a column called q_value
  if(annotate_significant_correlations){
    if(!"q_value" %in% names(correlation_table)){stop("Cannot annotate number of significant correlations if correlation_table does not contain q_values")}
  }
  
  # Create histogram
  plot = ggplot(correlation_table, aes(x = cor, y = (..count..)/sum(..count..))) +
    geom_histogram(color = "darkblue", fill = "lightblue", binwidth = binwidth) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), 
    axis.title = element_text(size = 20), axis.text = element_text(size = 18), axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(-1, 1, x_axis_breaks))) +
    labs(x = xlab, y = ylab, title = title)
  
  # Add number of significant correlations if specified
  if(annotate_significant_correlations){
      plot = plot + annotate("text", label = paste(prettyNum(sum(correlation_table$q_value < 0.05, na.rm = T), big.mark = ","),  "q-values < 0.05"), 
                      x = max(correlation_table$cor, na.rm = T), 
                      y = max(table(as.numeric(as.character(cut(correlation_table$cor, breaks = seq(-1.05, 1.05, binwidth), labels = seq(-1, 1, binwidth)))))/sum(!is.na(correlation_table$cor))), 
                      hjust = 1, vjust = 1, size = 8)
  }
  
  
  return(plot)
}

#' @param table1 A matrix or data.frame where rows are features and columns are samples.
#' @param table2 A matrix or data.frame where rows are features and columns are samples.
#' @param method Correlation method. Default is "pearson".
#' @param use See cor for details. Default is "p"
#' @param n_permutations Number of permutations for permutation test. Default is 10,000.
#' @return A data.frame with the correlation values for all pairs of features (or just those indicated) from the two input tables
#' @export
cor_permutation_test = function(table1, table2, use = "p", method = "p", n_permutations = 10000){
  
  # Convert table1 or table2 into matrices if they are vectors
  if(is.vector(table1)){
    table1 = matrix(table1, ncol = 1)
    colnames(table1) = deparse(substitute(table1))
    }
  if(is.vector(table1)){
    table1 = matrix(table2, ncol = 1)
    colnames(table2) = deparse(substitute(table2))
  }
  
  # Check that table1 and table2 have colnames
  if(is.null(colnames(table1)) | is.null(colnames(table2))){stop("table1 and table2 must have column names")}
  
  # Check that table1 and table2 have the same number of rows
  if(nrow(table1) != nrow(table2)){stop("table1 and table2 have different number of rows")}
  
  # Calculate correlations between table1 and table2 and turn results into a data.frame
  cor_results = cor(table1, table2, use = use, method = method, na.rm = T)
  cor_df = data.frame(table1_feature = rep(row.names(cor_results), times = ncol(cor_results)), 
        table2_feature = rep(colnames(cor_results), each = nrow(cor_results)), cor = c(cor_results))
  
  # Create a list of tables with the shuffled values for each clumn in table2
  shuffled_tables2 = setNames(lapply(1:ncol(table2), function(x) 
    matrix(unlist(lapply(1:n_permutations, function(y) sample(table2[, x]))), ncol = n_permutations)
  ), colnames(table2))
  
  # Calculate correlations with shuffled columns and combine into a single matrix
  shuffled_table2_combined = do.call(rbind, lapply(shuffled_table2_list, function(x) cor(table1, x, , use = use, method = method)))

  # Calculate p-value, the number of time shuffled correlations are equal or greater than the observed correlation
  cor_df$p_value = (rowSums(abs(shuffled_table2_combined) >= (abs(cor_df[["cor"]]) -0.001)) + 1)/(n_permutations + 1)
  
  return(cor_df)
  
}

#' Calculate correlation values for random pairs of features from one table or with each feature coming from two different tables
#'
#' @param table1 A matrix or data.frame where rows are features and columns are samples. Can also be a named vector.
#' @param table2 An optional  matrix or data.frame where rows are features and columns are samples. Can also be a named vector. 
#' If table2 not provided, calculates correlation values for random pairs of features from table1.
#' @param n_pairs Number of random pairs to calculate correlation values for 
#' @param prefix An integer specifying the length of the prefixes to use for names of \code{x} and \code{y}
#' @param method A character string indicating which correlation coefficient is to be computed. Identical to methods from cor()
#' @param sort_by_correlation A character string indicating whether to sort results by correlation values in decreasing order ("decreasing"), in increasing order ("increasing") or by absolute correlation ("absolute"). Default is not to sort results. 
#' @param calc_significance A logical value indicating whether a p-value should be calculated for the correlation using cor.test(). A q-value will also be returned using FRD correction
#' @return A data.frame with the correlation values for all pairs of features (or just those indicated) from the two input tables
#' @export
random_cor_tables = function(table1, table2 = NULL, n_pairs, prefix = NULL, method = "pearson", sort_by_correlation = NULL, calc_significance = F){

    # Check that allowed values for correlation method and sort_by_correlation are supplied
    method = match.arg(arg = method, choices = c("pearson", "spearman"), several.ok = F)
    match.arg(arg = sort_by_correlation, choices = c("increasing", "decreasing", "absolute"), several.ok = F)
    
    # If table2 is not provided, set table1 as table2 and note that the tables are the same
    if(is.null(table2)){
      table2 = table1
      ntables = 1
    } else {
      ntables = 2
    }

    # If it is specified to use a prefix of sample names rather than full names, column names are set to the prefixes of specified length
    if(!is.null(prefix)){
      colnames(table1) = substr(colnames(table1), 1, prefix)
      colnames(table2) = substr(colnames(table2), 1, prefix)}

    # Find sample names common to both tables and subset tables so that they contain these samples
    common_names = intersect(colnames(table1), colnames(table2))
    if(length(common_names) == 0){stop("There are no samples in common between table1 and table2")}
    table1 = table1[, common_names]
    table2 = table2[, common_names]
    
    # Create a data.frame with all pairs of features from two tables 
    feature_matches = data.frame(
      table1_features = sample(row.names(table1), size = n_pairs, replace = T),
      table2_features = sample(row.names(table2), size = n_pairs, replace = T)
      )
    
    # Make sure no feature pairs duplicated
    feature_matches = unique(feature_matches)
    
    # Remove rows where both features are the same if just using one table
    if(ntables == 1){
      feature_matches = dplyr::filter(feature_matches, table1_features != table2_features)
    }

    # Calculate correlation for random feature pairs
    cor_results = apply(feature_matches, 1, function(match)
      tryCatch(cor(unlist(table1[match[[1]], ]), unlist(table2[match[[2]], ]), method = method, use = "p"), error = function(err) NA))
    cor_df = data.frame(table1_feature = feature_matches[,1], table2_feature = feature_matches[,2],
      cor = cor_results)
    if(calc_significance){
      cor_df$p_value = apply(feature_matches, 1, function(match)
        tryCatch(cor.test(unlist(table1[match[[1]], ]), unlist(table2[match[[2]], ]), method = method, use = "p")$p.value, error = function(err) NA))
      cor_df$q_value = p.adjust(cor_df$p_value, method = "fdr")
    }

    # Sort by the absolute correlation value if specified then return results
    if(!is.null(sort_by_correlation)){
      if(sort_by_correlation == "decreasing"){cor_df = dplyr::arrange(cor_df, dplyr::desc(cor))
      } else if(sort_by_correlation == "increasing"){cor_df = dplyr::arrange(cor_df, cor)
      } else if(sort_by_correlation == "absolute"){cor_df = dplyr::arrange(cor_df, dplyr::desc(abs(cor)))}
    }
    return(cor_df)

}
