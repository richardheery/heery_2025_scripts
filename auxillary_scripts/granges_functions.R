count_covered_bases = function(gr){
  
  return(sum(width(reduce(gr, ignore.strand = T))))

}


calculate_regions_intersections <- function(gr1, gr2, ignore.strand = TRUE, overlap_measure = "absolute"){
  
  # Check that inputs have the correct data type
  stopifnot(is(gr1, "GRanges"), is(gr2, "GRanges"), 
    S4Vectors::isTRUEorFALSE(ignore.strand), is(overlap_measure, "character"))
  
  # Check allowed value provided for overlap_measure
  match.arg(overlap_measure, c("absolute", "proportion", "jaccard"))
  
  # Create GRanges with the intersection and union of gr1 and gr2
  intersection <- GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
  union <- c(gr1, gr2)
  
  # Calculate proportion, Jaccard index or absolute overlap depending on overlap_measure
  if(overlap_measure == "proportion"){
    return(count_covered_bases(intersection)/count_covered_bases(gr1))
  } else if(overlap_measure == "jaccard"){
    return(count_covered_bases(intersection)/count_covered_bases(union))
  } else {
    return(count_covered_bases(intersection))
  }

}

#' Find locations of ranges relative to a reference positions
#'
#' @param gr A GRanges object
#' @param reference_positions A GRanges object. Each range should have width 1. Upstream and downstream refer to reference_positions
#' @param gr_category A category that ranges belong to
#' @return A GRanges object
#' @export
relative_ranges = function(gr, reference_positions, gr_category = NULL){
  
  # Check that all reference positions have width 1
  if(!all(width(reference_positions) == 1)){stop("All reference_positions should have width 1")}

  # Get distances from start and end of ranges in gr from reference_positions
  relative_start = stranded_distance(resize(gr, 1, fix = "start"), reference_positions)
  relative_end = stranded_distance(resize(gr, 1, fix = "end"), reference_positions)
  
  # Create an IRanges with the relative distances
  relative_iranges = IRanges::IRanges(pmin(relative_start, relative_end), pmax(relative_start, relative_end), category = gr_category)
  
  # Convert IRanges to GRanges with "relative" as seqnames
  relative_granges = GenomicRanges::GRanges(seqnames = "relative", ranges = relative_iranges)
  
  # Add metadata from gr to relative_granges
  mcols(relative_granges) = mcols(gr)
  
  # Return relative_granges
  return(relative_granges)
  
}

#' Create relative bins and count the overlaps of relative ranges with relative bins
#'
#' @param relative_ranges A GRanges object returned by relative_ranges
#' @param bin_start Start of bins
#' @param bin_end End of bins
#' @param bin_step Size of bins
#' @param category A category associated with relative_ranges
#' @return A data.frame with the the number of ranges of each category in each bin
#' @export
bin_relative_ranges = function(relative_ranges, bin_start, bin_end, bin_step, category){
  
  bin_ranges = IRanges(
    start = seq(bin_start, bin_end - bin_step, bin_step), 
    end = seq(bin_start + bin_step, bin_end, bin_step) - 1
    )
  
  if(!is.null(category)){
    ranges_list = 
    split(relative_ranges, category)
  } else {ranges_list = list(relative_ranges)}
  
  overlap_df = data.frame(lapply(ranges_list, function(x)
    setNames(countOverlaps(bin_ranges, ranges(x)), start(bin_ranges) + width(bin_ranges)/2)))
  
  overlap_df = tibble::rownames_to_column(overlap_df, "bin_center")
  overlap_df$bin_center = as.numeric(overlap_df$bin_center)
  
  return(overlap_df)

}

#' Calculate the proportion of each region in a query set overlapping a subject set
#'
#' @param query_gr A GRanges object. The proportion of each region overlapping regions from subject_gr will be calculated. 
#' @param subject_gr A GRanges object
#' @param ignore.strand A logical value indicating whether to ignore strand when calculating overlaps. Default is TRUE. 
#' @return A matrix with the size of the overlaps between all pairs of GRanges in grl
#' @export
individual_proportion_overlaps = function(query_gr, subject_gr, ignore.strand = T){
  
  # Set names for query_gr
  names(query_gr) = paste0("region_", seq_along(query_gr))
  
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

#' Calculate the number of ranges of each GRanges object in a list overlapping the others with a minimum overlap proportion threshold
#'
#' @param grl A list of GRanges objects
#' @param ignore.strand A logical value indicating whether strand should be ignored when calculating overlaps. Default is TRUE.
#' @param overlap_threshold Minimum overlap proprtion of the query ranges with the subject to consider them overlapping. Default is 0.5. 
calculate_regions_overlap_list = function(grl, ignore.strand = T, overlap_threshold = 0.5){
  
  overlap_matrix = matrix(NA, nrow = length(grl), ncol = length(grl), dimnames = list(names(grl), names(grl)))
  
  for(i in names(grl)){
    for(j in names(grl)){
      overlap_matrix[i, j] = 
        sum(individual_proportion_overlaps(grl[[i]], grl[[j]], ignore.strand = ignore.strand)$overlap > overlap_threshold)/length(grl[[i]])
    }
  }
  
  return(overlap_matrix)
  
}
