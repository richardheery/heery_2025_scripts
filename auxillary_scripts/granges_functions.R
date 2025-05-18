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
  relative_start = genomeTools::stranded_distance(resize(gr, 1, fix = "start"), reference_positions)
  relative_end = genomeTools::stranded_distance(resize(gr, 1, fix = "end"), reference_positions)
  
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
