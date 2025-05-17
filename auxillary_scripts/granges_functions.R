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
