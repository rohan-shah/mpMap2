fromMpMap <- function(mpcross)
{
  isOldMpMapMpcrossArgument(mpcross)
}
setMethod(f = "+", signature = c("mpcross", "mpcross"), definition = function(e1, e2)
{
  if(!is.null(e1@map) || !is.null(e2@map))
  {
    warning("Combining mpcross objects which had pre-existing maps. These maps have been removed. ")
    e1@map <- NULL
    e2@map <- NULL
  }
  allMarkers <- unique(markers(e1), markers(e2))
  e1ExpandedGeneticData <- expand(e1, newMarkers = allMarkers)
  e2ExpandedGeneticData <- expand(e2, newMarkers = allMarkers)
  return(new("mpcross", geneticData = c(e1ExpandedGeneticData@geneticData, e2ExpandedGeneticData@geneticData)))
  
})
#Generally we drop the RF data, unless the sets of markers are disjoint
setMethod(f = "+", signature = c("mpcrossRF", "mpcross"), definition = function(e1, e2)
{
  #If the sets of markers are disjoint, then we re-run the RF computation (the alternative is to drop the existing RF computations which seems computationall wasteful)
  if(length(intersect(markers(e1), markers(e2))) == 0)
  {
    combined <- as(e1, "mpcross") + e2
    warning("Implicitly setting lineWeights parameter to 1")
    keepLod <- !is.null(e1@rf@lod)
    keepLkhd <- !is.null(e1@rf@lkhd)
    extraRFData <- estimateRFInternal(object = combined, recombValues = e1@rf@r, lineWeights = rep(1, nLines(e2)), marker1Range = marker1Range, marker2Range = marker2Range, keepLod = keepLod, keepLkhd = keepLkhd)
    stop("Need to check this section")
  }
  #If the markers are all the same, keep them in the same order
  if(length(markers(e1)) == length(markers(e2)) && all(markers(e1) == markers(e2)))
  {
    allMarkers <- markers(e1)
  }
  else
  {
    #Otherwise I'm not sure how these will be ordered
    allMarkers <- unique(c(markers(e1), markers(e2)))
  }
  e1ExpandedGeneticData <- expand(e1, newMarkers = allMarkers)
  e2ExpandedGeneticData <- expand(e2, newMarkers = allMarkers)
  return(new("mpcross", geneticData = c(e1ExpandedGeneticData@geneticData, e2ExpandedGeneticData@geneticData)))
})