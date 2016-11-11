#' @export
fromMpMap <- function(mpcross, selfing = "infinite", fixCodingErrors = FALSE)
{
  isOldMpMapMpcrossArgument(mpcross)
  
  oldPedigree <- mpcross$pedigree
  pedigreeLineNames <- rownames(oldPedigree)
  if(is.null(pedigreeLineNames) || length(unique(pedigreeLineNames)) != length(pedigreeLineNames))
  {
    stop("Pedigree of the input object must have unique row names")
  }
  #Attempt the reordering call, which requires building the package with the boost library
  reorderedPedigree <- reorderPedigree(lineNames = pedigreeLineNames, mother = as.integer(oldPedigree[, "Female"]), father = as.integer(oldPedigree[, "Male"]), selfing = selfing, warnImproperFunnels = TRUE)
  if(!is.null(reorderedPedigree))
  {
    newPedigree <- reorderedPedigree
  }
  else
  {
    newPedigree <- pedigree(lineNames = pedigreeLineNames, mother = oldPedigree[, "Female"], father = oldPedigree[, "Male"], selfing = selfing, warnImproperFunnels = TRUE)
  }

  finalsMarkerNames <- colnames(mpcross$finals)
  foundersMarkerNames <- colnames(mpcross$finals)
  if(!all.equal(finalsMarkerNames, foundersMarkerNames))
  {
	  stop("Founder and final marker names must be identical")
  }

  newHetDataList <- lapply(as.list(1:ncol(mpcross$founders)), function(x)
    {
      uniqueAlleles <- unique(mpcross$founders[, x])
      retVal <- cbind(uniqueAlleles, uniqueAlleles, uniqueAlleles)
      colnames(retVal) <- NULL
      return(retVal)
    })
  names(newHetDataList) <- foundersMarkerNames
  newHetData <- new("hetData", newHetDataList)

  if(fixCodingErrors)
  {
    codingErrors <- .Call("listCodingErrors", mpcross$founders, mpcross$finals, newHetData)
    uniqueMarkers <- unique(codingErrors$finals[,"Column"])
    #The results of listCodingErrors are zero-based indexing
    mpcross$finals[, uniqueMarkers+1] <- NA
    warning(paste0("Removing data for ", length(uniqueMarkers), " markers, because fixCodingErrors = TRUE was specified"))
  }
  geneticData <- new("geneticData", finals = mpcross$finals, founders = mpcross$founders, pedigree = newPedigree, hetData = newHetData)
  geneticDataList <- new("geneticDataList", list(geneticData))
  if("map" %in% names(mpcross))
  {
    return(new("mpcrossMapped", geneticData = geneticDataList, map = mpcross$map))
  }
  else
  {
    return(new("mpcross", geneticData = geneticDataList))
  }
}
setMethod(f = "+", signature = c("mpcross", "mpcross"), definition = function(e1, e2)
{
  if(pryr::address(e1) == pryr::address(e2))
  {
    stop("Cannot combine an object with itself")
  }
  allMarkers <- sort(unique(c(markers(e1), markers(e2))))
  #If e1 and e2 have exactly the same markers, this is trivial. 
  if(nMarkers(e1) == nMarkers(e1) && all(markers(e1) == markers(e2)))
  {
    return(new("mpcross", geneticData = new("geneticDataList", c(e1@geneticData, e2@geneticData))))
  }
  #If all the markers are already in e1, then use that order and only expand e2
  else if(nMarkers(e1) == length(allMarkers) && all(allMarkers == sort(markers(e1))))
  {
    e2ExpandedGeneticData <- expand(e2, newMarkers = allMarkers)
    return(new("mpcross", geneticData = new("geneticDataList", c(e1@geneticData, e2ExpandedGeneticData@geneticData))))
  }
  #If all the markers are already in e2, then use that order and only expand e1
  else if(nMarkers(e2) == length(allMarkers) && all(allMarkers == sort(markers(e2))))
  {
    e1ExpandedGeneticData <- expand(e1, newMarkers = allMarkers)
    return(new("mpcross", geneticData = new("geneticDataList", c(e1ExpandedGeneticData@geneticData, e2@geneticData))))
  }
  else
  {
    #Put all the markers from e1 first
    allMarkers <- c(markers(e1), setdiff(allMarkers, markers(e1)))
    e1ExpandedGeneticData <- expand(e1, newMarkers = allMarkers)
    e2ExpandedGeneticData <- expand(e2, newMarkers = allMarkers)
    return(new("mpcross", geneticData = new("geneticDataList", c(e1ExpandedGeneticData@geneticData, e2ExpandedGeneticData@geneticData))))
  }
})
setMethod(f = "+", signature = c("mpcrossRF", "mpcrossRF"), definition = function(e1, e2)
{
  if(length(e1@rf@theta@levels) != length(e2@rf@theta@levels) || any(e1@rf@theta@levels != e2@rf@theta@levels))
  {
    stop("Different recombination values were used for numerical maximum likelihood in two objects")
  }
  levels <- e1@rf@theta@levels
  keepLod <- !is.null(e1@rf@lod) && !is.null(e2@rf@lod)
  keepLkhd <- !is.null(e1@rf@lkhd) && !is.null(e2@rf@lkhd)

  newLod <- newLkhd <- NULL

  combined <- as(e1, "mpcross") + as(e2, "mpcross")
  lineWeights <- lapply(as.list(nLines(combined)), function(x) rep(1, x))
  marker1Indices <- match(markers(e1), markers(combined))
  marker2Indices <- match(markers(e2), markers(combined))
  intersectionIndices <- intersect(marker1Indices, marker2Indices)

  newGbLimit <- min(e1@rf@gbLimit, e2@rf@gbLimit)

  #Make a new rawSymmetricMatrix
  warning("Implicitly setting lineWeights parameter to 1")
  dataLengths <- nMarkers(combined) *(nMarkers(combined)+1)/2
  newTheta <- new("rawSymmetricMatrix", data = raw(dataLengths), levels = levels, markers = markers(combined))
  #Copy over all the existing data
  .Call("assignRawSymmetricMatrixDiagonal", newTheta, marker1Indices, e1@rf@theta@data, PACKAGE = "mpMap2")
  .Call("assignRawSymmetricMatrixDiagonal", newTheta, marker2Indices, e2@rf@theta@data, PACKAGE = "mpMap2")
  if(keepLod)
  {
    newLod <- new("dspMatrix", x = vector(mode="numeric", length = dataLengths), Dim = c(nMarkers(combined), nMarkers(combined)))
    colnames(newLod) <- rownames(newLod) <- markers(combined)
    newLod[marker1Indices, marker1Indices] <- e1@rf@lod
    newLod[marker2Indices, marker2Indices] <- e2@rf@lod
  }
  if(keepLkhd)
  {
    newLkhd <- new("dspMatrix", x = vector(mode="numeric", length = dataLengths), Dim = c(nMarkers(combined), nMarkers(combined)))
    colnames(newLkhd) <- rownames(newLkhd) <- markers(combined)
    newLkhd[marker1Indices, marker1Indices] <- e1@rf@lkhd
    newLkhd[marker2Indices, marker2Indices] <- e2@rf@lkhd
  }
  complementIntersectionIndices <- setdiff(1:nMarkers(combined), intersectionIndices)
  if(length(intersectionIndices) > 0)
  {
    reEstimatedPart1 <- estimateRFInternal(object = combined, recombValues = levels, lineWeights = lineWeights, keepLod = keepLod, keepLkhd = keepLkhd, markerRows = 1:nMarkers(combined), markerColumns = intersectionIndices, gbLimit = newGbLimit, verbose = list(verbose = FALSE, progressStyle = 1L))
    .Call("assignRawSymmetricMatrixFromEstimateRF", newTheta, 1:nMarkers(combined), intersectionIndices, reEstimatedPart1$theta)

    if(length(complementIntersectionIndices) > 0)
    {
      reEstimatedPart2 <- estimateRFInternal(object = combined, recombValues = levels, lineWeights = lineWeights, keepLod = keepLod, keepLkhd = keepLkhd, markerRows = intersectionIndices, markerColumns = complementIntersectionIndices, gbLimit = newGbLimit, verbose = list(verbose = FALSE, progressStyle = 1L))
      .Call("assignRawSymmetricMatrixFromEstimateRF", newTheta, intersectionIndices, complementIntersectionIndices, reEstimatedPart2$theta)
    }

    if(keepLod)
    {
      .Call("assignDspMatrixFromEstimateRF", newLod, 1:nMarkers(combined), intersectionIndices, reEstimatedPart1$lod, PACKAGE="mpMap2")
      if(length(complementIntersectionIndices) > 0) .Call("assignDspMatrixFromEstimateRF", newLod, intersectionIndices, complementIntersectionIndices, reEstimatedPart2$lod, PACKAGE="mpMap2")
    }
    if(keepLkhd)
    {
      .Call("assignDspMatrixFromEstimateRF", newLkhd, 1:nMarkers(combined), intersectionIndices, reEstimatedPart1$lkhd, PACKAGE="mpMap2")
      if(length(complementIntersectionIndices) > 0) .Call("assignDspMatrixFromEstimateRF", newLkhd, intersectionIndices, complementIntersectionIndices, reEstimatedPart2$lkhd, PACKAGE="mpMap2")
    }
  }
  rectangularRows <- setdiff(marker1Indices, intersectionIndices)
  rectangularColumns <- setdiff(marker2Indices, intersectionIndices)
  if(length(rectangularRows) > 0 && length(rectangularColumns) > 0)
  {
    rectangularPart <- estimateRFInternal(object = combined, recombValues = levels, lineWeights = lineWeights, keepLod = keepLod, keepLkhd = keepLkhd, markerRows = rectangularRows, markerColumns = rectangularColumns, gbLimit = newGbLimit, verbose = list(verbose = FALSE, progressStyle = 1L))
    .Call("assignRawSymmetricMatrixFromEstimateRF", newTheta, rectangularRows, rectangularColumns, rectangularPart$theta)
    if(keepLod)
    {
      .Call("assignDspMatrixFromEstimateRF", newLod, rectangularRows, rectangularColumns, rectangularPart$lod, PACKAGE="mpMap2")
    }
    if(keepLkhd)
    {
      .Call("assignDspMatrixFromEstimateRF", newLkhd, rectangularRows, rectangularColumns, rectangularPart$lkhd, PACKAGE="mpMap2")
    }
  }
  newRF <- new("rf", theta = newTheta, lod = newLod, lkhd = newLkhd, gbLimit = newGbLimit)
  return(new("mpcrossRF", combined, rf = newRF))
})
#Generally we drop the RF data, unless the sets of markers are disjoint
setMethod(f = "+", signature = c("mpcrossRF", "mpcross"), definition = function(e1, e2)
{
  #If the sets of markers are disjoint, then we re-run the RF computation (the alternative is to drop the existing RF computations which seems computationally wasteful)
  if(length(intersect(markers(e1), markers(e2))) == 0)
  {
    combined <- as(e1, "mpcross") + e2
    warning("Implicitly setting lineWeights parameter to 1")
    keepLod <- !is.null(e1@rf@lod)
    keepLkhd <- !is.null(e1@rf@lkhd)
    marker1Range <- 1:nMarkers(e1)
    marker2Indices <- match(markers(e2), markers(combined))
    if(any(diff(sort(marker2Indices))) != 1)
    {
      stop("Internal error: Markers should have been combined as two blocks")
    }
    marker2Range <- range(match(markers(e2), markers(combined)))
    extraRFData <- estimateRFInternal(object = combined, recombValues = e1@rf@r, lineWeights = rep(1, nLines(e2)), markerRows = marker1Range, markerColumns = marker2Range, keepLod = keepLod, keepLkhd = keepLkhd, gbLimit = e1@rf@gbLimit, verbose = list(verbose = FALSE, progressStyle = 1L))
    stop("Need to check this section")
  }
  #If the markers are all the same, keep them in the same order
  if(nMarkers(e1) == nMarkers(e2) && all(markers(e1) == markers(e2)))
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
