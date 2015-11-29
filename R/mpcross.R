fromMpMap <- function(mpcross)
{
  isOldMpMapMpcrossArgument(mpcross)
}
setMethod(f = "+", signature = c("mpcross", "mpcross"), definition = function(e1, e2)
{
  allMarkers <- unique(c(markers(e1), markers(e2)))
  #Put all the markers from e1 first
  allMarkers <- c(markers(e1), setdiff(allMarkers, markers(e1)))
  e1ExpandedGeneticData <- expand(e1, newMarkers = allMarkers)
  e2ExpandedGeneticData <- expand(e2, newMarkers = allMarkers)
  return(new("mpcross", geneticData = c(e1ExpandedGeneticData@geneticData, e2ExpandedGeneticData@geneticData))) 
})
setMethod(f = "+", signature = c("mpcrossRF", "mpcrossRF"), definition = function(e1, e2)
{
  if(length(e1@rf@theta@levels) != length(e2@rf@theta@levels) || any(e1@rf@theta@levels != e2@rf@theta@levels))
  {
    stop("Different recombination values were used for numerical maximum likelihood in two objects")
  }
  existingLod <- !is.null(e1@rf@lod) || !is.null(e2@rf@lod)
  existingLkhd <- !is.null(e1@rf@lkhd) || !is.null(e2@rf@lkhd)

  newLod <- newLkhd <- NULL

  combined <- as(e1, "mpcross") + as(e2, "mpcross")
  marker1Indices <- match(markers(e1), markers(combined))
  marker2Indices <- match(markers(e2), markers(combined))

  if(any(diff(sort(marker1Indices)) != 1))
  {
    stop("Internal error: Markers from object 1 should be in one consecutive block")
  }
  #If the markers are disjoint we keep the RF data - The resulting matrix is block-diagonal
  if(length(intersect(markers(e1), markers(e2))) == 0)
  {
    if(any(diff(sort(marker2Indices)) != 1))
    {
      stop("Internal error: Markers from object 2 should be in one consecutive block")
    }

    newTheta <- matrix(NA, nMarkers(combined), nMarkers(combined))
    colnames(newTheta) <- rownames(newTheta) <- markers(combined)
    
    if(existingLod) 
    {
      newLod <- newTheta
      newLod[marker1Indices, marker1Indices] <- e1@rf@lod
      newLod[marker2Indices, marker2Indices] <- e2@rf@lod
    }
    if(existingLkhd)
    {
      newLkhd <- newTheta
      newLkhd[marker1Indices, marker1Indices] <- e1@rf@lkhd
      newLkhd[marker2Indices, marker2Indices] <- e2@rf@lkhd
    }

    newTheta[marker1Indices, marker1Indices] <- e1@rf@theta
    newTheta[marker2Indices, marker2Indices] <- e2@rf@theta

    newRF <- new("rf", lod = newLod, lkhd = newLkhd, theta = newTheta)
    return(new("mpcrossRF", combined, rf = newRF))
  }
  #If markers for one are strictly contained within the other
  else if(all(markers(e1)%in% markers(e2)) || all(markers(e2) %in% markers(e1)))
  {
    warning("Implicitly setting lineWeights parameter to 1")
    #We need both otherings, so combine the objects in the other way
    combinedOther <- as(e2, "mpcross") + as(e1, "mpcross")
    marker1IndicesOther <- match(markers(e1), markers(combinedOther))
    marker2IndicesOther <- match(markers(e2), markers(combinedOther))
   
    if(all(markers(e1)%in% markers(e2)))
    {
      newTheta <- e2@rf@theta

      marker1Range <- range(marker1Indices)
      marker2Range <- c(1, nMarkers(combined))
      extraRFData <- estimateRFInternal(object = combined, recombValues = e1@rf@r, lineWeights = lapply(as.list(nLines(combined)), function(x) rep(1, x)), marker1Range = marker1Range, marker2Range = marker2Range, keepLod = existingLod, keepLkhd = existingLkhd)
      indicesX <- match(rownames(extraRFData$theta), rownames(newTheta))
      indicesY <- match(colnames(extraRFData$theta), colnames(newTheta))
      
      newTheta[indicesX, indicesY] <- extraRFData$theta
      newTheta[indicesY, indicesX] <- t(extraRFData$theta)

      if(existingLkhd)
      {
        newLkhd <- e2@rf@lkhd
        newLkhd[indicesX, indicesY] <- extraRFData$lkhd
        newLkhd[indicesY, indicesX] <- t(extraRFData$lkhd)
      }
      if(existingLod)
      {
        newLod <- e2@rf@lod
        newLod[indicesX, indicesY] <- extraRFData$lod
        newLod[indicesY, indicesX] <- t(extraRFData$lod)
      }
      newRF <- new("rf", theta = newTheta, lkhd = newLkhd, lod = newLod)
      return(new("mpcrossRF", combinedOther, rf = newRF))
    }
    else
    {
      newTheta <- e1@rf@theta
      if(any(diff(sort(marker2IndicesOther)) != 1))
      {
        stop("Internal error: Markers from object 2 should be in one consecutive block")
      }
      marker1Range <- c(1, nMarkers(combinedOther))
      marker2Range <- range(marker2IndicesOther)
      extraRFData <- estimateRFInternal(object = combinedOther, recombValues = e1@rf@r, lineWeights = lapply(as.list(nLines(combined)), function(x) rep(1, x)), marker1Range = marker1Range, marker2Range = marker2Range, keepLod = existingLod, keepLkhd = existingLkhd)
      indicesX <- match(rownames(extraRFData$theta), rownames(newTheta))
      indicesY <- match(colnames(extraRFData$theta), colnames(newTheta))

      newTheta[indicesX, indicesY] <- extraRFData$theta
      newTheta[indicesY, indicesX] <- t(extraRFData$theta)

      if(existingLkhd)
      {
        newLkhd <- e1@rf@lkhd
        newLkhd[indicesX, indicesY] <- extraRFData$lkhd
        newLkhd[indicesY, indicesX] <- t(extraRFData$lkhd)
      }
      if(existingLod)
      {
        newLod <- e1@rf@lod
        newLod[indicesX, indicesY] <- extraRFData$lod
        newLod[indicesY, indicesX] <- t(extraRFData$lod)
      }
      newRF <- new("rf", theta = newTheta, lkhd = newLkhd, lod = newLod)
      return(new("mpcrossRF", combined, rf = newRF))
    }
  }
  else
  {
    return(combined)
  }
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
    browser()
    extraRFData <- estimateRFInternal(object = combined, recombValues = e1@rf@r, lineWeights = rep(1, nLines(e2)), marker1Range = marker1Range, marker2Range = marker2Range, keepLod = keepLod, keepLkhd = keepLkhd)
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
