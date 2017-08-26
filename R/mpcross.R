#' @include mpcross-class.R
#' @include geneticData-class.R
#' @export
fromMpMap <- function(mpcross, selfing = "infinite", fixCodingErrors = FALSE)
{
  isOldMpMapMpcrossArgument(mpcross)
  
  oldPedigree <- mpcross$pedigree
  pedigreeLineNames <- rownames(oldPedigree)

  #The row names of the founders may not be lines that are named in the pedigree. In that case, rename them to follow the pedigree and issue a warning
  if(any(!(rownames(mpcross$founders) %in% pedigreeLineNames)))
  {
    warning("The row names of object$founders were not named in the pedigree. These row names are being changed to match those given at the start of the pedigree")
    rownames(mpcross$founders) <- pedigreeLineNames[1:nrow(mpcross$founders)]
  }

  if(!isTRUE(all.equal(rownames(mpcross$finals), pedigreeLineNames[mpcross$id])))
  {
    warning("The row names of object$finals may have been incorrect. These have been changed to match the row names of the pedigree and object$id")
    rownames(mpcross$finals) <- pedigreeLineNames[mpcross$id]
  }

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
  codingErrors <- listCodingErrors(founders = mpcross$founders, finals = mpcross$finals, hetData = newHetData)
  if(fixCodingErrors)
  {
    uniqueMarkers <- unique(codingErrors$finals[,"Column"])
    if(length(uniqueMarkers) > 0)
    {
      mpcross$finals[, uniqueMarkers] <- NA
      warning(paste0("Removing data for ", length(uniqueMarkers), " markers, because fixCodingErrors = TRUE was specified, and these markers had invalid alleles. For less aggressive removal, use listCodingErrors"))
    }
  }
  if(length(codingErrors$null))
  {
    newHetData[codingErrors$null] <- list(matrix(0L, 0, 3))
    mpcross$finals[,codingErrors$null] <- NA
    warning(paste0("Removing data for ", length(codingErrors$null), " markers, because these markers have NA founder alleles"))
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
#' @export
mpcross <- function(founders, finals, pedigree, hetData = infiniteSelfing, fixCodingErrors = FALSE)
{
	if(!isS4(pedigree) || !inherits(pedigree, "pedigree"))
	{
		stop("Input pedigree must be an S4 object of class peigree")
	}
	if(is.character(founders) || is.null(dim(founders)) || length(dim(founders)) != 2)
	{
		stop("Input founders must be a numeric matrix")
	}
	if(is.data.frame(founders)) founders <- as.matrix(founders)
	mode(founders) <- "integer"
	if(is.character(finals) || is.null(dim(finals)) || length(dim(finals)) != 2)
	{
		stop("Input finals must be a numeric matrix")
	}
	if(is.data.frame(finals)) finals <- as.matrix(finals)
	if(is.null(rownames(finals)) || is.null(colnames(finals)) || is.null(rownames(founders)) || is.null(colnames(founders)))
	{
		stop("Inputs founders and finals must have row and column names")
	}
	mode(finals) <- "integer"
	if(is.function(hetData))
	{
		hetData <- hetData(founders, finals, pedigree)
	}
	else if(!isS4(hetData) || !inherits(hetData, "hetData"))
	{
		stop("Input hetData must be an object of class hetData")
	}
	if(ncol(founders) != ncol(finals) || ncol(finals) != length(hetData))
	{
		stop("Inputs hetData, founders and finals must have the same number of markers")
	}
	sortedFounderMarkers <- sort(colnames(founders))
	sortedFinalMarkers <- sort(colnames(finals))
	sortedHetDataMarkers <- sort(names(hetData))
	if(any(sortedFounderMarkers != sortedFinalMarkers) || any(sortedFinalMarkers != sortedHetDataMarkers))
	{
		stop("Inputs founders, finals and hetData must have the same markers")
	}
	#Standardise marker order, if required
	if(any(colnames(founders) != colnames(finals)) || any(colnames(finals) != names(hetData)))
	{
		finals <- finals[,sortedFounderMarkers]
		hetData <- hetData[sortedFounderMarkers]
	}
	codingErrors <- listCodingErrors(founders = founders, finals = finals, hetData = hetData)
	if(length(codingErrors$null))
	{
		hetData[codingErrors$null] <- list(matrix(0L, 0, 3))
		finals[,codingErrors$null] <- NA
		warning(paste0("Removing data for ", length(codingErrors$null), " markers, because these markers have NA founder alleles"))
	}

	if(fixCodingErrors)
	{
		uniqueMarkers <- unique(codingErrors$finals[,"Column"])
		finals[, uniqueMarkers] <- NA
		warning(paste0("Removing data for ", length(uniqueMarkers), " markers, because fixCodingErrors = TRUE was specified. For less aggressive removal, use listCodingErrors"))
	}
	geneticData <- new("geneticData", founders = founders, hetData = hetData, finals = finals, pedigree = pedigree)
	mpcross <- new("mpcross", geneticData = new("geneticDataList", list(geneticData)))
	return(mpcross)
}
#' @export
mpcrossMapped <- function(cross, map, rf=NULL)
{
	if(inherits(cross, "mpcrossRF"))
	{
		if(!is.null(rf))
		{
			stop("Two objects of class rf were specified")
		}
		return(new("mpcrossMapped", as(cross, "mpcross"), rf = cross@rf, map = map))
	}
	else
	{
		return(new("mpcrossMapped", cross, map = map, rf = rf))
	}
}
