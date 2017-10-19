setClass("assignFounderPatternPrototype", slots = list(data = "matrix"))
assignFounderPatternPrototype <- function(founderMatrix)
{
	return(new("assignFounderPatternPrototype", data = founderMatrix))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("geneticData", "assignFounderPatternPrototype"), definition = function(e1, e2)
{
	assignFounderPatternPrototypeInternal(geneticData = e1, founderPattern = e2@data)
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "assignFounderPatternPrototype"), definition = function(e1, e2)
{
	if(length(e1@geneticData) != 1)
	{
		stop("Can only apply assignFounderPatternPrototype to an mpcross object with a single experiment")
	}
	result <- assignFounderPatternPrototypeInternal(geneticData = e1@geneticData[[1]], founderPattern = e2@data)
	geneticData <- new("geneticData", founders = result@founders, finals = result@finals, pedigree = e1@geneticData[[1]]@pedigree, hetData = result@hetData)
	geneticDataList <- new("geneticDataList", list(geneticData))
	return(new("mpcross", geneticData = geneticDataList))
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcrossMapped", "assignFounderPatternPrototype"), definition = function(e1, e2)
{
	if(length(e1@geneticData) != 1)
	{
		stop("Can only apply assignFounderPatternPrototype to an mpcross object with a single experiment")
	}
	result <- assignFounderPatternPrototypeInternal(geneticData = e1@geneticData[[1]], founderPattern = e2@data)
	geneticData <- new("geneticData", founders = result@founders, finals = result@finals, pedigree = e1@geneticData[[1]]@pedigree, map = e1@map, hetData = result@hetData)
	geneticDataList <- new("geneticDataList", list(geneticData))
	return(new("mpcrossMapped", geneticData = geneticDataList,  map = e1@map))
})
assignFounderPatternPrototypeInternal <- function(geneticData, founderPattern)
{
	nFounders <- nFounders(geneticData)
	if(!identical(dim(founderPattern), dim(geneticData@founders)))
	{
		stop("Input founderPattern had the wrong dimensions")
	}
	expectedFounders <- rep(1:nFounders, times = ncol(founderPattern))
	dim(expectedFounders) <- c(nFounders, ncol(founderPattern))
	if(any(expectedFounders != geneticData@founders))
	{
		stop("Can only apply assignFounderPattern to an object with fully informative founders")
	}
	#Create new object
	newObject <- geneticData
	newObject@founders <- founderPattern
	newObject@finals <- geneticData@finals

	newHetData <- vector(mode = "list", length = nMarkers(geneticData))
	names(newHetData) <- markers(geneticData)

	#Convert the hetData and finals data across. 
	sapply(1:ncol(founderPattern), function(x)
	{
		currentHetData <- geneticData@hetData[[x]]
		#Create a translation vector for the finals data. The hetData will be done separately. 
		translationVector <- vector(mode = "integer", length = nFounders*(nFounders+1)/2)
		#initially a value of -1 for hets. 
		translationVector[] <- -1L
		translationVector[1:nFounders] <- founderPattern[,x]

		allAlleles <- unique(currentHetData[,3])
		if(any(allAlleles > nFounders*(nFounders+1)/2) || any(allAlleles < 1))
		{
			stop("An encoding for an allele was out of range. All alleles should be between 0 and (nFounders*(nFounders+1)/2 + 1)")
		}
		nNewHomozygotes <- length(unique(founderPattern[,x]))
		newHet <- as.integer(max(unique(founderPattern[,x])) + 1)

		currentNewHetData <- matrix(0L, nrow = (length(allAlleles) - nFounders)*2 + nFounders, ncol = 3)
		currentNewHetData[,1] <- translationVector[currentHetData[,1]]
		currentNewHetData[,2] <- translationVector[currentHetData[,2]]

		#Every time we encounter a new heterozygote, put it into the translation vector. 
		for(y in 1:nrow(currentHetData))
		{
			if(translationVector[currentHetData[y, 3]] == -1)
			{
				firstRow <- head(which((currentNewHetData[,1] == currentNewHetData[y,1]) & (currentNewHetData[,2] == currentNewHetData[y,2])), 1)
				if(translationVector[currentHetData[firstRow, 3]] == -1)
				{
					translationVector[currentHetData[y, 3]] <- newHet
					newHet <- newHet + 1L
				}
				else
				{
					translationVector[currentHetData[y, 3]] <- translationVector[currentHetData[firstRow, 3]]
				}
			}
		}
		currentNewHetData[,3] <- translationVector[currentHetData[,3]]
		currentNewHetData <- currentNewHetData[!duplicated(currentNewHetData),]
		newObject@finals[,x] <<- translationVector[newObject@finals[,x]]
		newHetData[[x]] <<- currentNewHetData
	})
	newHetData <- new("hetData", newHetData)
	newObject@hetData <- newHetData
	return(newObject)
}
