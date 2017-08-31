#' @include geneticData-class.R
#' @include mpcross-class.R
setClass("multiparentSNPPrototype", slots = list(keepHets = "logical"))
multiparentSNPPrototype <- function(keepHets)
{
	return(new("multiparentSNPPrototype", keepHets = keepHets))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("geneticData", "multiparentSNPPrototype"), definition = function(e1, e2)
{
	nFounders <- nFounders(e1)
	if(nFounders == 2)
	{
		stop("multiparentSNPPrototype cannot be applied to a biparental design. Did you mean to use biparentDominant?")
	}
	keepHets <- e2@keepHets
	if(keepHets)
	{
		return(multiparentSNPPrototypeKeepHets(e1))
	}
	else return(multiparentSNPPrototypeRemoveHets(e1))
})
multiparentSNPPrototypeRemoveHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(apply(e1@founders, 2, function(x) length(unique(x))) != nFounders))
	{
		stop("Can only apply multiparentSNPPrototype to a fully informative multiparent design (one that has a different genotypes for every homozygote, at every marker)")
	}

	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		numberOfOnes <- sample(nFounders-1, 1)
		oneAlleles <- sample(1:nFounders, numberOfOnes)
		oneAllelesOldValues <- copied@founders[oneAlleles, x] 
		zeroAllelesOldValues <- copied@founders[-oneAlleles, x]

		currentHetData <- e1@hetData[[x]]
		oneAllelesOldValues <- currentHetData[currentHetData[,1] %in% oneAllelesOldValues & currentHetData[,2] %in% oneAllelesOldValues,3]
		zeroAllelesOldValues <- currentHetData[currentHetData[,1] %in% zeroAllelesOldValues & currentHetData[,2] %in% zeroAllelesOldValues,3]

		becomesOne <- copied@finals[,x] %in% oneAllelesOldValues
		becomesZero <- copied@finals[,x] %in% zeroAllelesOldValues
		hetValues <- currentHetData[xor(currentHetData[,1] %in% oneAllelesOldValues, currentHetData[,2] %in% oneAllelesOldValues), 3]
		isHet <- copied@finals[,x] %in% hetValues

		copied@finals[becomesOne, x] <<- 1L
		copied@finals[becomesZero, x] <<- 0L
		copied@finals[isHet, x] <<- NA
		copied@founders[,x] <<- 0L
		copied@founders[oneAlleles,x] <<- 1L
		copied@hetData[[x]] <<- rbind(c(0L,0L,0L), c(1L,1L,1L))
	})
	return(copied)

}
multiparentSNPPrototypeKeepHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(unlist(lapply(e1@hetData, function(x) length(unique(x[,3])))) != nFounders + nFounders*(nFounders-1)/2))
	{
		stop("Can only apply multiparentSNPPrototype to a fully informative multiparent design (one that has nFounders*nFounders possible genotypes at every marker)")
	}

	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		numberOfOnes <- sample(nFounders-1, 1)
		oneAlleles <- sample(1:nFounders, numberOfOnes)
		oneAllelesFounderValues <- copied@founders[oneAlleles, x]
		zeroAllelesFounderValues <- copied@founders[-oneAlleles, x]

		oneAllelesFinalValues <- copied@hetData[[x]][(copied@hetData[[x]][,1] %in% oneAllelesFounderValues) & (copied@hetData[[x]][,2] %in% oneAllelesFounderValues),3]
		zeroAllelesFinalValues <- copied@hetData[[x]][(copied@hetData[[x]][,1] %in% zeroAllelesFounderValues) & (copied@hetData[[x]][,2] %in% zeroAllelesFounderValues),3]

		becomesOne <- copied@finals[,x] %in% oneAllelesFinalValues
		becomesZero <- copied@finals[,x] %in% zeroAllelesFinalValues
		hetValues <- copied@hetData[[x]][xor(copied@hetData[[x]][,1] %in% oneAllelesFounderValues, copied@hetData[[x]][,2] %in% oneAllelesFounderValues),3]
		becomesHetValues <- copied@finals[,x] %in% hetValues

		copied@finals[becomesOne, x] <<- 1L
		copied@finals[becomesZero, x] <<- 0L
		copied@finals[becomesHetValues, x] <<- 2L
		copied@founders[,x] <<- 0L
		copied@founders[oneAlleles,x] <<- 1L
		copied@hetData[[x]] <<- rbind(c(0L,0L,0L), c(1L,1L,1L), c(0L, 1L, 2L), c(1L, 0L, 2L))
	})
	return(copied)
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "multiparentSNPPrototype"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Assigning SNP marker patterns will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	if(length(e1@geneticData) > 1)
	{
		stop("Attempting to change an object containing multiple data sets. Please change each dataset individually")
	}
	if(nFounders(e1) == 2)
	{
		stop("multiparentSNP cannot be applied to a biparental design. Did you mean to use biparentDominant?")
	}
	e1@geneticData[[1]] <- e1@geneticData[[1]]+e2
	return(e1)
})
