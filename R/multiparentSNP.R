#' @include geneticData-class.R
#' @include mpcross-class.R
setClass("multiparentSNP", slots = list(keepHets = "logical"))
#' @export
multiparentSNP <- function(keepHets)
{
	return(new("multiparentSNP", keepHets = keepHets))
}
setMethod(f = "+", signature = c("geneticData", "multiparentSNP"), definition = function(e1, e2)
{
	nFounders <- nFounders(e1)
	if(nFounders == 2)
	{
		stop("multiparentSNP cannot be applied to a biparental design. Did you mean to use biparentDominant?")
	}
	keepHets <- e2@keepHets
	if(keepHets)
	{
		return(multiparentSNPKeepHets(e1))
	}
	else return(multiparentSNPRemoveHets(e1))
})
multiparentSNPRemoveHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(apply(e1@founders, 2, function(x) length(unique(x))) != nFounders))
	{
		stop("Can only apply multiparentSNP to a fully informative multiparent design (one that has a different genotypes for every homozygote, at every marker)")
	}

	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		numberOfOnes <- sample(nFounders-1, 1)
		oneAlleles <- sample(1:nFounders, numberOfOnes)
		oneAllelesOldValues <- copied@founders[oneAlleles, x]
		zeroAllelesOldValues <- copied@founders[-oneAlleles, x]

		becomesOne <- copied@finals[,x] %in% oneAllelesOldValues
		isHet <- !(copied@finals[,x] %in% copied@founders[,x])

		copied@finals[becomesOne, x] <<- 1
		copied@finals[!becomesOne, x] <<- 0
		copied@finals[isHet, x] <<- NA
		copied@founders[,x] <<- 0
		copied@founders[oneAlleles,x] <<- 1
		copied@hetData[[x]] <<- rbind(c(0,0,0), c(1,1,1))
	})
	return(copied)

}
multiparentSNPKeepHets <- function(e1)
{
	nFounders <- nFounders(e1)
	if(any(unlist(lapply(e1@hetData, function(x) length(unique(x[,3])))) != nFounders + nFounders*(nFounders-1)/2))
	{
		stop("Can only apply multiparentSNP to a fully informative multiparent design (one that has nFounders*nFounders possible genotypes at every marker)")
	}

	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		numberOfOnes <- sample(nFounders-1, 1)
		oneAlleles <- sample(1:nFounders, numberOfOnes)
		oneAllelesOldValues <- copied@founders[oneAlleles, x]
		zeroAllelesOldValues <- copied@founders[-oneAlleles, x]

		becomesOne <- copied@finals[,x] %in% oneAllelesOldValues
		becomesHetValues <- xor(copied@hetData[[x]][,1] %in% oneAllelesOldValues, copied@hetData[[x]][,2] %in% oneAllelesOldValues)

		copied@finals[becomesOne, x] <<- 1
		copied@finals[!becomesOne, x] <<- 0
		copied@finals[becomesHetValues, x] <<- 2
		copied@founders[,x] <<- 0
		copied@founders[oneAlleles,x] <<- 1
		copied@hetData[[x]] <<- rbind(c(0,0,0), c(1,1,1), c(0, 1, 2), c(1, 0, 2))
	})
	return(copied)
}
setMethod(f = "+", signature = c("mpcross", "multiparentSNP"), definition = function(e1, e2)
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
