#' @include geneticData-class.R
#' @include mpcross-class.R
setClass("biparentalDominant", contains="NULL")
#' Make markers in a biparental cross dominant
#'
#' Change the markers in a biparental cross from fully informative to dominant. The dominant founder is chosen randomly for every marker. The transformation is applied to an object using the addition operator, see the example below for details.
#' @return An object of internal type \code{biparentalDominant}, which can be combined with an object of class \code{mpcross} using the addition operator. 
#' @examples
#' #Simulate an F2 design
#' f2Pedigree <- f2Pedigree(1000)
#' map <- qtl::sim.map(len = 100, n.mar = 11, include.x=FALSE)
#' cross <- simulateMPCross(map = map, pedigree = f2Pedigree, mapFunction = haldane, seed = 1)
#' founders(cross)
#' finals(cross)[1:10,]
#' #The heterozygotes are initially coded as 3
#' hetData(cross)[[1]]
#' #Make all markers dominant
#' dominantCross <- cross + biparentalDominant()
#' founders(dominantCross)
#' finals(dominantCross)[1:10,]
#' #The heterozygotes are now coded the same as one of the homozygotes
#' hetData(dominantCross)[1:4]
#' @export
biparentalDominant <- function()
{
	return(new("biparentalDominant"))
}
#' @rdname internalOperators
setMethod(f = "+", signature = c("geneticData", "biparentalDominant"), definition = function(e1, e2)
{
	#We only want to be able to apply this to a dataset that's fully informative. This is equivalent to the hetData having three unique values in the third column.
	if(any(unlist(lapply(e1@hetData, function(x) length(unique(x[,3])))) != 3))
	{
		stop("Can only apply biparentalDominant to a fully informative biparental design (one that has three possible genotypes at every marker)")
	}
	copied <- e1
	nMarkers <- nMarkers(copied)
	sapply(1:nMarkers, function(x)
	{
		hetDataThisMarker <- copied@hetData[[x]]
		isHetRow <- hetDataThisMarker[,1] != hetDataThisMarker[,2]
		hetCode <- unique(hetDataThisMarker[isHetRow, 3])
		homozygoteCodes <- hetDataThisMarker[!isHetRow, 3]
		dominant <- sample(homozygoteCodes, 1)

		copied@hetData[[x]][isHetRow,3] <<- dominant
		copied@finals[copied@finals[,x] == hetCode,x] <<- dominant
	})
	return(copied)
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "biparentalDominant"), definition = function(e1, e2)
{
	if(class(e1) != "mpcross")
	{
		warning("Assigning random dominant marker patterns will remove all data except genetic data")
	}
	e1 <- as(e1, "mpcross")
	if(length(e1@geneticData) > 1)
	{
		stop("Attempting to change an object containing multiple data sets. Please change each dataset individually")
	}
	if(nFounders(e1) != 2)
	{
		stop("biparentalDomimant can only be applied to a biparental design")
	}
	e1@geneticData[[1]] <- e1@geneticData[[1]]+e2
	return(e1)
})
