#' @export
listCodingErrors <- function(founders, finals, hetData)
{
	errors <- .Call("listCodingErrors", founders, finals, hetData, PACKAGE="mpMap2")
	errors$finals <- errors$finals + 1
	errors$hetData <- errors$hetData + 1
	errors$null <- errors$null + 1
	return(errors)
}
#' @export
listCodingErrorsMpMap <- function(founders, finals)
{
	newHetDataList <- lapply(as.list(1:ncol(founders)), function(x)
	{
		if(any(is.na(founders[,x])))
		{
			finals[,x] <<- NA
			return(matrix(0L, 0, 3))
		}
		else
		{
			uniqueAlleles <- unique(founders[, x])
			retVal <- cbind(uniqueAlleles, uniqueAlleles, uniqueAlleles)
			colnames(retVal) <- NULL
			return(retVal)
		}
	})
	names(newHetDataList) <- colnames(founders)

	newHetData <- new("hetData", newHetDataList)
	errors <- .Call("listCodingErrors", founders, finals, newHetData, PACKAGE="mpMap2")
	errors$finals <- errors$finals + 1
	errors$hetData <- errors$hetData + 1
	errors$null <- errors$null + 1
	return(errors)
}
