#' @export
listCodingErrors <- function(founders, finals, hetData)
{
	errors <- .Call("listCodingErrors", founders, finals, hetData, PACKAGE="mpMap2")
	errors$finals <- errors$finals + 1
	errors$hetData <- errors$hetData + 1
	errors$null <- errors$null + 1
	return(errors)
}
