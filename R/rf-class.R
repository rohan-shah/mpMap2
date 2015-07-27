setClassUnion("matrixOrNULL", c("matrix", "NULL"))
checkRF <- function(object)
{
	errors <- c()
	if(storage.mode(object@theta) != "double")
	{
		errors <- c(errors, "storage.mode of slot theta must be double")
	}
	if(!is.null(object@lod) && storage.mode(object@lod) != "double")
	{
		errors <- c(errors, "storage.mode of slot lod must be double")
	}
	if(!is.null(object@lkhd) && storage.mode(object@lkhd) != "double")
	{
		errors <- c(errors, "storage.mode of slot lkhd must be double")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.rf <- setClass("rf", slots = list(r = "numeric", theta = "matrix", lod = "matrixOrNULL", lkhd = "matrixOrNULL"), validity = checkRF)