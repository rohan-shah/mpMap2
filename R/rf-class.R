#' @include rawSymmetricMatrix.R
setClassUnion("dspMatrixOrNULL", c("dspMatrix", "NULL"))
checkRF <- function(object)
{
	errors <- c()
	thetaMarkers <- object@theta@markers
	if(!is.null(object@lod))
	{
		if(is.null(colnames(object@lod)))
		{
			errors <- c(errors, "Slot @lod must have column names")
		}
		else if(nrow(object@lod) != length(thetaMarkers))
		{
			errors <- c(errors, "Dimensions of @lod were inconsistent with those of @theta")
		}
		else if(any(thetaMarkers != colnames(object@lod)))
		{
			errors <- c(errors, "Column names of @lod were inconsistent with those of @theta")
		}
	}
	if(!is.null(object@lkhd))
	{
		if(is.null(colnames(object@lkhd)))
		{
			errors <- c(errors, "Slot @lkhd must have column names")
		}
		else if(nrow(object@lkhd) != length(thetaMarkers))
		{
			errors <- c(errors, "Dimensions of @lkhd were inconsistent with those of @theta")
		}
		else if(any(thetaMarkers != colnames(object@lkhd)))
		{
			errors <- c(errors, "Column names of @lkhd were inconsistent with those of @theta")
		}

	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.rf <- setClass("rf", slots = list(theta = "rawSymmetricMatrix", lod = "dspMatrixOrNULL", lkhd = "dspMatrixOrNULL", gbLimit = "numeric"), validity = checkRF)
