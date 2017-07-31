nonNegativeIntegerArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !((is.integer(x) || (abs(x - round(x)) < 1e-6)) && x>=0))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a non-negative integer"))
	}
}
positiveIntegerArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !((is.integer(x) || (abs(x - round(x)) < 1e-6)) && x>0.5))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a positive integer"))
	}
}
isPedigreeArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !is(x, "pedigree"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a pedigree object"))
	}
}
isDetailedPedigreeArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !is(x, "detailedPedigree"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a detailedPedigree object"))
	}
}
isMapArgument <- function(map)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(map) || !is(map, "map"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a map object"))
	}
}
isNumericVectorListArgument <- function(argument)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(argument))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a list of numeric vectors"))
	}
	lapply(argument, function(x)
	{
		if(storage.mode(x) != "double")
		{
			stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a list of numeric vectors"))
		}
	})
}
isIntegerMatrixArgument <- function(matrix)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(matrix) || typeof(matrix) != "integer" || length(dim(matrix)) != 2)
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an integer matrix"))
	}	
}
isNumericMatrixArgument <- function(matrix)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(matrix) || !is.numeric(matrix) || length(dim(matrix)) != 2)
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a numeric matrix"))
	}	
}
isOldMpMapMpcrossArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !inherits(x, "mpcross"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcross object"))
	}		
}
inheritsNewMpcrossArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !isS4(x) || !inherits(x, "mpcross"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcross object"))
	}
}
isNewMpcrossLGArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || (class(x) != "mpcrossLG" && !canCoerce(x, "mpcrossLG")))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcrossLG object"))
	}
}
isNewMpcrossMappedArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !inherits(x, "mpcrossMapped"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcrossMapped object"))
	}
}
isNewMpcrossRFArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || (class(x) != "mpcrossLG" && class(x) != "mpcrossRF" && class(x) != "mpcrossMapped"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcrossRF object"))
	}
	if(class(x) == "mpcrossMapped" && is.null(x@rf))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " cannot be automatically converted to an object of class mpcrossRF"))
	}
	if(class(x) == "mpcrossLG" && is.null(x@rf))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " cannot be automatically converted to an object of class mpcrossRF"))
	}
}
isMpMapPedigreeArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || class(x) != "pedigree")
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a pedigree object"))
	}		
}
isHetDataArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || class(x) != "hetData")
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a hetData object"))
	}
}
