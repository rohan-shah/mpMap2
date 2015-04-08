nonNegativeIntegerArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || !((is.integer(x) || (abs(x - round(x)) < 1e-6)) && x>=0))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a non-negative integer"))
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
isMapArgument <- function(map)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(map) || !is(map, "map"))
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be a map object"))
	}
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
isMpMapMpcrossArgument <- function(x)
{
	call <- sys.call(sys.parent(0))
	parentCall <- sys.call(sys.parent(1))
	if(missing(x) || class(x) != "mpcross")
	{
		stop(paste0("Argument ", deparse(call[[2]]), " of ", deparse(parentCall[[1]]), " must be an mpcross object"))
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