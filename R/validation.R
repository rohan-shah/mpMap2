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