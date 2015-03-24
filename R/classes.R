check_pedigree <- function(object)
{
	nLines <- length(object@lineNames)
	errors <- c()
	if(length(object@mother) != nLines || length(object@father) != nLines)
	{
		errors <- c(errors, "Inconsistent lengths for slots lineNames, mother and father")
	}

	if(any(is.na(object@initial)))
	{
		errors <- c("errors", "Slot initial contained NA values")
	}
	if(any(is.na(object@lineNames)))
	{
		errors <- c("errors", "Slot lineNames contained NA values")
	}
	if(any(is.na(object@mother)))
	{
		errors <- c("errors", "Slot mother contained NA values")
	}
	if(any(is.na(object@father)))
	{
		errors <- c("errors", "Slot father contained NA values")
	}

	if(length(object@initial) == 0)
	{
		errors <- c(errors, "Slot initial must be non-empty")
	}

	if(any(object@initial < 0 | object@initial > nLines))
	{
		errors <- c(errors, "Values in slot initial had invalid values")
	}
	if(any(object@mother < 0 | object@mother > nLines))
	{
		errors <- c(errors, "Values in slot mother had invalid values")
	}
	if(any(object@father < 0 | object@father > nLines))
	{
		errors <- c(errors, "Values in slot father had invalid values")
	}

	if(any(object@mother[object@initial] != 0))
	{
		errors <- c(errors, "An entry in slot mother for a line named in slot initial cannot have a mother")
	}
	if(any(object@father[object@initial] != 0))
	{
		errors <- c(errors, "An entry in slot father for a line named in slot initial cannot have a father")
	}

	if(length(unique(object@initial)) != length(object@initial))
	{
		errors <- c(errors, "Slot initial cannot contain duplicate values")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.pedigree <- setClass("pedigree", representation(lineNames = "character", mother = "integer", father = "integer", initial = "integer"), validity = check_pedigree)
check_mpcross <- function(object)
{
	errors <- c()
	if(!is.integer(object@founders))
	{
		errors <- c(errors, "Slot founders must be an integer matrix")
	}
	if(!is.integer(object@finals))
	{
		errors <- c(errors, "Slot finals must be an integer matrix")
	}
}
.mpcross <- setClass("mpcross", representation(finals = "matrix", founders = "matrix", pedigree = "pedigree"))