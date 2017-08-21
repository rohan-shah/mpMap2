checkMap <- function(object)
{
	if(length(object) == 0)
	{
		return("Map object must contain at least one chromosome")
	}
	if(any(unlist(lapply(object, length)) == 0))
	{
		return("All chromosomes of a map must contain at least one marker")
	}
	allNumeric <- unlist(lapply(object, is.numeric))
	if(!all(allNumeric)) return("A map object must be a list of numeric vectors")
	
	if(length(unique(names(object))) != length(object))
	{
		return("Chromosome names must be unique")
	}

	if(length(unique(unlist(lapply(object, names)))) != length(unlist(object)))
	{
		return("Marker names must be unique")
	}
	if(any(unlist(lapply(object, function(x) order(x) != 1:length(x)))))
	{
		return("Marker positions must be in ascending order within every chromosome")
	}
	return(TRUE)
}
.map4 <- setClass("map4", contains = "list", validity = checkMap)
setOldClass("map", S4Class = "map4")
removeClass("map4")
