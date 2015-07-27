checkPedigreeGraph <- function(object)
{
	errors <- c()
	if(class(object@graph) != "igraph")
	{
		errors <- c(errors, "Slot graph must have class igraph")
	}
	if(!is.numeric(object@layout) || ncol(object@layout) != 2)
	{
		errors <- c(errors, "Slot layout must be a numeric matrix with two columns")
	}
	if(length(errors) > 0) return(errors)
	if(length(V(object@graph)) != nrow(object@layout))
	{
		errors <- c(errors, "Slot layout had the wrong number of rows")
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
.pedigreeGraph <- setClass("pedigreeGraph", slots = list(graph = "ANY", layout = "matrix"), validity = checkPedigreeGraph)