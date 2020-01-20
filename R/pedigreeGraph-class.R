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
#' @title Graph for a pedigree
#' @description Graph for a pedigree
#' @details This class contains the directed graph corresponding to a pedigree, and data for laying out the graph on a plane.
#' @slot graph An object of class igraph. 
#' @slot layout A matrix where each row gives the position of a graph vertex in the plane. 
.pedigreeGraph <- setClass("pedigreeGraph", slots = list(graph = "ANY", layout = "matrix"), validity = checkPedigreeGraph)
