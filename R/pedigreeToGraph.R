#' @title Convert pedigree to a graph
#' @description Convert pedigree to a graph
#' @details It is often useful for visualisation purposes to generate the pedigree graph. In this graph, every genetic line is a vertex in a graph, and there is an edge from every parent to all the offspring. This function generates the graph, and lays the graph out in the plane in a way that tends to make the structure of the graph as clear as possible. 
#' @param pedigree The pedigree to convert into a graph
#' @return An object of class \code{pedigreeGraph}, containing the graph and a planar layout for the graph. 
#' @export
pedigreeToGraph <- function(pedigree)
{
	founders <- nFounders(pedigree)
	nPedigreeLines <- length(pedigree@lineNames)

	motherEdges <- sapply((founders+1):nPedigreeLines, 
		function(x) 
		{
			if(pedigree@mother[x] == pedigree@father[x]) return(c())
			else c(pedigree@lineNames[pedigree@mother[x]], pedigree@lineNames[x])
		})
	fatherEdges <- sapply((founders+1):nPedigreeLines, 
		function(x)
		{
			if(pedigree@mother[x] == pedigree@father[x]) return(c())
			else return(c(pedigree@lineNames[pedigree@father[x]], pedigree@lineNames[x]))
		})
	selfingEdges <- sapply((founders+1):nPedigreeLines, 
		function(x)
		{
			if(pedigree@mother[x] == pedigree@father[x]) return(c(pedigree@lineNames[pedigree@father[x]], pedigree@lineNames[x]))
			else return(c())
		})
	
	#uneven lengths make these lists
	fatherEdges <- unlist(fatherEdges)
	motherEdges <- unlist(motherEdges)
	selfingEdges <- unlist(selfingEdges)

	graph <- graph.empty() + vertices(pedigree@lineNames) + edges(c(motherEdges, fatherEdges, selfingEdges))
	type <- c(rep("maternal", length(motherEdges)/2), rep("paternal", length(fatherEdges)/2), rep("selfing", length(selfingEdges)/2))
	E(graph)$type <- type

	laidOut <- layout.reingold.tilford(graph, root = pedigree@lineNames[1:founders])
	return(new("pedigreeGraph", graph = graph, layout = laidOut))
}
