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
