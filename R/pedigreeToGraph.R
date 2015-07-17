pedigreeToGraph <- function(pedigree)
{
	founders <- nFounders(pedigree)
	nPedigreeLines <- length(pedigree@lineNames)
	motherEdges <- sapply((founders+1):nPedigreeLines, function(x) c(pedigree@lineNames[pedigree@mother[x]], pedigree@lineNames[x]))
	#We don't want to repeat edges if there's selfing
	fatherEdges <- sapply((founders+1):nPedigreeLines, 
		function(x)
		{
			if(pedigree@mother[x] == pedigree@father[x]) return(c())
			else return(c(pedigree@lineNames[pedigree@father[x]], pedigree@lineNames[x]))
		})
	#uneven lengths makes fatherEdges a list
	fatherEdges <- unlist(fatherEdges)
	graph <- graph.empty() + vertices(pedigree@lineNames) + edges(c(motherEdges, fatherEdges))
	laidOut <- layout.reingold.tilford(graph, root = pedigree@lineNames[1:founders])
	return(new("pedigreeGraph", graph = graph, layout = laidOut))
}