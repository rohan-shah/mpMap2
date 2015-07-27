#' @include pedigreeGraph-class.R
setMethod(f = "plot", signature = "pedigreeGraph", definition = function(x, y, ...)
{
	plot(x@graph, layout = x@layout, vertex.frame.color=NA, edge.arrow.mode="-", frame = 0, margin=0, ...)
})