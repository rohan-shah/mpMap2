#' @export
exportMapToPretzl <- function(inputObject, name)
{
	if(inherits(inputObject, "mpcrossMapped"))
	{
		map <- inputObject@map
	}
	else if(class(inputObject) == "map")
	{
		map <- inputObject
	}
	else
	{
		stop("inputObject must have class mpcrossMapped or map")
	}
	#Entry markers in the final json object
	output <- lapply(map, function(x)
	{
		data.frame(name = names(x), position = as.vector(x), stringsAsFactors = FALSE)
	})
	#Entry chromosomes in the final json object
	chromosomes <- data.frame(name = names(map), stringsAsFactors = FALSE)
	chromosomes$markers <- output
	#Entry geneticMap in the final json object
	output <- list(chromosomes = chromosomes, name = name)
	#Final json object
	output <- list(geneticmap = output)
	toJson <- jsonlite::toJSON(output)
	return(toJson)
}
