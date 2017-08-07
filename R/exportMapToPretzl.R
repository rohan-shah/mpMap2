#' @export
exportMapToPretzl <- function(inputObject, name, separateChromosomes = FALSE)
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
	if(separateChromosomes)
	{
		#Entry markers in the final json object
		toJSON <- lapply(1:length(map), function(x)
		{
			markerData <- list(data.frame(name = names(map[[x]]), position = as.vector(map[[x]]), stringsAsFactors = FALSE))
			names(markerData) <- names(map)[x]
			#Entry chromosomes in the final json object
			chromosomes <- data.frame(name = names(map)[x], stringsAsFactors = FALSE)
			chromosomes$markers <- markerData
			#Entry geneticMap in the final json object
			output <- list(chromosomes = chromosomes, name = paste0(name, "_", names(map)[x]))
			#Final json object
			output <- list(geneticmap = output)
			return(jsonlite::toJSON(output))
		})
		names(toJSON) <- names(map)
		return(toJSON)
	}
	else
	{
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
}
