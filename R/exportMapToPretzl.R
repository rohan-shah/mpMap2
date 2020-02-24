#' @title Export genetic map to Pretzl
#' @description Export genetic map to Pretzl
#' @details Convert the genetic map from an object of class \code{mpcrossMapped} to the JSON format used by Pretzl. Pretzl is a web app for visualising and comparing genetic maps. 
#' @param inputObject The object of class \code{mpcrossMapped} containing the genetic map
#' @param name If a single JSON object is being exported, the name of the exported map. 
#' @param separateChromosomes If \code{TRUE}, separate exports will be generated for each chromosome. The name associated with each chromosome map will contain the chromosome name as a suffix. 
#' @return A list containing JSON, suitable for import into Pretzl. 
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
