#' @export
plotProbabilities <- function(inputObject, positions)
{
	if(inherits(inputObject, "mpcrossMapped"))
	{
		if(length(inputObject@geneticData) != 1)
		{
			stop("Input object contained more than one data set")
		}
		geneticData <- inputObject@geneticData[[1]]
	}
	else if(inherits(inputObject, "geneticData"))
	{
		geneticData <- inputObject
	}
	else
	{
		stop("Input object must have class mpcrossMapped or geneticData")
	}
	if(is.null(geneticData@probabilities))
	{
		stop("No probabilities were found")
	}
	if(!missing(positions))
	{
		geneticData@probabilities <- subset(geneticData@probabilities, positions = positions)
	}
	positionNames <- unlist(lapply(geneticData@probabilities@map, names))
	positions <- unlist(geneticData@probabilities@map)
	nPositions <- length(positionNames)
	
	transposed <- transposeProbabilities(geneticData)
	nLines <- nrow(transposed)
	nFounders <- nFounders(geneticData)
	offsetVector <- c(0, head(unlist(lapply(geneticData@probabilities@map, max)), -1))
	cumsumOffsetVector <- cumsum(offsetVector)
	offsetVector <- rep(cumsumOffsetVector, times = unlist(lapply(geneticData@probabilities@map, length)))
	if(geneticData@pedigree@selfing == "infinite")
	{
		dataSets <- list()
		for(founder in 1:nFounders)
		{
			averages <- apply(transposed[,(0:(nPositions-1))*nFounders + founder], 2, mean)
			dataSets[[founder]] <- data.frame(value = averages, positionName = positionNames, position = positions + offsetVector, chromosome = rep(names(geneticData@probabilities@map), times = unlist(lapply(geneticData@probabilities@map, length))), founder = rownames(founders(geneticData))[founder], stringsAsFactors = FALSE)
		}
		combined <- do.call(rbind, dataSets)
		combined$founder <- factor(combined$founder)
		lineData <- data.frame(position = tail(cumsumOffsetVector, -1))
		ggplot(combined, aes(position, value, colour = founder)) + geom_line() + xlab("Distance (cM)") + ylab("Probability") + geom_vline(data = lineData, mapping = aes(xintercept = position)) + facet_wrap(~ founder, ncol = 1) + theme(legend.position="none")
	}
	else if(geneticData@pedigree@selfing == "finite")
	{
		dataSets <- list()
		nAlleles <- length(unique(geneticData@probabilities@key[,3]))
		for(allele in 1:nAlleles)
		{
			averages <- apply(transposed[,(0:(nPositions-1))*nAlleles + allele], 2, mean)
			dataSets[[allele]] <- data.frame(value = averages, positionName = positionNames, position = positions + offsetVector, chromosome = rep(names(geneticData@probabilities@map), times = unlist(lapply(geneticData@probabilities@map, length))), allele = allele, stringsAsFactors = FALSE)
		}
		combined <- do.call(rbind, dataSets)
		combined$allele <- factor(combined$allele)
		lineData <- data.frame(position = tail(cumsumOffsetVector, -1))
		ggplot(combined, aes(position, value, colour = allele)) + geom_line() + xlab("Distance (cM)") + ylab("Probability") + geom_vline(data = lineData, mapping = aes(xintercept = position)) + facet_wrap(~ allele, ncol = 1) + theme(legend.position="none")
	}
	else
	{
		stop("Internal error")
	}
}
