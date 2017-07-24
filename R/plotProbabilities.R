#' @export
plotProbabilities <- function(inputObject, positions, alleles, chromosomes)
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
	if(!missing(positions) && !missing(chromosomes))
	{
		stop("Cannot use inputs positions and chromosomes together")
	}
	if(is.null(geneticData@probabilities))
	{
		stop("No probabilities were found")
	}
	if(!missing(chromosomes))
	{
		subsettedMap <- inputObject@map[names(inputObject@map) %in% chromosomes]
		if(length(subsettedMap) == 0 || any(!(chromosomes %in% names(inputObject@map)))) stop("Please enter valid chromosome names")
		positions <- unlist(lapply(subsettedMap, names))
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
		if(missing(alleles)) alleles <- 1:8
		dataSets <- list()
		for(founder in 1:nFounders)
		{
			averages <- apply(transposed[,(0:(nPositions-1))*nFounders + founder,drop=FALSE], 2, mean)
			dataSets[[founder]] <- data.frame(value = averages, positionName = positionNames, position = positions + offsetVector, chromosome = rep(names(geneticData@probabilities@map), times = unlist(lapply(geneticData@probabilities@map, length))), founder = rownames(founders(geneticData))[founder], stringsAsFactors = FALSE)
		}
		combined <- do.call(rbind, dataSets[alleles])
		combined$founder <- factor(combined$founder)
		lineData <- data.frame(position = tail(cumsumOffsetVector, -1))
		ggplot(combined, aes(position, value, colour = founder)) + geom_line() + xlab("Distance (cM)") + ylab("Probability") + geom_vline(data = lineData, mapping = aes(xintercept = position)) + facet_wrap(~ founder, ncol = 1) + theme(legend.position="none")
	}
	else if(geneticData@pedigree@selfing == "finite")
	{
		nAlleles <- length(unique(geneticData@probabilities@key[,3]))
		if(missing(alleles)) alleles <- 1:nAlleles
		dataSets <- list()
		for(allele in 1:nAlleles)
		{
			averages <- apply(transposed[,(0:(nPositions-1))*nAlleles + allele,drop=FALSE], 2, mean)
			if(allele <= nFounders) alleleName <- rownames(founders(geneticData))[allele]
			else
			{
				legendRow <- match(allele, geneticData@probabilities@key[,3])
				alleleName <- paste0(rownames(founders(geneticData))[geneticData@probabilities@key[legendRow, 1]], " - ", rownames(founders(geneticData))[geneticData@probabilities@key[legendRow, 2]])
			}
			dataSets[[allele]] <- data.frame(value = averages, positionName = positionNames, position = positions + offsetVector, chromosome = rep(names(geneticData@probabilities@map), times = unlist(lapply(geneticData@probabilities@map, length))), allele = alleleName, stringsAsFactors = FALSE)
		}
		combined <- do.call(rbind, dataSets[alleles])
		combined$allele <- factor(combined$allele)
		lineData <- data.frame(position = tail(cumsumOffsetVector, -1))
		ggplot(combined, aes(position, value, colour = allele)) + geom_line() + xlab("Distance (cM)") + ylab("Probability") + geom_vline(data = lineData, mapping = aes(xintercept = position)) + facet_wrap(~ allele, ncol = 1) + theme(legend.position="none")
	}
	else
	{
		stop("Internal error")
	}
}
