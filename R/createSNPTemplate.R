createSNPTemplate <- function(inputObject, newData, hetEncoding, markerName)
{
	if(inherits(inputObject, "geneticData"))
	{
		geneticData <- inputObject
	}
	else if(inherits(inputObject, "mpcross"))
	{
		if(length(inputObject@geneticData[[1]]) > 1)
		{
			stop("If inputObject has class mpcross, it must contain only a single dataset")
		}
		geneticData <- inputObject@geneticData[[1]]
	}
	else
	{
		stop("Input object must have class mpcross or geneticData")
	}
	if(is.matrix(newData) || is.data.frame(newData))
	{
		if(ncol(newData) != 1) stop("Input newData must have a single column")
		newData <- newData[,1]
	}
	if(is.null(names(newData)))
	{
		stop("Input data must have names")
	}
	if(!is.integer(newData))
	{
		stop("Input newData must contain integer values")
	}
	founders <- cbind(newData[rownames(founders(geneticData))])
	colnames(founders) <- markerName
	alleles <- na.omit(unique(founders))
	if(missing(hetEncoding) || is.null(hetEncoding))
	{
		hetData <- rbind(cbind(alleles, alleles, alleles))
		warning("Input hetEncoding was missing, so no heterozygote encodings are generated. Was this desired?")
	}
	else if(length(hetEncoding) == 1)
	{
		hetData <- rbind(cbind(alleles, alleles, alleles), c(alleles[1], alleles[2], hetEncoding), c(alleles[2], alleles[1], hetEncoding))
	} 
	else
	{
		if(is.null(dim(hetEncoding)) || ncol(hetEncoding) != 3)
		{
			stop("Input hetEncoding must be a single number, or a matrix with three columns")
		}
		hetData <- rbind(cbind(alleles, alleles, alleles), hetEncoding)
	}
	dimnames(hetData) <- NULL
	hetData <- list(hetData)
	names(hetData) <- markerName
	hetData <- new("hetData", hetData)

	finals <- cbind(newData[lineNames(geneticData)])
	rownames(finals) <- lineNames(geneticData)
	colnames(finals) <- markerName
	
	founders <- cbind(newData[rownames(founders(geneticData))])
	colnames(founders) <- markerName
	result <- mpcross(finals = finals, founders = founders, hetData = hetData, pedigree = geneticData@pedigree)
	return(result)
}
