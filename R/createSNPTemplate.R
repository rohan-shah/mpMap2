#' @export
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
	founders <- cbind(newData[rownames(founders(geneticData))])
	colnames(founders) <- markerName
	alleles <- na.omit(unique(founders))
	if(length(alleles) > 2)
	{
		stop("Cannot call createSNPTemplate with more than two homozygous alleles")
	}
	hetData <- rbind(cbind(alleles, alleles, alleles), c(alleles[1], alleles[2], hetEncoding), c(alleles[2], alleles[1], hetEncoding))
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
