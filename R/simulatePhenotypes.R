setClass("normalPhenotype", contains="NULL", slots = list("means" = "numeric", "standardDeviations" = "numeric", phenotypeName = "character", marker = "character"))
#' @rdname internalOperators
setMethod(f = "+", signature = c("geneticData", "normalPhenotype"), definition = function(e1, e2)
{
	if(!(e2@marker %in% colnames(e1@finals)))
	{
		stop("Marker not found")
	}
	if(!is.null(e1@pheno) && e2@phenotypeName %in% colnames(e1@pheno))
	{
		stop(paste0("Phenotype ", e2@phenotypeName, " was already present"))
	}
	relevantHetData <- e1@hetData[[e2@marker]]
	relevantAlleles <- as.character(unique(relevantHetData[,3]))
	if(!isTRUE(all.equal(sort(relevantAlleles), sort(names(e2@means)))))
	{
		#Didn't give means/standard deviations for all alleles
		if(length(setdiff(sort(relevantAlleles), sort(names(e2@means)))) > 0)
		{
			stop(paste0("A value must be given for every possible allele of marker ", e2@marker, "; See hetData[[\"", e2@marker, "\"]][,3] for possible alleles"))
		}
		#Gave means/standard deviations for too many alleles
		else
		{
			stop(paste0("Means and standardDeviations were given for impossible alleles of marker ", e2@marker, "; See hetData[[\"", e2@marker, "\"]][,3] for possible alleles"))
		}
	}
	#Don't warn about NAs
	suppressWarnings(extraColumn <- rnorm(n = nrow(e1@finals), mean = e2@means[as.character(e1@finals[,e2@marker])], sd = e2@standardDeviations[as.character(e1@finals[,e2@marker])]))
	if(is.null(e1@pheno))
	{
		e1@pheno <- data.frame(extraColumn)
		rownames(e1@pheno) <- rownames(e1@finals)
		colnames(e1@pheno) <- e2@phenotypeName
	}
	else
	{
		argsList <- list(e1@pheno, extraColumn)
		names(argsList)[2] <- e2@phenotypeName
		e1@pheno <- do.call(cbind, argsList)
		rownames(e1@pheno) <- rownames(e1@finals)
	}
	return(e1)
})
#' @rdname internalOperators
setMethod(f = "+", signature = c("mpcross", "normalPhenotype"), definition = function(e1, e2)
{
	lapply(1:length(e1@geneticData), function(x) e1@geneticData[[x]] <<- e1@geneticData[[x]] + e2)
	return(e1)
})
#' @title Simulate normally distributed phenotype
#' @description Add a normally distributed phenotype 
#' @details Add a normally distributed phenotype to a given populations
#' @param means The means of the phenotype for all the different founder alleles
#' @param standardDeviations The standard deviations of the phenotype for all the different founder alleles
#' @param phenotypeName The name of the new phenotype
#' @param marker The name of the marker which controls this phenotype
#' @return An object of class \code{normalPhenotype} representing the phenotype. 
#' @export
normalPhenotype <- function(means, standardDeviations, phenotypeName, marker)
{
	if(missing(means) || missing(standardDeviations) || missing(phenotypeName) || missing(marker))
	{
		stop("Inputs means, standardDeviations, phenotypeName and marker are required")
	}
	if(length(means) != length(standardDeviations))
	{
		stop("The number of means and standardDeviations entered must be the same")
	}
	if(is.null(names(means))) names(means) <- as.character(1:length(means))
	if(is.null(names(standardDeviations))) names(standardDeviations) <- as.character(1:length(standardDeviations))
	if(!isTRUE(all.equal(sort(names(means)), sort(names(standardDeviations)))))
	{
		stop("The names of inputs means and standardDeviations must match") 
	}
	return(new("normalPhenotype", means = means, standardDeviations = standardDeviations, phenotypeName = phenotypeName, marker = marker))
}
