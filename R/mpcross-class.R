#' @include rf-class.R
#' @include lg-class.R
#' @include map-class.R
#' @include geneticData-class.R
#' @include canSkipValidity.R
checkCompatibleGeneticData <- function(objects)
{
	expectedMarkers <- markers(objects[[1]])
	errors <- c()
	for(i in 1:length(objects))
	{
		x <- objects[[i]]
		if(nMarkers(x) != length(expectedMarkers))
		{
			errors <- c(errors, paste0("Wrong number of markers in slot geneticData[[", i, "]]"))
		}
		else if(any(markers(x) != expectedMarkers))
		{
			errors <- c(errors, paste0("Wrong markers in slot geneticData[[", i, "]]"))
		}
	}

	return(errors)
}
checkMpcross <- function(object)
{
	errors <- c()
	if(length(object@geneticData) == 0)
	{
		errors <- c(errors, "Must contain at least one set of genetic data")
	}
	else
	{
		for(i in 1:length(object@geneticData))
		{
			x <- object@geneticData[[i]]
			if(class(x) != "geneticData")
			{
				errors <- c(errors, paste0("Value of slot geneticData[[", i, "]] must be a geneticData object"))
			}
			else
			{
				geneticErrors <- checkGeneticData(x)
				if(mode(geneticErrors) != "logical")
				{
					errors <- c(errors, paste0("Error in geneticData[[", i, "]]: ", geneticErrors))
				}
			}
		}
		errors <- c(errors, checkCompatibleGeneticData(object@geneticData))
	}
	allLineNames <- unlist(lapply(object@geneticData, function(x) rownames(x@finals)))
	if(length(allLineNames) != length(unique(allLineNames)))
	{
		errors <- c(errors, "Line names must be unique")
	}
	#Check that the line names are unique
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
#' @title A collection of multi-parent populations without a genetic map
#' @description A collection of multi-parent populations without a genetic map
#' @details An object of class \code{mpcross} contains data about one or more multi-parent populations, without a genetic map. As there is no genetic map, there is no information about IBD imputed genotypes or IBD genotype probabilities. There is also no information about estimated recombination fractions. 
#' 
#' A \code{mpcross} object must contain (at a minimum) genetic data about the founding lines of the population, genetic lines about the final lines of the population, a pedigree with information about how the final lines were generated from the founding lines, and information about how heterozygotes have been encoded. See \code{\link{geneticData-class}} for further information. See \code{\link{mpcross}} for the constructor function. 
#' @slot geneticData A list of objects of class \code{\link{geneticData-class}}, each representing a population.
.mpcross <- setClass("mpcross", slots = list(geneticData = "geneticDataList"), validity=checkMpcross, contains = "canSkipValidity")

checkMpcrossRF <- function(object)
{
	errors <- c()
	if(length(object@rf@theta@markers) != length(markers(object)))
	{
		errors <- "Inconsistent number of markers in slots @geneticData[[1]] and @rf@theta"
	}
	else if(any(object@rf@theta@markers != markers(object)))
	{
		errors <- "Inconsistent markers in slots @geneticData[[1]] and @rf@theta"
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
setClassUnion("rfOrNULL", c("rf", "NULL"))
#' @title A collection of multi-parent populations with recombination fraction estimates
#' @description A collection of multi-parent populations with recombination fraction estimates
#' @details An object of class \code{mpcrossRF} contains data about one or more multi-parent populations, without a genetic map, but with recombination fraction estimates. As there is no genetic map, there is no information about IBD imputed genotypes or IBD genotype probabilities. 
#' @slot geneticData A list of objects of class \code{geneticData}, each representing a population.
#' @slot rf Estimates of recombination fractions between every pair of genetic markers. 
.mpcrossRF <- setClass("mpcrossRF", contains = "mpcross", slots = list(rf = "rf"), validity=checkMpcrossRF)

checkMpcrossLG <- function(object)
{
	if(is.null(names(object@lg@groups)))
	{
		return("Slot lg@groups must have names")
	}
	if(any(is.na(names(object@lg@groups)) | names(object@lg@groups) != markers(object)))
	{
		return("Marker names implied by names of slots lg@groups and founders were different")
	}
	if(!is.null(object@lg@imputedTheta) && !is.null(object@rf))
	{
		if(length(object@lg@imputedTheta) > 0)
		{
			if(!identical(object@lg@imputedTheta[[1]]@levels, object@rf@theta@levels))
			{
				return("Slots lg@imputedTheta and rf@theta must have the same levels")
			}
		}
	}
	return(TRUE)
}
.mpcrossLG <- setClass("mpcrossLG", contains = "mpcross", slots = list(lg = "lg", rf = "rfOrNULL"), validity=checkMpcrossLG)

checkMpcrossMapped <- function(object)
{
	errors <- c()
	nMapMarkers <- sum(unlist(lapply(object@map, length)))
	if(nMarkers(object) != nMapMarkers)
	{
		errors <- c(errors, "Number of markers in map is different from the number of markers in slot founders")
	}
	else
	{
		markerNames <- unlist(lapply(object@map, names))
		names(markerNames) <- NULL
		for(i in 1:length(object@geneticData))
		{
			geneticDataMarkers <- markers(object@geneticData[[i]])
			names(geneticDataMarkers) <- NULL
			if(!identical(geneticDataMarkers, markerNames))
			{
				errors <- c(errors, "Marker names in the genetic data must be the same as those in the map, and must occur in the same order as unlist(lapply(map, names))")
			}
		}
	}
	if(!is.null(object@rf))
	{
		if(!identical(markers(object), markers(object@rf)))
		{
			errors <- c(errors, "Slot @rf had markers that were inconsistent with the genetic data")
		}
	}
	if(length(errors) > 0) return(errors)
	return(TRUE)
}
#' @title A collection of multi-parent populations with a genetic map
#' @description A collection of multi-parent populations with a genetic map
#' @details An object of class \code{mpcrossMapped} contains genetic data for one or more populations, and a genetic map. It might also contain data about the recombination fractions between markers (or it might not). 
#' @slot map The genetic map for all populations
#' @slot rf The recombination fraction data (which might be NULL). 
#' @slot geneticData A list of objects of class \code{geneticData}, each representing a population. 
.mpcrossMapped <- setClass("mpcrossMapped", contains = "mpcross", slots = list(map = "map", rf = "rfOrNULL"), validity=checkMpcrossMapped)

setAs("mpcrossMapped", "mpcrossRF", def = function(from, to)
	{
		if(is.null(from@rf))
		{
			stop("Must have recombination fraction data, to convert from class mpcrossMapped to mpcross")
		}
		return(new(to, as(from, "mpcross"), rf = from@rf))
	})
setAs("mpcrossMapped", "mpcrossLG", def = function(from, to)
	{
		groups <- rep(1:length(from@map), times = unlist(lapply(from@map, length)))
		names(groups) <- markers(from)
		allGroups <- unique(groups)
		lg <- new("lg", allGroups = allGroups, groups = groups)
		return(new(to, as(from, "mpcross"), lg = lg, rf = from@rf))
	})
setAs("mpcrossLG", "mpcrossRF", def = function(from, to)
	{
		if(is.null(from@rf) && !is.null(from@lg@imputedTheta) && length(from@lg@allGroups) == 1)
		{
			return(new(to, as(from, "mpcross"), rf = new("rf", theta = from@lg@imputedTheta[[1]])))
		}
		else if(is.null(from@rf))
		{
			stop("As no RF data is present, this object of class mpcrossLG cannot be automatically converted to an object of class mpcrossRF. Please call estimateRF to re-estimate recombination fractions")
		}
		else return(new(to, as(from, "mpcross"), rf = from@rf))
	})
#' @title Create allele encoding corresponding to infinite generations of selfing
#' @description Create allele encoding corresponding to infinite generations of selfing
#' @details In many experiments (particularly those that are significantly inbred), only marker homozygotes are observed, which means that the relationship between marker genotypes and marker alleles is particularly simple. In such cases, generally a marker genotype of some value (say 0) indicates that the individual is homozygous for marker allele 0. 
#' 
#' This function takes in genetic data for the founding lines, genetic data for the final population, and the pedigree. It returns an encoding for marker genotypes where every genotype is homozygous for the marker allele with the same value. 
#' @param founders The genetic data for the founding lines, which are assumed to be inbred
#' @param finals The genetic data for the lines genotyped at the end of the experiment. 
#' @param pedigree The pedigree for the experiment
#' @return An object of class \code{hetData}, which encodes only the marker homozygotes.
#' @examples
#' map <- qtl::sim.map()
#' pedigree <- f2Pedigree(1000)
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
#' #Initially the object contains markers that are fully informative.
#' #The final genetic data contains values 1, 2 and 3, while the genetic data for the founding 
#' #    lines contains only values 1 and 2. 
#' #A value of 1 or 2 in the final genetic data indicates a homozygote for the 
#' #	corresponding marker allele. 
#' #A value of 3 in the final genetic data indicates a heterozygote for the marker allele.
#' #Information about this encoding is stored in the hetData slot.
#' hetData(cross, "D1M1")
#' cross <- cross + biparentalDominant()
#' #Now we have converted all markers to dominant.
#' #The final genetic data contains values 1 and 2, and the genetic data for the founding 
#' #    lines contains only values 1 and 2. 
#' #A value of 2 indicates a homozygote for the corresponding marker allele, OR a 
#' #	marker heterozygote.
#' hetData(cross, "D1M1")
#' #But under infinite generations of selfing, the encoding is simpler. 
#' simpleEncoding <- infiniteSelfing(founders = founders(cross), finals = finals(cross), 
#' 	pedigree = pedigree)
#' simpleEncoding[["D1M1"]]
#' @export
infiniteSelfing <- function(founders, finals, pedigree)
{
	hetData <- lapply(1:ncol(founders), function(x)
		{
			alleles <- unique(founders[,x])
			unname(cbind(alleles,alleles,alleles))
		})
	names(hetData) <- colnames(founders)
	hetData <- new("hetData", hetData)
	return(hetData)
}
#' @title Create heterozygote encodings for SNP markers
#' @description Create encoding which assumes that the single non-homozygote value for a SNP marker is the heterozygote
#' @details This function takes in genotype data for the founding lines and the final poulation. It returns an encoding for hetorozygotes for all markers, where multiallelic markers are assumed to have no heterozygotes. For biallelic markers with three observed alleles in the final population, the extra allele is assumed to be the heterozygote. 
#' @param founders Genetic data for the founding lines of the population
#' @param finals Genetic data for the final genotyped lines of the population
#' @param pedigree Pedigree for the population. Unused by this particular function.
#' @return An object of class \code{hetData}, which contains encodings for the marker heterozygotes and the (unique) marker heterozygote 
#' @export
hetsForSNPMarkers <- function(founders, finals, pedigree)
{
	hetData <- lapply(1:ncol(founders), function(x)
	{
		founderCol <- founders[,x]
		finalCol <- finals[,x]
		uniqueFinals <- unique(finalCol)
		uniqueFounders <- unique(founderCol)
		if(length(uniqueFounders) == 2 && length(uniqueFinals) == 3)
		{
			hetEncoding <- setdiff(uniqueFinals, uniqueFounders)
			matrix <- rbind(cbind(uniqueFounders, uniqueFounders, uniqueFounders), c(uniqueFounders, hetEncoding), c(rev(uniqueFounders), hetEncoding))
		}
		else
		{
			matrix <- cbind(uniqueFounders, uniqueFounders, uniqueFounders)
		}
		dimnames(matrix) <- NULL
		return(matrix)
	})
	names(hetData) <- colnames(founders)
	hetData <- new("hetData", hetData)
	return(hetData)
}
