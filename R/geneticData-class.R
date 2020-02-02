#' @include hetData-class.R
#' @include pedigree-class.R
#' @include map-class.R
#' @include canSkipValidity.R
checkGeneticData <- function(object)
{
	errors <- c()
	if(!is.numeric(object@founders))
	{
		errors <- c(errors, "Slot founders must be a numeric matrix")
	}
	#Error if all founders are NA (because the next condition after this involves max(,na.rm=T), which requires at least one non-NA val)
	else if(all(is.na(object@founders)))
	{
		errors <- c(errors, "Slot founders cannot contain only NA")
	}
	else if(max(abs(round(object@founders) - object@founders), na.rm=TRUE) > 0)
	{
		errors <- c(errors, "Slot founders must contain integer values")
	}
	
	if(!is.numeric(object@finals))
	{
		errors <- c(errors, "Slot finals must be a numeric matrix")
	}
	#Don't do the next two checks if there is no data.
	else if(length(object@finals) > 0)
	{
		if(all(is.na(object@finals)))
		{
			errors <- c(errors, "Slot finals cannot contain only NA")
		}
		else if(max(abs(round(object@finals) - object@finals), na.rm=TRUE) > 0)
		{
			errors <- c(errors, "Slot finals must contain integer values")
		}
	}

	if(length(dim(object@finals)) != 2)
	{
		errors <- c(errors, "Slot finals must be a matrix")
	}
	if(length(dim(object@founders)) != 2)
	{
		errors <- c(errors, "Slot founders must be a matrix")
	}
	#Check that row and column names of finals and founders exist
	if(is.null(dimnames(object@founders)) || any(unlist(lapply(dimnames(object@founders), is.null))))
	{
		errors <- c(errors, "Slot founders must have row and column names")
	}
	#Don't check for row and column names if there is no actual data.
	if(length(object@finals) > 0 && (is.null(dimnames(object@finals)) || any(unlist(lapply(dimnames(object@finals), is.null)))))
	{
		errors <- c(errors, "Slot finals must have row and column names")
	}

	#Check for NA's in row and column names of finals and founders
	if(any(unlist(lapply(dimnames(object@founders), function(x) any(is.na(x))))))
	{
		errors <- c(errors, "Row and column names of slot founders cannot be NA")
	}
	#Skip this if there is no data (this avoids a warnings)
	if(length(object@finals) > 0 && any(unlist(lapply(dimnames(object@finals), function(x) any(is.na(x))))))
	{
		errors <- c(errors, "Row and column names of slot finals cannot be NA")
	}

	#Check that row and column names of finals and founders are unique
	if(any(unlist(lapply(dimnames(object@finals), function(x) length(x) != length(unique(x))))))
	{
		errors <- c(errors, "Row and column names of slot finals cannot contain duplicates")
	}
	if(any(unlist(lapply(dimnames(object@founders), function(x) length(x) != length(unique(x))))))
	{
		errors <- c(errors, "Row and column names of slot founders cannot contain duplicates")
	}

	nMarkers <- ncol(object@founders)
	markers <- colnames(object@finals)
	if(ncol(object@finals) != nMarkers || length(object@hetData) != nMarkers)
	{
		errors <- c(errors, "Slots finals, founders and hetData had different numbers of markers")
		#If they have the wronh dimensions then the colnames / names checks have the wrong lengths. So just return from here. 
		return(errors)
	}
	if(any(colnames(object@founders) != markers))
	{
		errors <- c(errors, "Slot finals must have the same colnames as slot founders")
	}
	if(any(names(object@hetData) != markers))
	{
		errors <- c(errors, "Slot hetData refers to different markers to slot finals")
	}
	
	#Check that each founder and final is in the pedigree
	if(!all(rownames(object@founders) %in% object@pedigree@lineNames))
	{
		errors <- c(errors, "Not all founder lines were named in the pedigree")
	}
	if(!all(rownames(object@finals) %in% object@pedigree@lineNames))
	{
		errors <- c(errors, "Not all final lines were named in the pedigree")
	}

	#If we have information on the founders in the pedigree, check that the founders ARE in fact those lines
	if(inherits(object@pedigree, "detailedPedigree"))
	{
		if(!all(rownames(object@founders) %in% object@pedigree@lineNames[object@pedigree@initial]) || nrow(object@founders) != length(object@pedigree@initial))
		{
			errors <- c(errors, "Founder lines did not match those specified in the pedigree")
		}
		#Only check this if there are actually final lines
		if(length(object@finals) > 0)
		{
			if(!all(rownames(object@finals) %in% object@pedigree@lineNames[object@pedigree@observed]) || nrow(object@finals) != sum(object@pedigree@observed))
			{
				errors <- c(errors, "Final lines did not match those specified in the pedigree")
			}
		}
	}
	#Don't continue to the C code if there are already problems - Especially if the founders / finals slots have the wrong types
	if(length(errors) > 0) return(errors)
	#This checks the relation between the het data, founder data and final data. It doesn't check that the het data is itself valid
	#It also checks that if any of the founders are NULL for a marker, ALL the founder alleles must be NULL, and all the finals alleles must be NULL too
	alleleDataErrors <- .Call("alleleDataErrors", object, 5, PACKAGE="mpMap2")
	if(length(alleleDataErrors) > 0)
	{
		errors <- c(errors, alleleDataErrors)
	}
	if(length(errors) > 0) return(errors)

	#Check imputed slot
	if(!is.null(object@imputed))
	{
		if(!identical(rownames(imputationData(object)), rownames(object@finals)))
		{
			return("Row names of slot imputed@data must be the same as those of slot finals")
		}
		errors <- validObject(object@imputed)
		if(length(errors) > 0) return(errors)
	}
	#Check probabilities
	if(!is.null(object@probabilities))
	{
		if(!is.numeric(object@probabilities@data))
		{
			return("Slot probabilities@data must be a numeric matrix")
		}
		nGenotypes <- length(unique(object@probabilities@key[,3]))
		if(nrow(object@probabilities@data) != nLines(object) * nGenotypes)
		{
			return("Number of rows of probabilities@data must be consistent with probabilities@key and nrow(finals)")
		}
		errors <- validObject(object@probabilities)
		if(length(errors) > 0) return(errors)
	}

	#Check pheno data
	if(!is.null(object@pheno))
	{
		if(ncol(object@pheno) == 0 || nrow(object@pheno) == 0)
		{
			return("A phenotype dataset with zero columns or rows is not allowed. Please enter NULL if there is no phenotype data")
		}
		if(nrow(object@pheno) != nrow(object@finals))
		{
			return("The number of rows in slot pheno must match the number of rows in slot finals")
		}
		if(!isTRUE(all.equal(rownames(object@pheno), rownames(object@finals), check.attributes = FALSE)))
		{
			return("The row names of slot pheno must match the row names of slot finals")
		}
	}
	return(TRUE)
}
checkImputedData <- function(object)
{
	if(!is.numeric(object@data))
	{
		return("Slot data must be an integer matrix")
	}
	if(!is.numeric(object@key))
	{
		return("Slot key must be an integer matrix")
	}
	if(ncol(object@key) != 3L)
	{
		return("Slot key must have three columns")
	}
	allMapMarkers <- unlist(lapply(object@map, names))
	if(!isTRUE(all.equal(allMapMarkers, colnames(object@data), check.attributes = FALSE)))
	{
		return("Slot data must have marker names that match the markers in slot map")
	}
	if(!.Call("checkImputedBounds", object, PACKAGE="mpMap2"))
	{
		return("Slot imputed@data must contain values in imputed@key")
	}
	tmp <- unlist(lapply(object@map, names))
	names(tmp) <- NULL
	if(!isTRUE(all.equal(colnames(object@data), tmp, check.attributes = FALSE)))
	{
		return("Column names of imputed object did not match the associated map")
	}
	return(TRUE)
}
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
.imputed <- setClass("imputed", slots=list(data = "matrix", key = "matrix", map = "map", errors = "matrixOrNULL"), validity = checkImputedData)
setClassUnion("imputedOrNULL", c("imputed", "NULL"))
checkProbabilities <- function(object)
{
	if(!is.numeric(object@data))
	{
		return("Slot data must be an integer matrix")
	}
	if(!is.numeric(object@key))
	{
		return("Slot key must be an integer matrix")
	}
	if(ncol(object@key) != 3L)
	{
		return("Slot key must have three columns")
	}
	tmp <- unlist(lapply(object@map, names))
	names(tmp) <- NULL
	if(!identical(colnames(object@data), tmp))
	{
		return("Column names of probabilities object did not match the associated map")
	}
	return(TRUE)
}
#' @title Identity-by-descent genotype probabilities
#' @description Identity-by-descent genotype probabilities
#' @details This object contains the identify-by-descent genotype probabilities, as computed by \code{\link{computeGenotypeProbabilities}}. The slot \code{data} is a numeric matrix containing the actual computations, where columns correspond to genetic positions. 
#' 
#' Describing the rows of the data matrix is more complicated. The slot \code{key} is a matrix containing three columns, the first two being founder alleles, and the third being an encoding of that combination. If \code{k} is the number of rows in \code{key}, then the first \code{k} rows of the data matrix correspond to the first genetic line in the population. Specifically, the first row corresponds to genotype probabilities for the first line, for the combination of founder alleles encoding as 1. The second corresponds to genotype probabilities for the first line, for the combination encoded as 2, etc. 
#' @param data An integer matrix containing the IBD genotype probabilities. Rows correspond to combinations of genetic lines and founder lines, and columns correspond to genetic positions. 
#' @param key A matrix identifying how pairs of founder alleles are mapped to rows of the data slot. 
#' @param map The map of positions at which IBD genotype probabilities are computed. 
.probabilities <- setClass("probabilities", slots=list(data = "matrix", key = "matrix", map = "map"), validity = checkProbabilities)
setClassUnion("probabilitiesOrNULL", c("probabilities", "NULL"))
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
#' @title Object containing the genetic data for a population
#' @description Object containing the genetic data for a population
#' @details This object contians the genetic data for a population. Required data includes the genetic data for the founding lines of the poulation, the final lines of the population, information about the enoding of heterozygotes, and the pedigree used to generate the final genetic lines from the founding genetic line. 
#' 
#' Optional data includes IBD genotype imputations, a data.frame of phenotypes, and IBD genotype probabilities.
#' 
#' This class has extensive validity checking, to ensure that all the different inputs are compatible and meet the requirements. If an error is found, an informative error message should be produced. 
#' @slot founders The genetic data for the founding lines of the population. Must be an integer matrix, where rows correspond to genetic lines and columns correspond to genetic markers.
#' @slot finals The genetic data for the final lines of the population. Must be an integer matrix, where rows correspond to genetic lines and columns correspond to genetic markers.
#' @slot hetData Information about the encoding of marker heterozygotes.
#' @slot pedigree Object of class \code{pedigree} with information about how the final genetic lines are generated from the founding lines. 
#' @slot imputed Optional data about imputed IBD genotypes. Can be generated using \code{\link{imputeFounders}}, assuming there is a genetic map available.
#' @slot probabilities Optional data about IBD genotype probabilities. Can be generated using \code{\link{computeGenotypeProbabilities}}, assuming there is a genetic map available.
#' @slot pheno Optional 
.geneticData <- setClass("geneticData", slots=list(finals = "matrix", founders = "matrix", hetData = "hetData", pedigree = "pedigree", imputed = "imputedOrNULL", pheno = "data.frameOrNULL", probabilities = "probabilitiesOrNULL"), validity = checkGeneticData, contains = "canSkipValidity")
checkGeneticDataList <- function(object)
{
	if(any(unlist(lapply(object, class)) != "geneticData"))
	{
		return("An object of class geneticDataList must contain objects of class geneticData")
	}
	return(TRUE)
}
.geneticDataList <- setClass("geneticDataList", "list", validity = checkGeneticDataList)
#' @title Initialize method which can skip the validity check
#' @details Initialize method which can skip the validity check
#' @description This is an initialization method with an optional \code{skipValidity} argument. If this argument is set to \code{TRUE}, the validity check is skipped. This is used by some internal functions within the package, as the validity check can be slow, and internal code is (presumably) guaranteed to produce valid objects. 
#' @param .Object the object to initialize
#' @param ... Other arguments. Only \code{skipValidity} is used.
#' @rdname initialize
setMethod(f = "initialize", signature = "geneticDataList", definition = canSkipValidityInitialize)
