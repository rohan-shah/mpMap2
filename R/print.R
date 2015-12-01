#' @include mpcross-class.R
print.mpcrossLG <- function(x)
{
  cat("This mpcross object contains linkage groups\n")
  callNextMethod()
}
setMethod(f = "print", signature = "mpcrossLG", definition = print.mpcrossLG)
print.mpcrossRF <- function(x)
{
  cat("This mpcross object contains recombination fractions\n\n")
  callNextMethod()
}
setMethod(f = "print", signature = "mpcrossRF", definition = print.mpcrossRF)
print.mpcross <- function(x)
{
  nGeneticDatasets <- length(x@geneticData)
  nLines <- nLines(x)
  if(nGeneticDatasets == 1)
  {
    missingInFounders <- apply(x@geneticData[[1]]@founders, 2, function(y) any(is.na(y)))
    founderAlleles <- apply(x@geneticData[[1]]@founders, 2, function(y) length(unique(y)))
    proportionMissing <- apply(x@geneticData[[1]]@finals, 2, function(y) sum(is.na(y)) / nLines)
    cat("-------------------------------------------------------\n")
    cat("Summary of mpcross object\n")
    cat("-------------------------------------------------------\n")
    cat(sum(missingInFounders), " markers had missing values in founders\n")
    cat(sum(founderAlleles==1), " markers had non-polymorphic founder genotypes\n")
    cat("-------------------------------------------------------\n")
    cat(sum(founderAlleles==2), " markers were biallelic.\n")
    cat(sum(founderAlleles > 2), " markers were multiallelic.\n")
    cat("-------------------------------------------------------\n")
    cat(sum(proportionMissing>.05), " markers had >5% missing data.\n")
    cat(sum(proportionMissing>.10), " markers had >10% missing data.\n")
    cat(sum(proportionMissing>.20), " markers had >20% missing data.\n")
  } else
  {
    cat("This mpcross object contains", nGeneticDatasets, " genetic data sets\n")
    for(i in 1:nGeneticDatasets)
    {
      missingInFounders <- apply(x@geneticData[[i]]@founders, 2, function(y) any(is.na(y)))
      founderAlleles <- apply(x@geneticData[[i]]@founders, 2, function(y) length(unique(y)))
      proportionMissing <- apply(x@geneticData[[i]]@finals, 2, function(y) sum(is.na(y)) / nLines[i])
      cat("\n-------------------------------------------------------\n")
      cat("Summary of data set", i, "\n")
      cat("-------------------------------------------------------\n")
      cat(sum(missingInFounders), " markers had missing values in founders\n")
      cat(sum(founderAlleles==1), " markers had non-polymorphic founder genotypes\n")
      cat("-------------------------------------------------------\n")
      cat(sum(founderAlleles==2), " markers were biallelic.\n")
      cat(sum(founderAlleles > 2), " markers were multiallelic.\n")
      cat("-------------------------------------------------------\n")
      cat(sum(proportionMissing>.05), " markers had >5% missing data.\n")
      cat(sum(proportionMissing>.10), " markers had >10% missing data.\n")
      cat(sum(proportionMissing>.20), " markers had >20% missing data.\n")
    }
  }
}
setMethod(f = "print", signature = "mpcross", definition = print.mpcross)
