% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneticData-class.R
\docType{class}
\name{geneticData-class}
\alias{geneticData-class}
\alias{.geneticData}
\title{Object containing the genetic data for a population}
\description{
Object containing the genetic data for a population
}
\details{
This object contians the genetic data for a population. Required data includes the genetic data for the founding lines of the poulation, the final lines of the population, information about the enoding of heterozygotes, and the pedigree used to generate the final genetic lines from the founding genetic line. 

Optional data includes IBD genotype imputations, a data.frame of phenotypes, and IBD genotype probabilities.

This class has extensive validity checking, to ensure that all the different inputs are compatible and meet the requirements. If an error is found, an informative error message should be produced.
}
\section{Slots}{

\describe{
\item{\code{founders}}{The genetic data for the founding lines of the population. Must be an integer matrix, where rows correspond to genetic lines and columns correspond to genetic markers.}

\item{\code{finals}}{The genetic data for the final lines of the population. Must be an integer matrix, where rows correspond to genetic lines and columns correspond to genetic markers.}

\item{\code{hetData}}{Information about the encoding of marker heterozygotes.}

\item{\code{pedigree}}{Object of class \code{pedigree} with information about how the final genetic lines are generated from the founding lines.}

\item{\code{imputed}}{Optional data about imputed IBD genotypes. Can be generated using \code{\link{imputeFounders}}, assuming there is a genetic map available.}

\item{\code{probabilities}}{Optional data about IBD genotype probabilities. Can be generated using \code{\link{computeGenotypeProbabilities}}, assuming there is a genetic map available.}

\item{\code{pheno}}{Optional}
}}

