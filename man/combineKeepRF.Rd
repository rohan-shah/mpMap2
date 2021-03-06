% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineKeepRF.R
\name{combineKeepRF}
\alias{combineKeepRF}
\title{Combine mpcross objects, keeping recombination fraction data}
\usage{
combineKeepRF(
  object1,
  object2,
  verbose = TRUE,
  gbLimit = -1,
  callEstimateRF = TRUE,
  skipValidity = FALSE
)
}
\arguments{
\item{object1}{An object of class \code{mpcrossRF}}

\item{object2}{Another object of class \code{mpcrossRF}}

\item{verbose}{Passed straight through to estimateRF}

\item{gbLimit}{Passed straight through to estimateRF}

\item{callEstimateRF}{Should \code{estimateRF} be called, to compute any missing estimates?}

\item{skipValidity}{Should we skip the validity check for object construction, in this function? Running the validity checks can be expensive, and in theory internal package code is trusted to generate valid objects.}
}
\value{
A new object of class \code{mpcrossRF} containing the combined information of the two input objects.
}
\description{
Combine mpcross objects, keeping recombination fraction data
}
\details{
This function takes two objects containing disjoint sets of markers, each containing estimated recombination fractions for their individual sets of markers. A new object is returned that contains the combined set of markers, and also contains recombination fraction data. 

This function is more efficient than other ways of achieving this, as it keeps the recombination fraction data contained in the original objects. If \code{callEstimateRF} is \code{TRUE}, it also computes the missing recombination fraction estimates between markers in different objects, using a call to \code{estimateRF}.
}
