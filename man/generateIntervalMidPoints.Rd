% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateIntervalMidPoints.R
\name{generateIntervalMidPoints}
\alias{generateIntervalMidPoints}
\title{Specify interval midpoints}
\usage{
generateIntervalMidPoints(object)
}
\arguments{
\item{object}{The object of class \code{mpcrossMapped} from which to take the interval midpoints.}
}
\value{
A function which can be applied to an object of class \code{mpcrossMapped} by \code{imputeFounders} or \code{computeGenotypeProbabilities}.
}
\description{
Specify interval midpoints
}
\details{
Some functions, such as \code{imputeFounders} and \code{computeGenotypeProbabilities}, take in a set of genetic positions as one of the inputs. This function is an easy way to specify the midpoint of every marker interval. 

Note that you don't have to explicitly evaluate this function, it can be passed in directly (see examples).
}
\examples{
data(simulatedFourParentData)
#Create object that includes the correct map
mapped <- new("mpcrossMapped", simulatedFourParentData, map = simulatedFourParentMap)
#Estimate IBD genotypes at all the markers, and marker midpoints
imputed <- imputeFounders(mapped, errorProb = 0.02, 
	extraPositions = generateIntervalMidPoints(mapped))
#Alternatively we can explicitly evaluate the function. This is identical to above.
imputed <- imputeFounders(mapped, errorProb = 0.02, 
	extraPositions = generateIntervalMidPoints)
}
