% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpcross.R
\name{mpcrossMapped}
\alias{mpcrossMapped}
\title{Create object of class mpcrossMapped}
\usage{
mpcrossMapped(cross, map, rf = NULL)
}
\arguments{
\item{cross}{An object of class \code{mpcross}}

\item{map}{A genetic map, formatted as in the \code{qtl} package.}

\item{rf}{Optional recombination fraction data. Leave as NULL if there is no such data.}
}
\value{
An object of class \code{mpcrossMapped}, constructed from the arguments.
}
\description{
Create object of class mpcrossMapped
}
\details{
This function constructs an object of class \code{mpcrossMapped} representing a multi-parent population with a map. It takes in an object of class \code{mpcross}, a genetic map, and optional recombination fraction data.
}
