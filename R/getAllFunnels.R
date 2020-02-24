#' @title Get funnels
#' @rdname getAllFunnels
#' @description Get the order of the founding lines, as they contribute to each line in the final population
#' @details In multi-parent experimental designs, the founding lines of the population are combined together through initial mixing generations. For experiments without further intercrossing generations, the order in which these mixing crosses occur influences the genotypes of the final lines. It can be important to examine or visualise these orders, which are known as funnels. 
#'
#' This function returns a matrix, where each row corresponds to a genetic line in the final population, and each column corresponds to a position in the mixing step. So if a row of the returned matrix contains the values 4, 1, 2, 3, then the pedigee that generated the first individual in the experiment started by crossing founders 4 and 1 to give individual 41, and 2 and 3 to give individual 23. Then individuals 41 and 23 are crossed to generate individual 4123, which after inbreeding results in the first final genetic line. 
#' 
#' If sex is considered to be unimportant, then many orderings are equivalent. For example, the ordering 4, 1, 2, 3 of the initial founders is equivalent to 1, 4, 2, 3. In this case each funnel can be put into a standardised ordering, by setting \code{standardised} to \code{FALSE}. 
#'
#' Note that if there are generations of random interbreeding in the population (often referred to as maintenance generations), then there is no "funnel" associated with a genetic line, and values of NA are returned. In that case, see \code{\link{getAllFunnelsIncAIC}}.  
#'
#' Note that funnels for all pedigrees simulated by mpMap2 are already standardised. This will not generally be the case for realy experiments. 
#' @param cross The object of class \code{mpcross} containing the pedigree of interest
#' @param standardised Should the output funnels be standardised?
#' @return An integer matrix with rows representing genetic lines, and columns representing positions within the funnel. 
#' @examples
#' data(simulatedFourParentData)
#' #Funnels used to generate the first ten lines
#' #Because this is simulated data, they are already standardised, 
#' #' with the first founder in the first position in the mixing step. 
#' getAllFunnels(simulatedFourParentData)[1:10, ]
#' @export 
getAllFunnels <- function(cross, standardised = FALSE)
{
	if(!is.logical(standardised) || length(standardised) != 1)
	{
		stop("Input standardised must be TRUE or FALSE")
	}
	if(class(cross) == "geneticData")
	{
		return(.Call("getAllFunnels", cross, standardised, PACKAGE="mpMap2"))
	}
	else if(inherits(cross, "mpcross"))
	{
		if(length(cross@geneticData) == 1)
		{
			return(.Call("getAllFunnels", cross@geneticData[[1]], standardised, PACKAGE="mpMap2"))
		}
		else
		{
			return(lapply(cross@geneticData, function(x) .Call("getAllFunnels", x, standardised, PACKAGE="mpMap2")))
		}
	}
	else
	{
		stop("Input must be of class geneticData or mpcross")
	}
}
#' @title Get all funnels, including AIC lines
#' @description Get every order of the founding lines, which makes a contribution to the final population
#' @rdname getAllFunnelsIncAIC
#' @details This function is similar to \code{\link{getAllFunnels}}, but more useful for populations with maintenance (or AIC) generations. It returns a list of all the mixing orders in the initial generations, which make a genetic contribution to the final population. Unlike for \code{\link{getAllFunnels}}, rows of the returned matrix DO NOT refer to specific genetic lines. 
#' @param cross The object of class \code{mpcross} containing the pedigree of interest
#' @param standardised Should the output funnels be standardised?
#' @return Matrix of mixing orders that contribute to the final popluation. Rows DO NOT refer to specific genetic lines. 
#' @examples
#' set.seed(1)
#' pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 1000,
#'      selfingGenerations = 6, intercrossingGenerations = 1)
#' #Assume infinite generations of selfing in subsequent analysis
#' selfing(pedigree) <- "infinite"
#' #Generate random map
#' map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x = FALSE)
#' #Simulate data
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1L)
#' #Because we have maintenance in this experiment, we can't get out the funnels per genetic line
#' funnels <- getAllFunnels(cross)
#' dim(funnels)
#' funnels[1:10,]
#' #But we can get out a list of all the funnels that go into the experiment.
#' funnels <- getAllFunnelsIncAIC(cross)
#' dim(funnels)
#' funnels[1:10,]
#' @export
getAllFunnelsIncAIC <- function(cross, standardised = FALSE)
{
	if(!is.logical(standardised) || length(standardised) != 1)
	{
		stop("Input standardised must be TRUE or FALSE")
	}
	if(class(cross) == "geneticData")
	{
		return(.Call("getAllFunnels", cross, standardised, PACKAGE="mpMap2"))
	}
	else if(inherits(cross, "mpcross"))
	{
		if(length(cross@geneticData) == 1)
		{
			return(.Call("getAllFunnelsIncAIC", cross@geneticData[[1]], standardised, PACKAGE="mpMap2"))
		}
		else
		{
			return(lapply(cross@geneticData, function(x) .Call("getAllFunnelsIncAIC", x, standardised, PACKAGE="mpMap2")))
		}
	}
	else
	{
		stop("Input must be of class geneticData or mpcross")
	}
}
