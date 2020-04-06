#' @title Generate an eight-parent pedigree with improper funnels
#'
#' @description
#' Generate a eight-parent pedigree starting from inbred founders, where the founders in the funnels are not necessarily distinct.
#'
#' @seealso \code{\link{eightParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeRandomFunnels}}, \code{\link{twoParentPedigree}}
#'
#' @param initialPopulationSize The number of initially generated lines, whose genetic material is a mosaic of the eight founding lines. These lines are generated using three generations of structured mating. 
#' @param selfingGenerations The number of selfing generations at the end of the pedigree.
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each initially generated line (if no intercrossing is specified). These lines are then selfed according to selfingGenerations.
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
#' @details
#' Generate a eight-parent pedigree starting from inbred founders. The founders in the funnel for every line are chosen \emph{with replacement}. So for any line from the final population, it is likely that some founding lines are absent from the corresponding funnel, and some appear multiple times. 
#' @examples 
#' pedigree <- eightParentPedigreeImproperFunnels(initialPopulationSize = 10, 
#' 	selfingGenerations = 0, nSeeds = 1)
#' #Generate map
#' map <- qtl::sim.map()
#' #Simulate data
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
#' #Get out a list of funnels, which are rows of this matrix. Note that, of the values 1:8, 
#' #    some are missing within a row, and some are repeated. 
#' getAllFunnels(cross)
#' #convert the pedigree to a graph
#' pedigreeAsGraph <- pedigreeToGraph(pedigree)
#' #Plot it
#' \donttest{plot(pedigreeAsGraph)}
#' #Write it to a file in DOT format
# \donttest{write.graph(graph = pedigreeAsGraph@@graph, format = "dot", file = "./pedigree.dot")}
# This is written in C because otherwise it's just too damn slow (especially for generating the huge populations that we want to use to get numerically accurate results for unit testing)
eightParentPedigreeImproperFunnels <- function(initialPopulationSize, selfingGenerations, nSeeds)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  nonNegativeIntegerArgument(nSeeds)

  return(.Call("eightParentPedigreeImproperFunnels", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), PACKAGE="mpMap2"))
}
eightParentPedigreeImproperFunnelsPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  nonNegativeIntegerArgument(nSeeds)
  initialPopulationSize <- as.integer(initialPopulationSize)
  selfingGenerations <- as.integer(selfingGenerations)
  nSeeds <- as.integer(nSeeds)

  funnels <- sample(1:8, initialPopulationSize*8, replace=TRUE)
  dim(funnels) <- c(initialPopulationSize, 8)

  entries <- 8L + 7L*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize
  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:8] <- father[1:8] <- 0L

  #We have to generate new funnels for every line
  #Create four two-way crosses for funnels
  for (i in 1:initialPopulationSize)
  {
    mother[8L + (i-1)*4+1:4] <- funnels[i, c(1,3,5,7)]
    father[8L + (i-1)*4+1:4] <- funnels[i, c(2,4,6,8)]
  }
  #Create two four-way crosses for each funnel
  mother[8L + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 7L + 2L*(1:(2*initialPopulationSize))
  father[8L + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 8L + 2L*(1:(2*initialPopulationSize))
  #Create a line which mixes all eight founders, for each funnel
  mother[8L + initialPopulationSize*6 + 1:initialPopulationSize] <- 7L + initialPopulationSize*4L + 2L*(1:initialPopulationSize)
  father[8L + initialPopulationSize*6 + 1:initialPopulationSize] <- 8L + initialPopulationSize*4L + 2L*(1:initialPopulationSize)

  currentIndex <- 1L + 8L + initialPopulationSize*5L + initialPopulationSize
  #The next free spot in the pedigree
  nextFree <- currentIndex+initialPopulationSize
  #Now the selfing. 
  #First the case of one generation of selfing
  if(selfingGenerations == 1)
  {
    #The line that we're going to self
    for(lineCounter in currentIndex:(currentIndex+initialPopulationSize-1))
    {
      mother[nextFree:(nextFree+nSeeds-1)] <- father[nextFree:(nextFree+nSeeds-1)] <- lineCounter
      observed[nextFree+nSeeds-1] <- TRUE
      nextFree <- nextFree + nSeeds
    }
  }
  else if(selfingGenerations > 1)
  {
    for(lineCounter in currentIndex:(currentIndex+initialPopulationSize-1))
    {
      #And the number of selfed lines coming off this one
      for(seedCounter in 1:nSeeds)
      {
        father[nextFree:(nextFree+selfingGenerations-1)] <- mother[nextFree:(nextFree+selfingGenerations-1)] <- c(lineCounter, nextFree:(nextFree+selfingGenerations-2))
        observed[nextFree+selfingGenerations-1] <- TRUE
        nextFree <- nextFree + selfingGenerations
      }
    }
  }
  #No selfing
  else
  {
    observed[(nextFree-initialPopulationSize+1):length(observed)] <- TRUE
  }
  if(length(mother) != entries || length(father) != entries || length(observed) != entries || length(lineNames) != entries)
  {
	  stop("Internal error")
  }
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:8L, observed = observed, selfing = "infinite", warnImproperFunnels = FALSE))
}
