#' @title Generate a four-parent pedigree
#'
#' @description
#' Generate a four-parent pedigree starting from inbred founders, using a single funnel
#'
#' @details
#' Note that unlike \code{\link{fourParentPedigreeRandomFunnels}}, there is no intercrossing allowed in the single funnel case because the relevant haplotype probabilities assume randomly chosen funnels
#' @seealso \code{\link{fourParentPedigreeRandomFunnels}}, \code{\link{twoParentPedigree}}
#' @param initialPopulationSize The number of F1 lines generated
#' @param selfingGenerations The number of selfing generations at the end of the pedigree
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each F1 if no intercrossing is specified. These lines are then selfed according to selfingGenerations
#' @param intercrossingGenerations The number of generations of intercrossing, after each F2 line is generated. 
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
fourParentPedigreeSingleFunnel <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  positiveIntegerArgument(nSeeds)
  nonNegativeIntegerArgument(intercrossingGenerations)

  if(initialPopulationSize <= 2 && intercrossingGenerations > 0)
  {
    stop("Random mating is impossible with only two lines per generation")
    #....and more importantly it means that the sample command below gets screwed up, because we're calling sample(x) where length(x) == 1, which samples from 1:x
  }
  return(.Call("fourParentPedigreeSingleFunnel", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), as.integer(intercrossingGenerations), PACKAGE="mpMap2"))
}
fourParentPedigreeSingleFunnelPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  positiveIntegerArgument(nSeeds)
  intercrossingGenerations <- as.integer(intercrossingGenerations)
  initialPopulationSize <- as.integer(initialPopulationSize)

  if(initialPopulationSize <= 2 && intercrossingGenerations > 0)
  {
    stop("Random mating is impossible with only two lines per generation")
#....and more importantly it means that the sample command below gets screwed up, because we're calling sample(x) where length(x) == 1, which samples from 1:x
  }

  #Four founders, + 6 pairs + intialPopulationSize funnels
  entries <- 4L + 2L + initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize

  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:4] <- father[1:4] <- 0L

  #Crosses of the founders
  mother[5:6] <- c(1L, 3L)
  father[5:6] <- c(2L, 4L)
  #line 5 = 1,2
  #line 6 = 3,4
  mother[6+1:initialPopulationSize] <- 5L
  father[6+1:initialPopulationSize] <- 6L

  currentIndex <- 7L
  nextFree <- currentIndex + initialPopulationSize
  if(intercrossingGenerations > 0)
  {
    for(i in 1:intercrossingGenerations)
    {
      for(lineCounter in nextFree:(nextFree+initialPopulationSize-1))
      {
        mother[lineCounter] <- lineCounter - initialPopulationSize
        father[lineCounter] <- sample(setdiff((nextFree-initialPopulationSize):(nextFree-1), lineCounter - initialPopulationSize), 1)
      }
      currentIndex <- currentIndex + initialPopulationSize
      nextFree <- nextFree + initialPopulationSize
    }
  }

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
    observed[(nextFree-initialPopulationSize):(nextFree-1)] <- TRUE
  }
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:4L, observed = observed, selfing = "infinite", warnImproperFunnels = TRUE))
}
