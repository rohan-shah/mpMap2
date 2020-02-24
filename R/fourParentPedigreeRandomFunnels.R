#' @title Generate a four-parent pedigree
#'
#' @description
#' Generate a four-parent pedigree starting from inbred founders, using a random funnel
#'
#' @seealso \code{\link{fourParentPedigreeSingleFunnel}}, \code{\link{twoParentPedigree}}
#'
#' @param initialPopulationSize The number of F1 lines generated
#' @param selfingGenerations The number of selfing generations at the end of the pedigree
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each F1 if no intercrossing is specified. These lines are then selfed according to selfingGenerations
#' @param intercrossingGenerations The number of generations of random mating performed from the F1 generation. Population size is maintained at that specified by initialPopulationSize
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
# This is written in C because otherwise it's just too damn slow (especially for generating the huge populations that we want to use to get numerically accurate results for unit testing)
fourParentPedigreeRandomFunnels <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  return(.Call("fourParentPedigreeRandomFunnels", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), as.integer(intercrossingGenerations), PACKAGE="mpMap2"))
}
fourParentPedigreeRandomFunnelsPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  positiveIntegerArgument(nSeeds)
  nonNegativeIntegerArgument(intercrossingGenerations)
  intercrossingGenerations <- as.integer(intercrossingGenerations)
  initialPopulationSize <- as.integer(initialPopulationSize)

  if(initialPopulationSize <= 2 && intercrossingGenerations > 0)
  {
    stop("Random mating is impossible with only two lines per generation")
    #....and more importantly it means that the sample command below gets screwed up, because we're calling sample(x) where length(x) == 1, which samples from 1:x
  }
  
  #Four founders, + 6 pairs + intialPopulationSize funnels
  entries <- 4L + 6L + initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize

  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:4] <- father[1:4] <- 0L

  #Crosses of the founders
  mother[5:10] <- c(1L, 3L, 1L, 2L, 1L, 2L)
  father[5:10] <- c(2L, 4L, 3L, 4L, 4L, 3L)
  #line 5 = 1,2
  #line 6 = 3,4
  #line 7 = 1,3
  #line 8 = 2,4
  #line 9 = 1,4
  #line10 = 2,3

  #The three different funnels:
  #lines 5, 6
  #lines 7, 8
  #lines 9, 10
  funnelChoices <- sample(1:3, replace = TRUE, initialPopulationSize)
  mother[10+1:initialPopulationSize] <- 4L + 2L*funnelChoices - 1L
  father[10+1:initialPopulationSize] <- 4L + 2L*funnelChoices

  #Now intercrossing generations
  currentIndex <- 11L
  if(intercrossingGenerations > 0)
  {
    lastGenerationStart <- currentIndex
    lastGenerationEnd <- currentIndex-1L+initialPopulationSize
    for(i in 1:intercrossingGenerations)
    {
      for(lineCounter in lastGenerationStart:lastGenerationEnd)
      {
        mother[lineCounter + initialPopulationSize] <- lineCounter
        father[lineCounter + initialPopulationSize] <- sample(setdiff(lastGenerationStart:lastGenerationEnd, lineCounter), 1)
      }
      lastGenerationStart <- lastGenerationStart + initialPopulationSize
      lastGenerationEnd <- lastGenerationEnd + initialPopulationSize
    }
    currentIndex <- lastGenerationStart
  }

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
    #...and no intercrossing
    if(intercrossingGenerations == 0)
    {
      observed[11:length(observed)] <- TRUE
    }
    #If there's no selfing but there is intercrossing then mark the last set of intercrossing lines as observed
    else
    {
      observed[lastGenerationStart:lastGenerationEnd] <- TRUE
    }
  }
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:4L, observed = observed, selfing = "infinite", warnImproperFunnels = TRUE))
}
