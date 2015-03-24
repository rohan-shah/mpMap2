#' Generate a two-parent pedigree which starts from inbred founders
#'
#' @param initialPopulationSize The number of F1 lines generated
#' @param selfingGenerations The number of selfing generations at the end of the pedigree
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each F1 if no intercrossing is specified. These lines are then selfed according to selfingGenerations
#' @param intercrossingGenerations The number of generations of random mating performed from the F1 generation. Population size is maintained at that specified by initialPopulationSize
twoParentPedigree <- function(initialPopulationSize, selfingGenerations, nSeeds, intercrossingGenerations)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  nonNegativeIntegerArgument(nSeeds)
  nonNegativeIntegerArgument(intercrossingGenerations)

  if(initialPopulationSize <= 2 && intercrossingGenerations > 0)
  {
    stop("Random mating is impossible with only two lines per generation")
    #....and more importantly it means that the sample command below gets screwed up, because we're calling sample(x) where length(x) == 1, which samples from 1:x
  }
  
  entries <- 2L + initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize

  mother <- father <- rep(NA, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:2] <- father[1:2] <- 0L

  mother[3:(2+initialPopulationSize)] <- 1L
  father[3:(2+initialPopulationSize)] <- 2L

  currentIndex <- 3L
  if(intercrossingGenerations > 0 && initialPopulationSize > 2)
  {
    lastGenerationStart <- 3L
    lastGenerationEnd <- 2L+initialPopulationSize
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
        nextFree <- nextFree + selfingGenerations
      }
    }
  }
  return(new("pedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:2L))
}