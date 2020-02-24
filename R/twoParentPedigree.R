#' @title Generate a two-parent pedigree which starts from inbred founders
#'
#' @description
#' Generate a two-parent pedigree starting from inbred founders
#'
#' @param initialPopulationSize The number of F1 lines generated
#' @param selfingGenerations The number of selfing generations at the end of the pedigree
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each F1 if no intercrossing is specified. These lines are then selfed according to selfingGenerations
#' @param intercrossingGenerations The number of generations of random mating performed from the F1 generation. Population size is maintained at that specified by initialPopulationSize
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @examples
#' plotWOptions <- function(graph) 
#' 	plot(graph, vertex.size = 8, vertex.label.cex=0.6, edge.arrow.size=0.01, edge.width=0.2)
#' #F2 design
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 1, 
#' 	intercrossingGenerations = 0, nSeeds = 1)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#'
#' #An equivalent F2 design (if the founders really are inbred)
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 0, 
#' 	intercrossingGenerations = 1, nSeeds = 0)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#'
#' #Another equivalent F2 design (if the founders really are inbred)
#' pedigree <- twoParentPedigree(initialPopulationSize = 1, selfingGenerations = 1, 
#' 	intercrossingGenerations = 0, nSeeds=10)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#' 
#' #A RIL design (10 generations of inbreeding)
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 10, 
#' 	intercrossingGenerations = 0, nSeeds = 1)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#'
#' #Another RIL design (10 generations of inbreeding)
#' pedigree <- twoParentPedigree(initialPopulationSize = 1, selfingGenerations = 10, 
#' 	intercrossingGenerations = 0, nSeeds = 10)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#
#' #One generation of mixing followed by 10 generations of inbreeding
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 10, 
#' 	intercrossingGenerations = 1, nSeeds = 1)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#'
#' #Two generations of mixing and no inbreeding
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 0, 
#' 	intercrossingGenerations = 2, nSeeds = 0)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#'
#' #One generation of mixing, and then two selfed lines are generated (10 generations of selfing)
#' pedigree <- twoParentPedigree(initialPopulationSize = 10, selfingGenerations = 10, 
#' 	intercrossingGenerations = 1, nSeeds = 2)
#' graph <- pedigreeToGraph(pedigree)
#' plotWOptions(graph)
#' @export
twoParentPedigree <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:2] <- father[1:2] <- 0L

  mother[3:(2+initialPopulationSize)] <- 1L
  father[3:(2+initialPopulationSize)] <- 2L

  currentIndex <- 3L
  #If intercrossingGenerations == 0 skip this
  if(intercrossingGenerations > 0)
  {
    lastGenerationStart <- currentIndex
    lastGenerationEnd <- currentIndex-1+initialPopulationSize
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
      observed[nextFree:(nextFree+nSeeds-1)] <- TRUE
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
      observed[3:length(observed)] <- TRUE
    }
    #If there's no selfing but there is intercrossing then mark the last set of intercrossing lines as observed
    else
    {
      observed[lastGenerationStart:lastGenerationEnd] <- TRUE
    }
  }
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:2L, observed = observed, selfing = "infinite", warnImproperFunnels = TRUE))
}
