#' @title Generate an eight-parent pedigree, using a single funnel
#'
#' @description
#' Generate a eight-parent pedigree starting from inbred founders, using a single funnel. 
#'
#' @seealso \code{\link{eightParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeRandomFunnels}}, \code{\link{twoParentPedigree}}
#'
#' @param initialPopulationSize The number of initially generated lines, whose genetic material is a mosaic of the eight founding lines. These lines are generated using three generations of structured mating. 
#' @param selfingGenerations The number of selfing generations at the end of the pedigree. 
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each initially generated line (if no intercrossing is specified). These lines are then selfed according to selfingGenerations.
#' @param intercrossingGenerations The number of generations of random mating performed from the F1 generation. Population size is maintained at that specified by initialPopulationSize. 
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
#' @examples 
#' pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 10, 
#' 	selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 1)
#' map <- qtl::sim.map()
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
#' #Get out a list of funnels, which are rows of this matrix. For this pedigree, every funnel is 1:8. 
#' getAllFunnels(cross)
#' #convert the pedigree to a graph
#' pedigreeAsGraph <- pedigreeToGraph(pedigree)
#' #Plot it
#' \dontrun{plot(pedigreeAsGraph)}
#' #Write it to a file in DOT format
#' \dontrun{write.graph(graph = pedigreeAsGraph@@graph, format = "dot", file = "./pedigree.dot")}
# This is written in C because otherwise it's just too damn slow (especially for generating the huge populations that we want to use to get numerically accurate results for unit testing)
eightParentPedigreeSingleFunnel <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  return(.Call("eightParentPedigreeSingleFunnel", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), as.integer(intercrossingGenerations), PACKAGE="mpMap2"))
}
eightParentPedigreeSingleFunnelPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
{
  nonNegativeIntegerArgument(initialPopulationSize)
  nonNegativeIntegerArgument(selfingGenerations)
  nonNegativeIntegerArgument(nSeeds)
  nonNegativeIntegerArgument(intercrossingGenerations)
  intercrossingGenerations <- as.integer(intercrossingGenerations)
  initialPopulationSize <- as.integer(initialPopulationSize)
  selfingGenerations <- as.integer(selfingGenerations)
  nSeeds <- as.integer(nSeeds)

  if(initialPopulationSize <= 2 && intercrossingGenerations > 0)
  {
    stop("Random mating is impossible with only two lines per generation")
    #....and more importantly it means that the sample command below gets screwed up, because we're calling sample(x) where length(x) == 1, which samples from 1:x
  }
  entries <- 8L + 4L + 3L*initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize
  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:8] <- father[1:8] <- 0L
  mother[9:12] <- c(1L, 3L, 5L, 7L)
  father[9:12] <- c(2L, 4L, 6L, 8L)
  mother[13:(12+2*initialPopulationSize)] <- rep(c(9L, 11L), times = initialPopulationSize)
  father[13:(12+2*initialPopulationSize)] <- rep(c(10L, 12L), times = initialPopulationSize)
  mother[(13 + 2*initialPopulationSize):(12 + 3*initialPopulationSize)] <- seq(13L, 12L+2L*initialPopulationSize, by = 2L)
  father[(13 + 2*initialPopulationSize):(12 + 3*initialPopulationSize)] <- seq(14L, 12L+2L*initialPopulationSize, by = 2L)

  currentIndex <- 1L + 8L + 4L + initialPopulationSize*2L
  if(intercrossingGenerations > 0)
  {
    lastGenerationStart <- currentIndex
    lastGenerationEnd <- currentIndex-1L+initialPopulationSize
    for(i in 1:intercrossingGenerations)
    {
      for(lineCounter in lastGenerationStart:lastGenerationEnd)
      {
        mother[lineCounter + initialPopulationSize] <- lineCounter
        father[lineCounter + initialPopulationSize] <- sample(setdiff(lastGenerationStart:lastGenerationEnd, lineCounter), 1L)
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
    for(lineCounter in currentIndex:(currentIndex+initialPopulationSize-1L))
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
      observed[(nextFree-initialPopulationSize+1):length(observed)] <- TRUE
    }
    #If there's no selfing but there is intercrossing then mark the last set of intercrossing lines as observed
    else
    {
      observed[lastGenerationStart:lastGenerationEnd] <- TRUE
    }
  }
  if(length(mother) != entries || length(father) != entries || length(observed) != entries || length(lineNames) != entries)
  {
	  stop("Internal error")
  }
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:8L, observed = observed, selfing = "infinite", warnImproperFunnels = TRUE))
}
