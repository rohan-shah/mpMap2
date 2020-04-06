#' @title Generate an eight-parent pedigree, using random funnels
#'
#' @description
#' Generate a eight-parent pedigree starting from inbred founders, using a random funnel. 
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
#' pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 10, 
#' 	selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 10)
#' #Generate map
#' map <- qtl::sim.map()
#' #Simulate data
#' cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane)
#' #Get out a list of funnels, which are rows of this matrix. For this pedigree, every 
#' #    funnel is a random ordering of 1:8. 
#' getAllFunnels(cross)
#' #convert the pedigree to a graph
#' pedigreeAsGraph <- pedigreeToGraph(pedigree)
#' #Plot it
#' \donttest{plot(pedigreeAsGraph)}
#' #Write it to a file in DOT format
# \donttest{write.graph(graph = pedigreeAsGraph@@graph, format = "dot", file = "./pedigree.dot")}

# This is written in C because otherwise it's just too damn slow (especially for generating the huge populations that we want to use to get numerically accurate results for unit testing)
eightParentPedigreeRandomFunnels <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  return(.Call("eightParentPedigreeRandomFunnels", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), as.integer(intercrossingGenerations), PACKAGE="mpMap2"))
}
eightParentPedigreeRandomFunnelsPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  #All pairs
  pairs <- combn(1:8, 2)
  quads <- vector()
  #Generate all the quads. We do this by taking every pair, and combining it with all the other pairs which contain different founders, and which are strictly after the smallest founder of the first pair. 
  for (i in 1:(ncol(pairs)-3))
  {
    quads <- cbind(quads, rbind(matrix(rep(pairs[,i], choose(6-min(pairs[,i])+1, 2)), nrow=2, ncol=choose(6-min(pairs[,i])+1, 2)), combn(setdiff((min(pairs[,i])+1):8, pairs[,i]), 2)))
  }

  #Now we just take pairs of quads which are disjoint
  funnels <- vector()
  for (i in 1:(ncol(quads)-1))
  {
    addons <- which(lapply(apply(quads[,i:ncol(quads)], 2, function(x) return(intersect(quads[,i], x))), length)==0)
    funnels <- cbind(funnels, rbind(matrix(rep(quads[,i], length(addons)), nrow=4, ncol=length(addons)), quads[,(i:ncol(quads))[addons]]))
  }
  funnelNumbers <- sample(1:ncol(funnels), initialPopulationSize, replace=TRUE)

  entries <- 8L + 7L*initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize
  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:8] <- father[1:8] <- 0L

  #We have to generate new funnels for every line
  #Create four two-way crosses for funnels
  for (i in 1:initialPopulationSize)
  {
    mother[8L + (i-1)*4+1:4] <- funnels[c(1,3,5,7), funnelNumbers[i]]
    father[8L + (i-1)*4+1:4] <- funnels[c(2,4,6,8), funnelNumbers[i]]
  }
  #Create two four-way crosses for each funnel
  mother[8L + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 7L + 2L*(1:(2*initialPopulationSize))
  father[8L + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 8L + 2L*(1:(2*initialPopulationSize))
  #Create a line which mixes all eight founders, for each funnel
  mother[8L + initialPopulationSize*6 + 1:initialPopulationSize] <- 7L + initialPopulationSize*4L + 2L*(1:initialPopulationSize)
  father[8L + initialPopulationSize*6 + 1:initialPopulationSize] <- 8L + initialPopulationSize*4L + 2L*(1:initialPopulationSize)

  currentIndex <- 1L + 8L + initialPopulationSize*5L + initialPopulationSize
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
      observed[(nextFree-initialPopulationSize):length(observed)] <- TRUE
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
