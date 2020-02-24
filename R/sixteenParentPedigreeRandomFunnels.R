#' @title Generate a sixteen-parent pedigree
#'
#' @description
#' Generate a sixteen-parent pedigree starting from inbred founders, using a random funnel
#'
#' @seealso \code{\link{eightParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeSingleFunnel}}, \code{\link{fourParentPedigreeRandomFunnels}}, \code{\link{twoParentPedigree}}
#'
#' @param initialPopulationSize The number of F1 lines generated
#' @param selfingGenerations The number of selfing generations at the end of the pedigree
#' @param nSeeds The number of progeny taken from each intercrossing line, or from each F1 if no intercrossing is specified. These lines are then selfed according to selfingGenerations
#' @param intercrossingGenerations The number of generations of random mating performed from the F1 generation. Population size is maintained at that specified by initialPopulationSize
#' @return An object of class \code{detailedPedigree} representing the experimental design, suitable for simulation using simulateMPCross. 
#' @export
# This is written in C because otherwise it's just too damn slow (especially for generating the huge populations that we want to use to get numerically accurate results for unit testing)
sixteenParentPedigreeRandomFunnels <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
{
	nonNegativeIntegerArgument(initialPopulationSize)
	nonNegativeIntegerArgument(selfingGenerations)
	nonNegativeIntegerArgument(nSeeds)
	nonNegativeIntegerArgument(intercrossingGenerations)
	intercrossingGenerations <- as.integer(intercrossingGenerations)
	initialPopulationSize <- as.integer(initialPopulationSize)
	selfingGenerations <- as.integer(selfingGenerations)
	return(.Call("sixteenParentPedigreeRandomFunnels", as.integer(initialPopulationSize), as.integer(selfingGenerations), as.integer(nSeeds), as.integer(intercrossingGenerations), PACKAGE="mpMap2"))
}
sixteenParentPedigreeRandomFunnelsPrototype <- function(initialPopulationSize, selfingGenerations, nSeeds = 1L, intercrossingGenerations)
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
  pairs <- combn(1:16, 2)
  #Arrange all pairs so that the minimum is in the first row, and the maximum in the second row
  pairs <- apply(pairs, 2, function(x) c(min(x), max(x)))

  funnels <- matrix(1L, nrow = initialPopulationSize, ncol = 16)
  for (i in 1:initialPopulationSize)
  {
    funnels[i, ] <- sample(1:16, 16, replace = FALSE)
  }

  entries <- 16L + ncol(pairs) + 7L*initialPopulationSize + intercrossingGenerations*initialPopulationSize + nSeeds*selfingGenerations*initialPopulationSize
  mother <- father <- rep(NA, entries)
  observed <- rep(FALSE, entries)
  lineNames <- paste0("L", 1:entries)
  mother[1:16] <- father[1:16] <- 0L

  #Put in all pairs
  for(i in 1:ncol(pairs))
  {
    mother[16L+i] <- pairs[1, i]
    father[16L+i] <- pairs[2, i]
  }

  #We have to generate new funnels for every line
  for (i in 1:initialPopulationSize)
  {
    #Find the id of the pair among the generated pairs
    pairsThisFunnel <- rbind(funnels[i,c(1,3,5,7,9,11,13,15)], funnels[i,c(2,4,6,8,10,12,14,16)])
    pairIDsThisFunnel <- apply(pairsThisFunnel, 2, function(x) which(pairs[1,] == min(x) & pairs[2,] ==  max(x)))
    mother[16L + ncol(pairs) + (i-1)*4+1:4] <- pairIDsThisFunnel[c(1,3,5,7)]+16L
    father[16L + ncol(pairs) + (i-1)*4+1:4] <- pairIDsThisFunnel[c(2,4,6,8)]+16L
  }
  #Create two eight-way crosses for each funnel
  mother[16L + ncol(pairs) + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 15L + ncol(pairs) + 2L*(1:(2*initialPopulationSize))
  father[16L + ncol(pairs) + initialPopulationSize*4 + 1:(2*initialPopulationSize)] <- 16L + ncol(pairs) + 2L*(1:(2*initialPopulationSize))
  #Create a line which mixes all sixteen founders, for each funnel
  mother[16L + ncol(pairs) + initialPopulationSize*6 + 1:initialPopulationSize] <- 15L + ncol(pairs) + initialPopulationSize*4L + 2L*(1:initialPopulationSize)
  father[16L + ncol(pairs) + initialPopulationSize*6 + 1:initialPopulationSize] <- 16L + ncol(pairs) + initialPopulationSize*4L + 2L*(1:initialPopulationSize)

  currentIndex <- 1L + 16L + ncol(pairs) + initialPopulationSize*6L
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
  return(new("detailedPedigree", lineNames = lineNames, mother = mother, father = father, initial = 1L:16L, observed = observed, selfing = "infinite", warnImproperFunnels = TRUE))
}
