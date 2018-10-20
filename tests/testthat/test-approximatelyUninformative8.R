context("uninformative 8-parent markers")
test_that("Check that the eight-parent uninformative marker combination gives an RF estimate of NA in the right cases",
	{
		testInfiniteSelfing <- function(pedigree)
		{
			pedigree@selfing <- "infinite"
			map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel=TRUE, include.x=FALSE, sex.sp=FALSE, eq.spacing=TRUE)
			cross <- simulateMPCross(pedigree = pedigree, map = map, mapFunction = haldane)

			firstColumnFunction <- function(x)
			{
				if(x %in% c(1, 4, 5, 7, 8)) return(1)
				if(x %in% c(2, 3, 6)) return(0)
				return(NA)
			}
			secondColumnFunction <- function(x)
			{
				if(x %in% c(3, 7, 8)) return(1)
				if(x %in% c(1, 2, 4, 5, 6)) return(0)
				return(NA)
			}
			cross@geneticData[[1]]@founders[,1] <- sapply(cross@geneticData[[1]]@founders[,1], firstColumnFunction)
			cross@geneticData[[1]]@finals[,1] <- sapply(cross@geneticData[[1]]@finals[,1], firstColumnFunction)

			cross@geneticData[[1]]@founders[,2] <- sapply(cross@geneticData[[1]]@founders[,2], secondColumnFunction)
			cross@geneticData[[1]]@finals[,2] <- sapply(cross@geneticData[[1]]@finals[,2], secondColumnFunction)
			cross@geneticData[[1]]@hetData[[1]] <- cross@geneticData[[1]]@hetData[[2]] <- rbind(c(0,0,0), c(1,1,1))

			validObject(cross)

			return(estimateRF(cross))
		}
		#Infinite selfing, with or without intercrossing, is more uninformative for the single funnel design than for the random funnel design
		for(intercrossingGenerations in 0:1)
		{
			pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			rf <- testInfiniteSelfing(pedigree)
			expect_true(is.na(rf@rf@theta[1,2]))
		}
		intercrossingGenerations <- 0
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
		rf <- testInfiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))

		intercrossingGenerations <- 1
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
		rf <- testInfiniteSelfing(pedigree)
		expect_true(is.na(rf@rf@theta[1,2]))
	})
