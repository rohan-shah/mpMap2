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
				if(x %in% c(1, 3, 5, 7)) return(0)
				if(x %in% c(2, 4, 6, 8)) return(1)
				return(NA)
			}
			secondColumnFunction <- function(x)
			{
				if(x %in% c(3, 4, 7, 8)) return(0)
				if(x %in% c(1, 2, 5, 6)) return(1)
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
		testFiniteSelfing <- function(pedigree)
		{
			pedigree@selfing <- "finite"
			map <- qtl::sim.map(len = 10, n.mar = 2, anchor.tel=TRUE, include.x=FALSE, sex.sp=FALSE, eq.spacing=TRUE)
			cross <- simulateMPCross(pedigree = pedigree, map = map, mapFunction = haldane)
			firstHetData <- cross@geneticData[[1]]@hetData[[1]]
			secondHetData <- cross@geneticData[[1]]@hetData[[2]]

			firstColumnFunction <- function(x)
			{
				hetDataIndex <- match(x, firstHetData[,3])
				founder1 <- firstHetData[hetDataIndex, 1]
				founder2 <- firstHetData[hetDataIndex, 2]
				allele1 <- c(0, 1, 0, 1, 0, 1, 0, 1)[founder1]
				allele2 <- c(0, 1, 0, 1, 0, 1, 0, 1)[founder2]
				if(allele1 == allele2)
				{
					return(allele1)
				}
				return(2)
			}
			secondColumnFunction <- function(x)
			{
				hetDataIndex <- match(x, secondHetData[,3])
				founder1 <- secondHetData[hetDataIndex, 1]
				founder2 <- secondHetData[hetDataIndex, 2]
				allele1 <- c(1, 1, 0, 0, 1, 1, 0, 0)[founder1]
				allele2 <- c(1, 1, 0, 0, 1, 1, 0, 0)[founder2]
				if(allele1 == allele2)
				{
					return(allele1)
				}
				return(2)
			}
			cross@geneticData[[1]]@founders[,1] <- sapply(cross@geneticData[[1]]@founders[,1], firstColumnFunction)
			cross@geneticData[[1]]@finals[,1] <- sapply(cross@geneticData[[1]]@finals[,1], firstColumnFunction)

			cross@geneticData[[1]]@founders[,2] <- sapply(cross@geneticData[[1]]@founders[,2], secondColumnFunction)
			cross@geneticData[[1]]@finals[,2] <- sapply(cross@geneticData[[1]]@finals[,2], secondColumnFunction)
			cross@geneticData[[1]]@hetData[[1]] <- cross@geneticData[[1]]@hetData[[2]] <- rbind(c(0,0,0), c(1,1,1), c(0, 1, 2), c(1, 0, 2))

			validObject(cross)

			return(estimateRF(cross))
		}
		#Infinite selfing, with or without intercrossing, is uninformative
		for(intercrossingGenerations in 0:1)
		{
			pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			rf <- testInfiniteSelfing(pedigree)
			expect_true(is.na(rf@rf@theta[1,2]))
			pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = intercrossingGenerations, nSeeds = 1)
			rf <- testInfiniteSelfing(pedigree)
			expect_true(is.na(rf@rf@theta[1,2]))
		}
		#Now finite generations of selfing
		#1 generation of intercrossing and 0 generations of selfing is uninformative for the single funnel design
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(is.na(rf@rf@theta[1,2]))
		#1 generation of intercrossing and 0 generations of selfing is uninformative for the random funnel design
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 0, intercrossingGenerations = 1, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(is.na(rf@rf@theta[1,2]))

		#Zero generations of intercrossing and zero generations of selfing, the result is uninformative for the single funnel design
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(is.na(rf@rf@theta[1,2]))
		#Zero generations of intercrossing and zero generations of selfing, the result is informative for the random funnel design
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 0, intercrossingGenerations = 0, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))

		#Both selfing and intercrossing, the result is informative
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 1, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 1, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))
		
		#Selfing without intercrossing, the result is informative
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
		rf <- testFiniteSelfing(pedigree)
		expect_true(!is.na(rf@rf@theta[1,2]))
	})
