context("hetData validation")

#Set up a basic mpcross object
twoWayGeneticDataSetUp <- function()
{
	env <- sys.frame(sys.parent(1))
	pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
	map <- sim.map(len = rep(100, 1), n.mar = 11, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
	cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
	return(cross@geneticData[[1]])
}
test_that("geneticData rejects floating point values", 
	{
		geneticData <- twoWayGeneticDataSetUp()
		expect_identical(validObject(geneticData), TRUE)

		#floating point value in founders given an error
		copied <- geneticData
		copied@founders[1,1] <- 0.5
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#floating point value in finals gives an exception
		copied <- geneticData
		copied@finals[1,1] <- 0.5
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Setting these arrays to numeric gives no exception
		copied <- geneticData
		storage.mode(copied@founders) <- "numeric"
		storage.mode(copied@finals) <- "numeric"
		for(i in 1:length(copied@hetData)) storage.mode(copied@hetData[[i]]) <- "numeric"
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
test_that("geneticData rejects invalid values", 
	{
		geneticData <- twoWayGeneticDataSetUp()
		expect_identical(validObject(geneticData), TRUE)
		
		#Invalid final genotype throws error
		copied <- geneticData
		copied@finals[1,1] <- 100
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Invalid founder genotype
		copied <- geneticData
		copied@founders[1,1] <- 100
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#This causes the hetData to mismatch the founders, so it throws an error
		copied <- geneticData
		copied@founders[,1] <- 1
		copied@finals[,1] <- 1
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#But if we now fix the hetData, there's no error
		copied@hetData[[1]] <- rbind(c(1,1,1))
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
test_that("Missing data correctly validated in two-way geneticData", 
	{
		geneticData <- twoWayGeneticDataSetUp()
		expect_identical(validObject(geneticData), TRUE)

		#An NA founder throws an error
		copied <- geneticData
		copied@founders[1,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#All NA founders throw an error
		copied <- geneticData
		copied@founders[,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#All NA founders and finals throw an error
		copied <- geneticData
		copied@founders[,1] <- NA
		copied@finals[,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Deleting the hetData throws an error
		copied <- geneticData
		copied@hetData[[1]] <- matrix(0L, 0, 3)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Deleting the hetData and setting all the finals gives no error (There are no finals, so every final can be *found* in the hetData)
		copied <- geneticData
		copied@hetData[[1]] <- matrix(0L, 0, 3)
		copied@finals[,1] <- NA
		expect_identical(validObject(copied, complete=TRUE), TRUE)

		#Deleting the hetData and the founders gives an error
		copied <- geneticData
		copied@hetData[[1]] <- matrix(0L, 0, 3)
		copied@founders[,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#If all founders and finals are NA, and the hetData is a 0x3 matrix, then no error
		copied <- geneticData
		copied@founders[,1] <- NA
		copied@finals[,1] <- NA
		copied@hetData[[1]] <- matrix(0L, 0, 3)
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})