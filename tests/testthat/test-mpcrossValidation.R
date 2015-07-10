context("mpcross validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
map <- sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)

test_that("Simulated cross passes validation",
	{
		expect_identical(validObject(cross, complete=TRUE), TRUE)
	})
test_that("At least one set of geneticData is required",
	{
		copied <- cross
		copied@geneticData <- list()
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("All elements of geneticData must have correct class",
	{
		copied <- cross
		copied@geneticData <- c(copied@geneticData, 1)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- cross
		copied@geneticData <- c(copied@geneticData, cross)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
#Errors in checkCompatibleGeneticData
test_that("All geneticData entries must have the same markers",
	{
		#Having two sets of genetic data is fine
		copied <- cross
		copied@geneticData <- c(copied@geneticData, copied@geneticData)
		expect_identical(validObject(cross, complete=TRUE), TRUE)

		#Removing a marker is fine
		copied <- cross
		copied@geneticData[[1]]@founders <- copied@geneticData[[1]]@founders[,-1]
		copied@geneticData[[1]]@finals <- copied@geneticData[[1]]@finals[,-1]
		copied@geneticData[[1]]@hetData[[1]] <- NULL
		expect_identical(validObject(cross, complete=TRUE), TRUE)

		#Changing the name of a marker is fine
		copied <- cross
		colnames(copied@geneticData[[1]]@founders)[1] <- colnames(copied@geneticData[[1]]@finals)[1] <- names(copied@geneticData[[1]]@hetData)[1] <- "newMarker"
		expect_identical(validObject(cross, complete=TRUE), TRUE)

		#Having different numbers of markers is an error
		copied <- cross
		copied@geneticData <- c(copied@geneticData, copied@geneticData)
		copied@geneticData[[1]]@founders <- copied@geneticData[[1]]@founders[,-1]
		copied@geneticData[[1]]@finals <- copied@geneticData[[1]]@finals[,-1]
		copied@geneticData[[1]]@hetData[[1]] <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Having different marker names is an error
		copied <- cross
		copied@geneticData <- c(copied@geneticData, copied@geneticData)
		colnames(copied@geneticData[[1]]@founders)[1] <- colnames(copied@geneticData[[1]]@finals)[1] <- names(copied@geneticData[[1]]@hetData)[1] <- "newMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})