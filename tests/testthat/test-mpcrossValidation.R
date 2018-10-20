context("mpcross validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
map <- qtl::sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)

test_that("Simulated cross passes validation",
	{
		expect_identical(validObject(cross, complete=TRUE), TRUE)
	})
test_that("At least one set of geneticData is required",
	{
		copied <- cross
		copied@geneticData <- new("geneticDataList", list())
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("All elements of geneticData must have correct class",
	{
		copied <- cross
		expect_that(copied@geneticData <- new("geneticDataList", c(copied@geneticData, 1)), throws_error())

		copied <- cross
		expect_that(copied@geneticData <- new("geneticDataList", c(copied@geneticData, cross)), throws_error())
	})
test_that("If there are multiple genetic data sets, founder line names can be repeated",
	{
		copied <- cross
		copied2 <- cross

		rownames(copied2@geneticData[[1]]@finals) <- paste0(rownames(copied2@geneticData[[1]]@finals), ",2")
		pedigreeSubset <- copied2@geneticData[[1]]@pedigree@lineNames %in% rownames(copied@geneticData[[1]]@finals)
		copied2@geneticData[[1]]@pedigree@lineNames[pedigreeSubset] <- paste0(copied2@geneticData[[1]]@pedigree@lineNames[pedigreeSubset], ",2")
		combined <- copied + copied2
		expect_identical(validObject(combined, complete=TRUE), TRUE)
	})
test_that("If there are multiple genetic data sets, line names must still be unique",
	{
		copied <- cross
		copied2 <- cross

		rownames(copied2@geneticData[[1]]@finals) <- paste0(rownames(copied2@geneticData[[1]]@finals), ",2")
		rownames(copied2@geneticData[[1]]@founders) <- paste0(rownames(copied2@geneticData[[1]]@founders), ",2")
		copied2@geneticData[[1]]@pedigree@lineNames <- paste0(copied2@geneticData[[1]]@pedigree@lineNames, ",2")
		expect_identical(validObject(copied2, complete=TRUE), TRUE)

		combined <- copied + copied2
		expect_identical(validObject(combined, complete=TRUE), TRUE)
		
		copied2@geneticData[[1]]@pedigree@lineNames[copied2@geneticData[[1]]@pedigree@lineNames == rownames(copied2@geneticData[[1]]@finals)[1]] <- rownames(copied@geneticData[[1]]@finals)[1]
		rownames(copied2@geneticData[[1]]@finals)[1] <- rownames(copied@geneticData[[1]]@finals)[1]
		expect_identical(validObject(copied2, complete=TRUE), TRUE)
		expect_that(combined <- copied + copied2, throws_error())
	})
#Errors in checkCompatibleGeneticData
test_that("All geneticData entries must have the same markers",
	{
		#Having two sets of genetic data is fine
		copied <- cross
		copied@geneticData <- new("geneticDataList", c(copied@geneticData, copied@geneticData))
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
		copied@geneticData <- new("geneticDataList", c(copied@geneticData, copied@geneticData))
		copied@geneticData[[1]]@founders <- copied@geneticData[[1]]@founders[,-1]
		copied@geneticData[[1]]@finals <- copied@geneticData[[1]]@finals[,-1]
		copied@geneticData[[1]]@hetData[[1]] <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Having different marker names is an error
		copied <- cross
		copied@geneticData <- new("geneticDataList", c(copied@geneticData, copied@geneticData))
		colnames(copied@geneticData[[1]]@founders)[1] <- colnames(copied@geneticData[[1]]@finals)[1] <- names(copied@geneticData[[1]]@hetData)[1] <- "newMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Object with no final lines is valid",
	{
		pedigree <- f2Pedigree(10)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross@geneticData[[1]]@finals <-cross@geneticData[[1]]@finals[0,,drop=F]
		expect_identical(validObject(cross,complete=TRUE), TRUE)
	})
rm(pedigree, map, cross)
