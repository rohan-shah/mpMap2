test_that("Checking that addExtraMarkerFromRawCall works for an eight-parent pedigree", 
{
	set.seed(1)
	populationSize <- 1000

	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = populationSize, selfingGenerations = 5, nSeeds = 1L, intercrossingGenerations = 1L)
	originalMap <- qtl::sim.map(len = 100, n.mar = 201, anchor.tel = TRUE, include.x = FALSE, sex.sp = FALSE, eq.spacing = TRUE)
	cross <- simulateMPCross(map = originalMap, pedigree = pedigree, mapFunction = haldane, seed = 1)
        crossSNP <- cross + multiparentSNP(keepHets = FALSE)
	crossSNPMapped <- mpcrossMapped(crossSNP, map = originalMap)

	crossSNPImputed <- imputeFounders(crossSNPMapped, errorProb = 0.1)

	marker <- "D1M101"
	underlyingAlleles <- imputationData(crossSNPImputed)[, marker]

	rawData <- matrix(0.0, nrow = populationSize, ncol = 2)

	indices <- which(underlyingAlleles <= 4)
	rawData[indices, 1] <- rnorm(n = length(indices), mean = 0.25, sd = 1)
	rawData[indices, 2] <- rnorm(n = length(indices), mean = 0.5, sd = 1)

	indices <- which(underlyingAlleles >= 5)
	rawData[indices, 1] <- rnorm(n = length(indices), mean = 0.75, sd = 1)
	rawData[indices, 2] <- rnorm(n = length(indices), mean = 0.5, sd = 1)

	#Imputations without extra points
	rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = FALSE)
	newMarkerPosition <- originalMap[[1]][names(which.max(rawCall@data))]
	expect_equivalent(originalMap[[1]][marker], newMarkerPosition, tolerance = 2, scale = 1)
	expect_error(rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = TRUE))

	#Imputations with extra points (grid)
	crossSNPImputed <- imputeFounders(crossSNPMapped, errorProb = 0.1, extraPositions = generateGridPositions(0.25))
	rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = FALSE)
	newMarkerPosition <- imputationMap(crossSNPImputed)[[1]][names(which.max(rawCall@data))]
	expect_equivalent(originalMap[[1]][marker], newMarkerPosition, tolerance = 2, scale = 1)

	rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = TRUE)
	newMarkerPosition <- imputationMap(crossSNPImputed)[[1]][names(which.max(rawCall@data))]
	expect_equivalent(originalMap[[1]][marker], newMarkerPosition, tolerance = 2, scale = 1)

	#Imputations with extra points (mid-points)
	crossSNPImputed <- imputeFounders(crossSNPMapped, errorProb = 0.1, extraPositions = generateIntervalMidPoints)
	rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = FALSE)
	newMarkerPosition <- imputationMap(crossSNPImputed)[[1]][names(which.max(rawCall@data))]
	expect_equivalent(originalMap[[1]][marker], newMarkerPosition, tolerance = 2, scale = 1)

	rawCall <- addExtraMarkerFromRawCall(mpcrossMapped = crossSNPImputed, newMarker = rawData, useOnlyExtraImputationPoints = TRUE)
	newMarkerPosition <- imputationMap(crossSNPImputed)[[1]][names(which.max(rawCall@data))]
	expect_equivalent(originalMap[[1]][marker], newMarkerPosition, tolerance = 2, scale = 1)
})
