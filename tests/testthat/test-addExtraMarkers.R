test_that("Test validation",
{
	map <- qtl::sim.map(len = 20, n.mar = 21, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(200)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	cross <- new("mpcrossMapped", cross, map = map)

	subset1 <- subset(cross, markers = c(1:10, 16:21), keepMap = TRUE)
	subset2 <- subset(cross, markers = 11:15)
	#No imputation data
	expect_error(addExtraMarkers(subset1, newMarkers = subset2))

	#No RF data
	subset1 <- imputeFounders(subset1, extraPositions = generateGridPositions(1))
	expect_error(addExtraMarkers(subset1, newMarkers = subset2))
	subset1 <- estimateRF(subset1)
	
	#Now it works
	statistics <- addExtraMarkers(subset1, newMarkers = subset2, onlyStatistics = TRUE, verbose = FALSE)
	expect_true(is.numeric(statistics@data))

	#Check that the resulting objcet is correct. 
	capture.output(results <- addExtraMarkers(subset1, newMarkers = subset2, verbose = FALSE))
	permutation <- match(markers(results$object), markers(cross))
	expect_gt(abs(cor(permutation, 1:nMarkers(cross))), 0.91)
})
test_that("Test that markers go to the right place",
{
	map <- qtl::sim.map(len = 200, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(200)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	wholeRF <- estimateRF(cross)
	for(relevantMarker in c(1, 10, 50, 90, 101))
	{
		subset1 <- subset(cross, markers = (1:101)[-relevantMarker], keepMap = TRUE)
		subset2 <- subset(cross, markers = relevantMarker)
		subset1RF <- subset(wholeRF, markers = markers(subset1))

		grouped <- formGroups(subset1RF, groups = 1, clusterBy = "theta", method = "average")
		estimatedMap <- estimateMap(grouped, maxOffset = 15)
		mapped <- new("mpcrossMapped", grouped, map = estimatedMap, rf = grouped@rf)

		imputed <- imputeFounders(mapped, errorProb = 0.1, extraPositions = generateGridPositions(1))
		capture.output(added <- addExtraMarkers(imputed, newMarkers = subset2, reorder = FALSE, maxOffset = 15, reorderRadius = 40))
		permutation <- sapply(markers(added$object), function(x) match(x, markers(cross)))
		expect_gt(cor(1:101, permutation), 0.99)
		reestimatedMap <- estimateMap(added$object, maxOffset = 15)
		expect_gt(cor(reestimatedMap[[1]], added$object@map[[1]]), 0.99)

		#Check that adding the RF data has worked
		reorderedRF <- subset(wholeRF@rf, markers = markers(added$object))
		expect_identical(reorderedRF, added$object@rf)
	}
})
