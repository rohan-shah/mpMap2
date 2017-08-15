context("Test addExtraMarkers")

test_that("Test validation",
{
	map <- sim.map(len = 20, n.mar = 21, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(200)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	cross <- new("mpcrossMapped", cross, map = map)

	subset1 <- subset(cross, markers = c(1:10, 16:21), keepMap = TRUE)
	subset2 <- subset(cross, markers = 11:15)
	#No imputation data
	expect_error(addExtraMarkers(subset1, newMarkers = subset2, attemptMpMap2Interactive = FALSE))

	#No RF data
	subset1 <- imputeFounders(subset1, extraPositions = generateGridPositions(1))
	expect_error(addExtraMarkers(subset1, newMarkers = subset2, attemptMpMap2Interactive = FALSE))
	subset1 <- estimateRF(subset1)
	
	#Now it works
	statistics <- addExtraMarkers(subset1, newMarkers = subset2, attemptMpMap2Interactive = FALSE, onlyStatistics = TRUE, verbose = FALSE)
	expect_true(is.numeric(statistics))

	#Check that the resulting objcet is correct. 
	results <- addExtraMarkers(subset1, newMarkers = subset2, attemptMpMap2Interactive = FALSE, verbose = FALSE)
	permutation <- match(markers(results$object), markers(cross))
	expect_gt(abs(cor(permutation, 1:nMarkers(cross))), 0.95)
})
