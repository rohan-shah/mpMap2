context("Test estimateMap")
test_that("Check that argument maxMarkers works correctly", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- sim.map(len = rep(100, 1), n.mar = rep(401, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		estimated.map1 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 60)
		estimated.map2 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 40)
		estimated.map3 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 80)

		expect_equal(estimated.map1, estimated.map2, tolerance = 0.01)
		expect_equal(estimated.map1, estimated.map3, tolerance = 0.01)

		expect_identical(names(estimated.map1[[1]]), markers(cross))
		expect_identical(names(estimated.map2[[1]]), markers(cross))
		expect_identical(names(estimated.map3[[1]]), markers(cross))
	})

