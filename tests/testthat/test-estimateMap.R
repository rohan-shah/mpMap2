context("Test estimateMap")
test_that("Check that argument maxMarkers works correctly", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = rep(100, 1), n.mar = rep(101, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		estimated.map1 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 60)
		estimated.map2 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 40)
		estimated.map3 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 80)
		estimated.map4 <- estimateMap(grouped, maxOffset = 4, maxMarkers = 1000)

		expect_equal(estimated.map1, estimated.map2, tolerance = 0.02)
		expect_equal(estimated.map1, estimated.map3, tolerance = 0.02)
		expect_equal(estimated.map1, estimated.map4, tolerance = 0.02)

		expect_identical(names(estimated.map1[[1]]), markers(cross))
		expect_identical(names(estimated.map2[[1]]), markers(cross))
		expect_identical(names(estimated.map3[[1]]), markers(cross))
		expect_identical(names(estimated.map4[[1]]), markers(cross))
	})

