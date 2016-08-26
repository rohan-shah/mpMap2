context("Test subsetting of mpcrossMapped objects")
test_that("Subsetting of mapped objects by lines discards rf data",
{
	map <- sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")
	mapped <- new("mpcrossMapped", grouped, map = estimateMap(grouped), rf = rf@rf)
	expect_warning(subset(mapped, lines = rownames(finals(cross))[1:100]), "Subset function is discarding the recombination fraction data")

	subsettedMarkers <- subset(mapped, markers = names(mapped@map[[1]]))
	expect_is(subsettedMarkers, "mpcrossLG")
})
