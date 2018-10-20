context("Test subsetting of mpcrossLG objects")
test_that("Subsetting of mpcrossLG objects by lines discards rf data, and warns about retaining linkage group data",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	grouped <- formGroups(rf, groups = 2, clusterBy = "theta", method = "average")

	warnings <- testthat::capture_warnings(subsetted <- subset(grouped, lines = rownames(finals(cross))[1:100]))
	expect_match(warnings, "Retaining linkage group data", all = FALSE)
	expect_match(warnings, "Discarding rf data", all = FALSE)

	expect_warning(subsetted <- subset(grouped, lines = rownames(finals(cross))[1:100]), "Discarding rf data")
})
