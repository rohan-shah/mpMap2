context("Test subsetting of mpcrossRF objects")
test_that("Subsetting of mpcrossRF object by lines discards rf data",
{
	map <- qtl::sim.map(len = rep(100, 2), n.mar = rep(11, 2), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	expect_warning(subsetted <- subset(rf, lines = lineNames(cross)[1:10]), "Discarding rf")
	expect_equivalent(class(subsetted), "mpcross")
})
