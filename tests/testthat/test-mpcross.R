context("Test mpcross function")
test_that("Test that mpcross works without inputting a value for hetData",
	{
		rilPedigree <- rilPedigree(100, selfingGenerations = 8)
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=rilPedigree, mapFunction = haldane)

		expect_error(cross1 <- mpMap2::mpcross(founders = founders(cross), finals = finals(cross), pedigree = rilPedigree))
		expect_warning(cross1 <- mpMap2::mpcross(founders = founders(cross), finals = finals(cross), pedigree = rilPedigree, fixCodingErrors = TRUE))
		expect_warning(cross2 <- mpMap2::mpcross(founders = founders(cross), finals = finals(cross), pedigree = rilPedigree, hetData = infiniteSelfing, fixCodingErrors=TRUE))
		expect_identical(cross1, cross2)
	})
test_that("Test that mpcross works with hetData = infiniteSelfing",
	{
		rilPedigree <- rilPedigree(100, selfingGenerations = 8)
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=rilPedigree, mapFunction = haldane)

		#No need for fixCodingErrors here
		expect_error(cross2 <- mpMap2::mpcross(founders = founders(cross), finals = finals(cross), pedigree = rilPedigree, hetData = hetsForSNPMarkers), NA)
	})
