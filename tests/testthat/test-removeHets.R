context("removeHets")

test_that("Check that removeHets works in F2",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		withoutHets <- cross + removeHets()
		expect_equal(sum(is.na(withoutHets@geneticData[[1]]@finals)), sum(cross@geneticData[[1]]@finals == 3))
		expect_identical(dimnames(cross@geneticData[[1]]@finals), dimnames(withoutHets@geneticData[[1]]@finals))
		expect_identical(dimnames(cross@geneticData[[1]]@founders), dimnames(withoutHets@geneticData[[1]]@founders))
		expect_identical(names(cross@geneticData[[1]]@hetData), names(withoutHets@geneticData[[1]]@hetData))
	})
test_that("Check that removeHets correctly handles the biparental case where the encodings don't start at 1",
	{
		pedigree <- rilPedigree(100, 5)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()

		#Alter everything by adding 10
		cross@geneticData[[1]]@founders <- cross@geneticData[[1]]@founders + 10L
		cross@geneticData[[1]]@finals <- cross@geneticData[[1]]@finals + 10L
		newHetData <- lapply(cross@geneticData[[1]]@hetData, function(x) x+10L)
		names(newHetData) <- names(cross@geneticData[[1]]@hetData)

		cross@geneticData[[1]]@hetData <- new("hetData", newHetData)
		withoutHets <- cross + removeHets()
		expect_identical(cross, withoutHets)
	})
test_that("Check that removeHets correctly handles the 4-parent case where the encodings don't start at 1",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
		expect_identical(cross, cross + removeHets())

		#Alter everything by adding 10
		cross@geneticData[[1]]@founders <- cross@geneticData[[1]]@founders + 10L
		cross@geneticData[[1]]@finals <- cross@geneticData[[1]]@finals + 10L
		newHetData <- lapply(cross@geneticData[[1]]@hetData, function(x) x+10L)
		names(newHetData) <- names(cross@geneticData[[1]]@hetData)

		cross@geneticData[[1]]@hetData <- new("hetData", newHetData)
		withoutHets <- cross + removeHets()
		expect_identical(cross, withoutHets)
	})
