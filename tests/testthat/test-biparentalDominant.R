context("biparentalDominant")

test_that("Check that biparentalDominant works for an F2",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_true(all(cross@geneticData[[1]]@finals <= 2))
		expect_true(all(unlist(lapply(cross@geneticData[[1]]@hetData, function(x) length(unique(x[,3])) == 2))))
	})

test_that("Check that biparentalDominant works for a RIL",
	{
		pedigree <- rilPedigree(100, 2)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_true(all(cross@geneticData[[1]]@finals <= 2))
		expect_true(all(unlist(lapply(cross@geneticData[[1]]@hetData, function(x) length(unique(x[,3])) == 2))))
	})
test_that("Check that biparentalDominant doesn't work on a 4-way design",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		expect_that(cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biParentalDominant(), throws_error())
	})
test_that("Check that biparentalDominant can't be applied twice",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_that(cross+biparentalDominant(), throws_error())
	})
test_that("Check that biparentalDominant can't be applied to multiple datasets at once",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross1 <- subset(cross, markers = 1:5)
		cross2 <- subset(cross, markers = 6:10)

		rownames(cross2@geneticData[[1]]@finals) <- paste0(rownames(cross2@geneticData[[1]]@finals), ",2")
		rownames(cross2@geneticData[[1]]@founders) <- paste0(rownames(cross2@geneticData[[1]]@founders), ",2")
		cross2@geneticData[[1]]@pedigree@lineNames <- paste0(cross2@geneticData[[1]]@pedigree@lineNames, ",2")
		cross <- cross1 + cross2
		expect_that(cross + biparentalDominant(), throws_error())
	})
