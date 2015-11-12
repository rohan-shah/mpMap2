context("biparentalDominant")

test_that("Check that biparentalDominant works for an F2",
	{
		pedigree <- f2Pedigree(100)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_true(all(cross@geneticData[[1]]@finals <= 2))
		expect_true(all(unlist(lapply(cross@geneticData[[1]]@hetData, function(x) length(unique(x[,3])) == 2))))
	})

test_that("Check that biparentalDominant works for a RIL",
	{
		pedigree <- rilPedigree(100, 2)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_true(all(cross@geneticData[[1]]@finals <= 2))
		expect_true(all(unlist(lapply(cross@geneticData[[1]]@hetData, function(x) length(unique(x[,3])) == 2))))
	})
test_that("Check that biparentalDominant can't be applied twice",
	{
		pedigree <- f2Pedigree(100)
		map <- sim.map(len = 100, n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		expect_that(cross+biparentalDominant(), throws_error())
	})
