context("eightParentPedigreeRandomFunnels")

test_that("Check that the C and R versions are identical",
	{
		parameters <- expand.grid(selfingGenerations = 0:3, intercrossingGenerations = 0:3, nSeeds = 1:3)
		seed <- 1
		apply(parameters, 1, function(x)
			{
				set.seed(seed)
				pedigreeR <- mpMap2:::eightParentPedigreeRandomFunnelsPrototype(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				set.seed(seed)
				pedigreeC <- eightParentPedigreeRandomFunnels(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				expect_identical(pedigreeR, pedigreeC)
				seed <<- seed + 1
			})
	})

test_that("Correct number of lines are generated with zero generations of selfing",
{
	pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 10, selfingGenerations = 0, intercrossingGenerations = 0)
	map <- qtl::sim.map(len = 1, n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	cross <- simulateMPCross(pedigree = pedigree, map = map, mapFunction = haldane, seed =1)
	expect_equal(nLines(cross), 10)
})
