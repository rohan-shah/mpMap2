context("sixteenParentPedigreeRandomFunnels")

test_that("Check that the C and R versions are identical",
	{
		parameters <- expand.grid(selfingGenerations = 0:2, nSeeds = 1:4, intercrossingGenerations = 0:3)
		seed <- 1
		apply(parameters, 1, function(x)
			{
				set.seed(seed)
				pedigreeR <- mpMap2:::sixteenParentPedigreeRandomFunnelsPrototype(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				set.seed(seed)
				pedigreeC <- sixteenParentPedigreeRandomFunnels(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				expect_identical(pedigreeR, pedigreeC)
				seed <<- seed + 1
			})
	})
