context("eightParentPedigreeSingleFunnel")

test_that("Check that the C and R versions are identical",
	{
		parameters <- expand.grid(selfingGenerations = 0:3, intercrossingGenerations = 0:3, nSeeds = 1:3)
		seed <- 1
		apply(parameters, 1, function(x)
			{
				set.seed(seed)
				pedigreeR <- mpMap2:::eightParentPedigreeSingleFunnelPrototype(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				set.seed(seed)
				pedigreeC <- eightParentPedigreeSingleFunnel(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				expect_identical(pedigreeR, pedigreeC)
				seed <<- seed + 1
			})
	})
