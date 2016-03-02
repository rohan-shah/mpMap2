context("fourParentPedigreeSingleFunnel")

test_that("Check that the C and R versions are identical",
	{
		parameters <- expand.grid(selfingGenerations = 0:3, intercrossingGenerations = 0:3, nSeeds = 1:3)
		seed <- 1
		apply(parameters, 1, function(x)
			{
				set.seed(seed)
				pedigreeR <- mpMap2:::fourParentPedigreeSingleFunnelPrototype(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				set.seed(seed)
				pedigreeC <- fourParentPedigreeSingleFunnel(selfingGenerations = x["selfingGenerations"], intercrossingGenerations = x["intercrossingGenerations"], initialPopulationSize = 20, nSeeds = x["nSeeds"])
				expect_equal(pedigreeR, pedigreeC)
				seed <<- seed + 1
			})
	})
test_that("Input nSeeds must be positive",
	{
			expect_that(fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, intercrossingGenerations = 0, nSeeds = 0), throws_error("positive integer"))
	})
