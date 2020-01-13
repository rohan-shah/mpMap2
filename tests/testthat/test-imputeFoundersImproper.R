context("Founder imputation isn't implemented for improper funnels (funnels with repeated founders)")
test_that("Eight parent design",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			expect_error(result <- imputeFounders(mapped, errorProb = 0.1))
			expect_error(result <- computeGenotypeProbabilities(mapped, errorProb = 0.1))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 201, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 201, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- eightParentPedigreeImproperFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1)
		pedigree1@selfing <- "infinite"
		pedigree2 <- eightParentPedigreeImproperFunnels(initialPopulationSize = 500, selfingGenerations = 10, nSeeds = 1)
		pedigree2@selfing <- "infinite"
		pedigree3 <- eightParentPedigreeImproperFunnels(initialPopulationSize = 500, selfingGenerations = 5, nSeeds = 1)
		pedigree3@selfing <- "finite"
		pedigree4 <- eightParentPedigreeImproperFunnels(initialPopulationSize = 500, selfingGenerations = 5, nSeeds = 1)
		pedigree4@selfing <- "finite"
		
		pedigrees <- list(pedigree1, pedigree2, pedigree3, pedigree4)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}
	})
