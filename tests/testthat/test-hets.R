context("Test use of hets, when only *one* is informative for hets")

test_that("Test imputation",
{
	map <- qtl::sim.map(len = 5, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE, n.mar = 21)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
	pedigree@selfing <- "finite"
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 2)

	cross2 <- cross + multiparentSNP(keepHets = FALSE)
	cross2@geneticData[[1]]@finals[,11] <- cross@geneticData[[1]]@finals[,11]
	cross2@geneticData[[1]]@founders[,11] <- cross@geneticData[[1]]@founders[,11]
	cross2@geneticData[[1]]@hetData[[11]] <- cross@geneticData[[1]]@hetData[[11]]
	validObject(cross2, complete=TRUE)

	mapped <- new("mpcrossMapped", cross2, map = map)
	imputed <- imputeFounders(mapped)

	t <- table(imputed@geneticData[[1]]@imputed@data)
	expect_gt(sum(t[as.integer(names(t)) > 8]), length(imputed@geneticData[[1]]@imputed@data)/4)
})
