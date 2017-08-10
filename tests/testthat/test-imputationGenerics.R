context("imputationGenerics")

test_that("Test that flatImputationMapNames works",
{
	pedigree1 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped1_", lineNames(pedigree1))
	pedigree1@selfing <- "finite"
	map <- qtl::sim.map(len = c(10, 10), n.mar = 11, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE)
	cross1 <- simulateMPCross(map = map, pedigree = pedigree1, mapFunction = "haldane", seed = 1) + multiparentSNP(keepHets = TRUE)

	mapped1 <- new("mpcrossMapped", cross1, map = map)
	expect_error(flatImputationMapNames(mapped1))

	imputed1 <- imputeFounders(mapped1)
	expect_identical(unlist(lapply(imputed1@geneticData[[1]]@imputed@map, names)), flatImputationMapNames(imputed1))

	pedigree2 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped2_", lineNames(pedigree2))
	pedigree2@selfing <- "finite"
	cross2 <- simulateMPCross(map = map, pedigree = pedigree2, mapFunction = "haldane", seed = 2) + multiparentSNP(keepHets = TRUE)

	mapped2 <- new("mpcrossMapped", cross2, map = map)
	imputed2 <- imputeFounders(mapped2, extraPositions = generateGridPositions(1))
	expect_identical(unlist(lapply(imputed2@geneticData[[1]]@imputed@map, names)), flatImputationMapNames(imputed2))
	
	combined <- imputed1 + imputed2
	expect_identical(list(unlist(lapply(imputed1@geneticData[[1]]@imputed@map, names)), unlist(lapply(imputed2@geneticData[[1]]@imputed@map, names))), flatImputationMapNames(combined))
})
test_that("Test that imputationMap works",
{
	pedigree1 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped1_", lineNames(pedigree1))
	pedigree1@selfing <- "finite"
	map <- qtl::sim.map(len = c(10, 10), n.mar = 11, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE)
	cross1 <- simulateMPCross(map = map, pedigree = pedigree1, mapFunction = "haldane", seed = 1) + multiparentSNP(keepHets = TRUE)

	mapped1 <- new("mpcrossMapped", cross1, map = map)
	expect_error(imputationMap(mapped1))

	imputed1 <- imputeFounders(mapped1)
	expect_identical(imputed1@geneticData[[1]]@imputed@map, imputationMap(imputed1))

	pedigree2 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped2_", lineNames(pedigree2))
	pedigree2@selfing <- "finite"
	cross2 <- simulateMPCross(map = map, pedigree = pedigree2, mapFunction = "haldane", seed = 2) + multiparentSNP(keepHets = TRUE)

	mapped2 <- new("mpcrossMapped", cross2, map = map)
	imputed2 <- imputeFounders(mapped2, extraPositions = generateGridPositions(1))
	expect_identical(imputed2@geneticData[[1]]@imputed@map, imputationMap(imputed2))
	
	combined <- imputed1 + imputed2
	expect_identical(list(imputed1@geneticData[[1]]@imputed@map, imputed2@geneticData[[1]]@imputed@map), imputationMap(combined))
})
test_that("Test that imputationData works",
{
	pedigree1 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped1_", lineNames(pedigree1))
	pedigree1@selfing <- "finite"
	map <- qtl::sim.map(len = c(10, 10), n.mar = 11, anchor.tel = TRUE, include.x = FALSE, eq.spacing = TRUE)
	cross1 <- simulateMPCross(map = map, pedigree = pedigree1, mapFunction = "haldane", seed = 1) + multiparentSNP(keepHets = TRUE)

	mapped1 <- new("mpcrossMapped", cross1, map = map)
	expect_error(imputationData(mapped1))

	imputed1 <- imputeFounders(mapped1)
	expect_identical(imputed1@geneticData[[1]]@imputed@data, imputationData(imputed1))

	pedigree2 <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
	lineNames(pedigree1) <- paste0("ped2_", lineNames(pedigree2))
	pedigree2@selfing <- "finite"
	cross2 <- simulateMPCross(map = map, pedigree = pedigree2, mapFunction = "haldane", seed = 2) + multiparentSNP(keepHets = TRUE)

	mapped2 <- new("mpcrossMapped", cross2, map = map)
	imputed2 <- imputeFounders(mapped2, extraPositions = generateGridPositions(1))
	expect_identical(imputed2@geneticData[[1]]@imputed@data, imputationData(imputed2))
	
	combined <- imputed1 + imputed2
	expect_identical(list(imputed1@geneticData[[1]]@imputed@data, imputed2@geneticData[[1]]@imputed@data), imputationData(combined))
})
