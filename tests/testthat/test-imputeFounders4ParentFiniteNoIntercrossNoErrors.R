context("Founder imputation, four parents, finite selfing, no intercrossing, no errors")
test_that("Test zero generations of intercrossing, with marker heterozygotes, without selfing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
	
			expect_identical(nrow(tmp), ncol(tmp))
			expect_identical(rownames(tmp), colnames(tmp))

			endHet <- nrow(tmp)
			expect_gt(sum(diag(tmp)) / sum(tmp), 0.965)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "finite"
		pedigree2 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "finite"

		pedigrees <- list(pedigree1, pedigree2)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}

	})
test_that("Test zero generations of intercrossing, with marker heterozygotes, with selfing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=TRUE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
	
			expect_identical(nrow(tmp), ncol(tmp))
			expect_identical(rownames(tmp), colnames(tmp))

			endHet <- nrow(tmp)
			expect_gt(sum(diag(tmp)) / sum(tmp), 0.965)
			#If a het is called, it should be correct with at least 95% chance
			expect_gt(sum(diag(tmp[5:endHet,5:endHet])) / sum(tmp[5:endHet,]), 0.95)
			#If a homozygote is called, it should be correct with a higher probability
			expect_gt(sum(diag(tmp[1:4,1:4])) / sum(tmp[1:4,1:4]), 0.95)
			#If a het is called, it should actually be a het (any het) with at least 99% chance
			expect_gt(sum(tmp[5:endHet,5:endHet]) / sum(tmp[5:endHet,]), 0.99)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "finite"
		pedigree2 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "finite"
		pedigree3 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
		pedigree3@selfing <- "finite"
		pedigree4 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
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
test_that("Test zero generations of intercrossing, without marker heterozygotes, without selfing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, heterozygoteMissingProb = 1, homozygoteMissingProb = 0.01, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			#Correct imputation rate should be 0.945
			expect_gt(sum(diag(tmp)) / sum(tmp), 0.93)
			expect_identical(nrow(tmp), ncol(tmp))
			expect_identical(rownames(tmp), colnames(tmp))

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "finite"
		pedigree2 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 0, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "finite"

		pedigrees <- list(pedigree1, pedigree2)
		for(map in maps)
		{
			for(pedigree in pedigrees)
			{
				testFunc(pedigree, map)
			}
		}

	})
test_that("Test zero generations of intercrossing, without marker heterozygotes, with selfing",
	{
		testFunc <- function(pedigree, map)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
			mapped <- new("mpcrossMapped", cross, map = map)
			suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))

			#Hetrozygotes will be discarded in imputation, which means that the imputed version won't be EXACTLY the same as the original data
			expect_identical(result@geneticData[[1]]@imputed@data, result@geneticData[[1]]@finals)

			cross2 <- cross + multiparentSNP(keepHets=FALSE)
			mapped <- new("mpcrossMapped", cross2, map = map)
			result <- imputeFounders(mapped, heterozygoteMissingProb = 1, homozygoteMissingProb = 0.01, errorProb = 0)
			tmp <- table(result@geneticData[[1]]@imputed@data, cross@geneticData[[1]]@finals)
			#Correct imputation rate should be 0.945
			expect_gt(sum(diag(tmp)) / sum(tmp), 0.945)
			#If a homozygote is called, it should be correct with 95% probability
			expect_gt(sum(diag(tmp[1:4,1:4])) / sum(tmp[1:4,1:4]), 0.95)
			expect_identical(nrow(tmp), ncol(tmp))
			expect_identical(rownames(tmp), colnames(tmp))
			endHet <- nrow(tmp)
			#If a het is called, it should actually be a het (any het) with at least 99% chance
			expect_gt(sum(tmp[5:endHet,5:endHet]) / sum(tmp[5:endHet,]), 0.99)

			expect_true(all(result@geneticData[[1]]@imputed@errors == 0))
		}
		map1 <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map2 <- qtl::sim.map(len = c(100, 100), n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		maps <- list(map1, map2)

		pedigree1 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree1@selfing <- "finite"
		pedigree2 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
		pedigree2@selfing <- "finite"
		pedigree3 <- fourParentPedigreeRandomFunnels(initialPopulationSize = 500, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
		pedigree3@selfing <- "finite"
		pedigree4 <- fourParentPedigreeSingleFunnel(initialPopulationSize = 500, selfingGenerations = 3, nSeeds = 1, intercrossingGenerations = 0)
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


