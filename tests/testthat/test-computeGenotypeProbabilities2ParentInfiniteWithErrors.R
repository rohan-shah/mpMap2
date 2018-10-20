context("genotype probability computation, two parents, infinite selfing, with errors")
test_that("Test zero generations of intercrossing, with errors, no extra positions",
	{
		sampleSize <- 500
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- rilPedigree(populationSize = sampleSize, selfingGenerations = 6)
		pedigree@selfing <- "infinite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
		mapped <- new("mpcrossMapped", cross, map = map)
		mapped@geneticData[[1]]@finals[,"D1M26"] <- -mapped@geneticData[[1]]@finals[,"D1M26"] + 3
		suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
		genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
			{
				apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
			})
		genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))

		#The error model should compensate for the flipped marker
		expect_gt(cor(genotypesFromProbabilities[,"D1M25"], cross@geneticData[[1]]@finals[,"D1M25"], method = "spearman", use = "complete.obs"), 0.91)
		correct <- sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"])))
		expect_gt(correct / sampleSize, 0.91)
		expect_gt(cor(genotypesFromProbabilities[,"D1M27"], cross@geneticData[[1]]@finals[,"D1M27"], method = "spearman", use = "complete.obs"), 0.91)
		
		expect_gt(cor(genotypesFromProbabilities[,"D1M26"], genotypesFromProbabilities[,"D1M25"], method = "spearman"), 0.89)
		expect_gt(cor(genotypesFromProbabilities[,"D1M26"], genotypesFromProbabilities[,"D1M27"], method = "spearman"), 0.89)
		expect_true(all(result@geneticData[[1]]@probabilities@data[1:10,1:20] != 1 & result@geneticData[[1]]@probabilities@data[1:10,1:20] != 0))
	})
test_that("Test non-zero generations of intercrossing, with errors, no extra positions",
	{
		sampleSize <- 500
		testFunc <- function(map, pedigree)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- -mapped@geneticData[[1]]@finals[,"D1M26"] + 3
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))
			#The error model should compensate for the flipped marker
			expect_gt(cor(genotypesFromProbabilities[,"D1M25"], cross@geneticData[[1]]@finals[,"D1M25"], method = "spearman", use = "complete.obs"), 0.86)
			correct <- sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"])))
			expect_gt(correct / sampleSize, 0.91)
			expect_gt(cor(genotypesFromProbabilities[,"D1M27"], cross@geneticData[[1]]@finals[,"D1M27"], method = "spearman", use = "complete.obs"), 0.86)

			expect_gt(cor(genotypesFromProbabilities[,"D1M26"], genotypesFromProbabilities[,"D1M25"], method = "spearman"), 0.86)
			expect_gt(cor(genotypesFromProbabilities[,"D1M26"], genotypesFromProbabilities[,"D1M27"], method = "spearman"), 0.86)
			expect_true(all(result@geneticData[[1]]@probabilities@data[1:10,1:20] != 1 & result@geneticData[[1]]@probabilities@data[1:10,1:20] != 0))
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "infinite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "infinite"
		pedigrees <- list(pedigree1, pedigree2)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
test_that("Test zero generations of intercrossing, with errors, with extra positions",
	{
		sampleSize <- 500

		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- rilPedigree(populationSize = sampleSize, selfingGenerations = 6)
		pedigree@selfing <- "infinite"
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
		mapped <- new("mpcrossMapped", cross, map = map)
		mapped@geneticData[[1]]@finals[,"D1M26"] <- -mapped@geneticData[[1]]@finals[,"D1M26"] + 3
		suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
		genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
			{
				apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
			})
		genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
		colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))

		#The error model should compensate for the flipped marker
		expect_gt(cor(genotypesFromProbabilities[,"D1M25"], cross@geneticData[[1]]@finals[,"D1M25"], method = "spearman", use = "complete.obs"), 0.91)
		correct <- sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"])))
		expect_gt(correct / sampleSize, 0.91)
		expect_gt(cor(genotypesFromProbabilities[,"D1M27"], cross@geneticData[[1]]@finals[,"D1M27"], method = "spearman", use = "complete.obs"), 0.91)

		expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman"), 0.91)
		expect_gt(cor(genotypesFromProbabilities[,"D1M27"], genotypesFromProbabilities[,"extra"], method = "spearman"), 0.91)
	})
test_that("Test non-zero generations of intercrossing, with errors, with extra positions",
	{
		sampleSize <- 1000
		testFunc <- function(map, pedigree)
		{
			cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + removeHets()
			mapped <- new("mpcrossMapped", cross, map = map)
			mapped@geneticData[[1]]@finals[,"D1M26"] <- -mapped@geneticData[[1]]@finals[,"D1M26"] + 3
			suppressWarnings(result <- computeGenotypeProbabilities(mapped, extraPositions = list("1" = c("extra" = 25.5)), errorProb = 0.05))
			genotypesFromProbabilities <- lapply(1:nLines(result), function(x)
				{
					apply(result@geneticData[[1]]@probabilities@data[(2*x-1):(2*x),], 2, which.max)
				})
			genotypesFromProbabilities <- do.call(rbind, genotypesFromProbabilities)
			colnames(genotypesFromProbabilities) <- unlist(lapply(result@geneticData[[1]]@probabilities@map, names))

			#The error model should compensate for the flipped marker
			expect_gt(cor(genotypesFromProbabilities[,"D1M25"], cross@geneticData[[1]]@finals[,"D1M25"], method = "spearman", use = "complete.obs"), 0.85)
			correct <- sum(diag(table(genotypesFromProbabilities[,"D1M26"], cross@geneticData[[1]]@finals[,"D1M26"])))
			expect_gt(correct / sampleSize, 0.89)
			expect_gt(cor(genotypesFromProbabilities[,"D1M27"], cross@geneticData[[1]]@finals[,"D1M27"], method = "spearman", use = "complete.obs"), 0.85)

			expect_gt(cor(genotypesFromProbabilities[,"extra"], genotypesFromProbabilities[,"D1M26"], method = "spearman"), 0.85)
			expect_gt(cor(genotypesFromProbabilities[,"D1M27"], genotypesFromProbabilities[,"extra"], method = "spearman"), 0.85)
		}
		map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree1 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 1)
		pedigree1@selfing <- "infinite"
		pedigree2 <- twoParentPedigree(initialPopulationSize = sampleSize, selfingGenerations = 6, nSeeds = 1, intercrossingGenerations = 2)
		pedigree2@selfing <- "infinite"
		pedigrees <- list(pedigree1, pedigree2)
		for(pedigree in pedigrees)
		{
			testFunc(map, pedigree)
		}
	})
