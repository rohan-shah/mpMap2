context("Conversion of 8-way cross to mpwgaim object")
suppressWarnings(capture.output(couldLoadPackages <- require(mpwgaim, quietly=T) && require(mpMap, quietly=T)))
if(couldLoadPackages)
{
	test_that("Require probabilities, with type = mpMarker", 
	{
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, nSeeds = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len=rep(200,2), n.mar=rep(101,2), eq.spacing=TRUE, include.x=FALSE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + multiparentSNP(keepHets = FALSE)
		mapped <- new("mpcrossMapped", cross2, map = map)
		expect_error(suppressWarnings(newConverted <- as.mpInterval(mapped, positions = "all", type = "mpMarker")), "not have genotype probabilities")

		prob <- computeGenotypeProbabilities(mapped, extraPositions = generateGridPositions(10))
		suppressWarnings(newConverted <- as.mpInterval(prob, positions = "all", type = "mpMarker"))
		
		prob <- computeGenotypeProbabilities(mapped)
		suppressWarnings(newConverted <- as.mpInterval(prob, positions = "all", type = "mpMarker"))

		prob <- computeGenotypeProbabilities(mapped, extraPositions = generateGridPositions(10))
		newConverted <- as.mpInterval(prob, positions = "marker", type = "mpMarker")
		
		prob <- computeGenotypeProbabilities(mapped)
		newConverted <- as.mpInterval(prob, positions = "marker", type = "mpMarker")
	})
	test_that("Conversion of 8-way object agrees with old code, with type = mpMarker",
	{
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, nSeeds = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len=rep(200,2), n.mar=rep(101,2), eq.spacing=TRUE, include.x=FALSE)

		#Add duplicated positions
		map[[1]] <- c("Extra1" = 0, map[[1]], "Extra2" = 200)
		map[[2]] <- c("Extra3" = 0, map[[2]], "Extra4" = 200)
		
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + multiparentSNP(keepHets = FALSE)
		mapped <- new("mpcrossMapped", cross2, map = map)
		#This error probability is chosen to match the one in mpcross2int. It seems the one in mpcross2int can't be easily changed. 
		prob <- computeGenotypeProbabilities(mapped)
		newConverted <- as.mpInterval(prob, positions = "marker", type = "mpMarker")

		old <- list()
		old$finals <- finals(cross2)
		old$founders <- founders(cross2)
		old$pedigree <- sim.mpped(nfounders=8, nfunnels=1, nperfam=100, nssdgen = 10)
		old$map <- map
		old$id <- which(old$pedigree[,4] == 1)
		old$fid <- 1:8
		class(old) <- "mpcross"
		suppressWarnings(capture.output(oldConverted <- mpcross2int(old, gen.type="mpMarker", method = "qtl", geprob = 0)))

		#They may disagree on the marker encoding, so fix that up
		for(chr in 1:2)
		{
			sapply(1:ncol(newConverted$geno[[chr]]$founders), function(x)
			{
				if(newConverted$geno[[chr]]$founders[1, x] != oldConverted$geno[[chr]]$founders[1, x])
				{
					newConverted$geno[[chr]]$founders[,x] <<- -newConverted$geno[[chr]]$founders[,x]
					newConverted$geno[[chr]]$data[,x] <<- -newConverted$geno[[chr]]$data[,x]
				}
			})
		}
		expect_identical(newConverted$nfounders, 8L)
		for(chr in 1:2)
		{
			expect_identical(oldConverted$geno[[chr]]$dist, newConverted$geno[[chr]]$dist)
			expect_identical(oldConverted$geno[[chr]]$map, newConverted$geno[[chr]]$map)
			expect_true(all.equal(oldConverted$geno[[chr]]$founders, newConverted$geno[[chr]]$founders, check.attributes = FALSE))
			expect_true(all.equal(oldConverted$geno[[chr]]$data, newConverted$geno[[chr]]$data, check.attributes = FALSE))
		}
	})
	test_that("Conversion removes co-located markers, with type = mpInterval",
	{
		pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, nSeeds = 1, intercrossingGenerations = 0)
		map <- qtl::sim.map(len=rep(200,2), n.mar=rep(101,2), eq.spacing=TRUE, include.x=FALSE)
		map[[1]] <- c("added_1_1" = 0, "added_1_2" = 0, map[[1]])
		map[[2]] <- c("added_2_1" = 0, "added_2_2" = 0, "added_2_3" = 0, map[[2]])

		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + removeHets()
		mapped <- new("mpcrossMapped", cross2, map = map)
		#This error probability is chosen to match the one in mpcross2int. It seems the one in mpcross2int can't be easily changed. 
		expect_warning(newConverted <- as.mpInterval(mapped, type = "mpInterval", errorProb = 0.01, homozygoteMissingProb = 1, heterozygoteMissingProb = 1))

		expectedMarkers <- length(unique(map[[1]]))
		expectedIntervals <- length(unique(map[[1]])) - 1L
		expectedColumns <- expectedIntervals*8L
		expect_identical(ncol(newConverted$geno[[1]]$founders), expectedMarkers)
		expect_identical(length(newConverted$geno[[1]]$map), expectedMarkers)
		expect_identical(ncol(newConverted$geno[[1]]$intval), expectedColumns)
		expect_identical(ncol(newConverted$geno[[1]]$data), expectedMarkers)

		expect_identical(ncol(newConverted$geno[[2]]$founders), expectedMarkers)
		expect_identical(length(newConverted$geno[[2]]$map), expectedMarkers)
		expect_identical(ncol(newConverted$geno[[2]]$intval), expectedColumns)
		expect_identical(ncol(newConverted$geno[[2]]$data), expectedMarkers)
	})
}
