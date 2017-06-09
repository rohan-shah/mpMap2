context("Conversion of 8-way cross to mpwgaim object")
capture.output(couldLoadPackages <- suppressWarnings(require(mpwgaim, quietly=T) && require(mpMap, quietly=T)))
if(couldLoadPackages)
{
	test_that("Conversion of 4-way object agrees with old code, for SNP markers",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, nSeeds = 1, intercrossingGenerations = 0)
		map <- sim.map(len=rep(200,2), n.mar=rep(101,2), eq.spacing=TRUE, include.x=FALSE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + multiparentSNP(keepHets = FALSE)
		mapped <- new("mpcrossMapped", cross2, map = map)
		prob <- computeGenotypeProbabilities(mapped)
		newConverted <- as.mpInterval(prob)
		class(newConverted[[1]]) <- c("mpMarker", "cross", "interval")
		newConverted <- newConverted[[1]]

		old <- list()
		old$finals <- finals(cross2)
		old$founders <- founders(cross2)
		old$pedigree <- sim.mpped(nfounders=4, nfunnels=1, nperfam=100, nssdgen = 10)
		old$map <- map
		old$id <- which(old$pedigree[,4] == 1)
		old$fid <- 1:4
		class(old) <- "mpcross"
		capture.output(oldConverted <- mpcross2int(old, gen.type="mpMarker", method = "qtl"))

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
		expect_identical(newConverted$nfounders, 4L)
		for(chr in 1:2)
		{
			expect_identical(oldConverted$geno[[chr]]$dist, newConverted$geno[[chr]]$dist)
			expect_identical(oldConverted$geno[[chr]]$map, newConverted$geno[[chr]]$map)
			expect_true(all.equal(oldConverted$geno[[chr]]$founders, newConverted$geno[[chr]]$founders, check.attributes = FALSE))
			expect_true(all.equal(oldConverted$geno[[chr]]$data, newConverted$geno[[chr]]$data, check.attributes = FALSE))
			expect_true(all.equal(oldConverted$geno[[chr]]$imputed.dat, newConverted$geno[[chr]]$imputed.dat, check.attributes = FALSE, tolerance = 0.01))
		}
	})
	test_that("Conversion of 4-way object agrees with old code, for multiallelic markers",
	{
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 5, nSeeds = 1, intercrossingGenerations = 0)
		map <- sim.map(len=rep(200,2), n.mar=rep(101,2), eq.spacing=TRUE, include.x=FALSE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		cross2 <- cross + removeHets()
		mapped <- new("mpcrossMapped", cross2, map = map)
		prob <- computeGenotypeProbabilities(mapped)
		newConverted <- as.mpInterval(prob)
		class(newConverted[[1]]) <- c("mpMarker", "cross", "interval")
		newConverted <- newConverted[[1]]

		old <- list()
		old$finals <- finals(cross2)
		old$founders <- founders(cross2)
		old$pedigree <- sim.mpped(nfounders=4, nfunnels=1, nperfam=100, nssdgen = 10)
		old$map <- map
		old$id <- which(old$pedigree[,4] == 1)
		old$fid <- 1:4
		class(old) <- "mpcross"
		capture.output(oldConverted <- mpcross2int(old, gen.type="mpMarker", method = "qtl"))

		expect_identical(newConverted$nfounders, 4L)
		for(chr in 1:2)
		{
			expect_identical(oldConverted$geno[[chr]]$dist, newConverted$geno[[chr]]$dist)
			expect_identical(oldConverted$geno[[chr]]$map, newConverted$geno[[chr]]$map)
			expect_true(all.equal(oldConverted$geno[[chr]]$founders, newConverted$geno[[chr]]$founders + 1, check.attributes = FALSE))
			expect_true(all.equal(oldConverted$geno[[chr]]$data, newConverted$geno[[chr]]$data + 1, check.attributes = FALSE))
			expect_true(all.equal(oldConverted$geno[[chr]]$imputed.dat, newConverted$geno[[chr]]$imputed.dat, check.attributes = FALSE, tolerance = 0.01))
		}
	})
}
