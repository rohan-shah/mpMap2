context("toMpMap")

test_that("Checking that code to convert to mpMap works for 8-way object",
	{
		if(require(mpMap, quietly = TRUE))
		{
			map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
			pedigree <- sim.mpped(nfounders = 8, nfunnels = 1, nperfam = 100, nssdgen = 6, nseeds = 1, iripgen = 0, seed = 1)
			mpcross <- sim.mpcross(map = map, pedigree = pedigree, map.function = "haldane")

			#Warning about rownames of founders and finals being incorrect
			expect_warning(mpcross2 <- fromMpMap(mpcross))
			reverted <- toMpMap(mpcross2)

			expect_equivalent(reverted$founders, mpcross$founders)
			expect_equivalent(reverted$finals, mpcross$finals)
			expect_equivalent(reverted$pedigree, mpcross$pedigree)
			expect_identical(reverted$id, mpcross$id)
			expect_identical(reverted$fid, mpcross$fid)
			expect_identical(reverted$map, mpcross$map)
			expect_identical(attr(reverted, "type"), attr(mpcross, "type"))
		}
	})
test_that("Checking that code to convert to mpMap works for 4-way object",
	{
		if(require(mpMap, quietly = TRUE))
		{
			map <- qtl::sim.map(len = c(50, 50), n.mar = 51, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
			pedigree <- sim.mpped(nfounders = 4, nfunnels = 1, nperfam = 100, nssdgen = 6, nseeds = 1, iripgen = 0, seed = 1)
			mpcross <- sim.mpcross(map = map, pedigree = pedigree, map.function = "haldane")

			#Warning about rownames of founders and finals being incorrect
			expect_warning(mpcross2 <- fromMpMap(mpcross))
			reverted <- toMpMap(mpcross2)

			expect_equivalent(reverted$founders, mpcross$founders)
			expect_equivalent(reverted$finals, mpcross$finals)
			expect_equivalent(reverted$pedigree, mpcross$pedigree)
			expect_identical(reverted$id, mpcross$id)
			expect_identical(reverted$fid, mpcross$fid)
			expect_identical(reverted$map, mpcross$map)
			expect_identical(attr(reverted, "type"), attr(mpcross, "type"))
		}
	})
