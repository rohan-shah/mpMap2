context("Test formGroups")
test_that("Check that formGroups can be applied multiple times", 
	{
		f2Pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = rep(100, 2), n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)
		formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = "average", clusterBy = "theta")
		expect_warning(formGroups <- formGroups(mpcrossRF = formGroups, groups = 2, method = "average", clusterBy = "theta"), NA)
		formGroups@rf <- NULL
		expect_that(formGroups <- formGroups(mpcrossRF = formGroups, groups = 2, method = "average", clusterBy = "theta"), throws_error("cannot be automatically converted to an object of class mpcrossRF"))
	})
test_that("Check that formGroups requires lod entry if clusterBy is \"combined\" or \"lod\"",
	{
		f2Pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = rep(100, 2), n.mar = 10, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross)

		for(method in c("average", "complete", "single"))
		{
			for(clusterBy in c("combined", "lod"))
			{
				expect_that(formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy), throws_error("must have a @rf@lod entry"))
				expect_that(formGroupsPreCluster <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy, preCluster = TRUE), throws_error("must have a @rf@lod entry"))
			}
			formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = "theta")
			formGroupsPreCluster <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = "theta", preCluster = TRUE)
			expect_null(formGroups@lg@imputedTheta)
			expect_null(formGroupsPreCluster@lg@imputedTheta)
		}
	})
test_that("Check that formGroups works for an f2 design",
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = rep(100, 2), n.mar = 100, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)

		for(method in c("average", "complete", "single"))
		{
			for(clusterBy in c("combined", "theta", "lod"))
			{
				formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy)
				formGroupsPreCluster <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy, preCluster = TRUE)
				expect_identical(formGroups@lg, formGroupsPreCluster@lg)
				expect_identical(length(formGroups@lg@allGroups), 2L)
				chr1Group <- formGroups@lg@groups[names(map[[1]])[1]]
				chr2Group <- formGroups@lg@groups[names(map[[2]])[1]]
				expectedGroups <- c(rep(chr1Group, 100), rep(chr2Group, 100))
				names(expectedGroups) <- markers(cross)

				expect_identical(formGroups@lg@groups, expectedGroups)
				expect_null(formGroups@lg@imputedTheta)
			}
		}
	})
test_that("Check that formGroups works for an f2 design with small chromosomes",
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = rep(5, 2), n.mar = 100, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)

		for(method in c("average", "complete", "single"))
		{
			for(clusterBy in c("combined", "theta", "lod"))
			{
				formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy)
				formGroupsPreCluster <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy, preCluster = TRUE)
				expect_identical(formGroups@lg, formGroupsPreCluster@lg)
				expect_identical(length(formGroups@lg@allGroups), 2L)
				chr1Group <- formGroups@lg@groups[names(map[[1]])[1]]
				chr2Group <- formGroups@lg@groups[names(map[[2]])[1]]
				expectedGroups <- c(rep(chr1Group, 100), rep(chr2Group, 100))
				names(expectedGroups) <- markers(cross)

				expect_identical(formGroups@lg@groups, expectedGroups)
				expect_null(formGroups@lg@imputedTheta)
			}
		}
	})

test_that("Check that formGroups works for a ril design",
	{
		f2Pedigree <- rilPedigree(1000, selfingGenerations = 10)
		map <- qtl::sim.map(len = rep(100, 2), n.mar = 100, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		suppressWarnings(rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE))

		for(method in c("average", "complete", "single"))
		{
			for(clusterBy in c("combined", "theta", "lod"))
			{
				formGroups <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy)
				formGroupsPreCluster <- formGroups(mpcrossRF = rf, groups = 2, method = method, clusterBy = clusterBy, preCluster = TRUE)
				expect_identical(formGroups@lg, formGroupsPreCluster@lg)
				expect_identical(length(formGroups@lg@allGroups), 2L)
				chr1Group <- formGroups@lg@groups[names(map[[1]])[1]]
				chr2Group <- formGroups@lg@groups[names(map[[2]])[1]]
				expectedGroups <- c(rep(chr1Group, 100), rep(chr2Group, 100))
				names(expectedGroups) <- markers(cross)

				expect_identical(formGroups@lg@groups, expectedGroups)
				expect_null(formGroups@lg@imputedTheta)
			}
		}
	})
