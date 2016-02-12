context("Test openmp support")
test_that("Check that we have the same results for an f2, with and without openmp ",
	{
		f2Pedigree <- f2Pedigree(10)
		map <- sim.map(len = 100, n.mar = 4000, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		.Call("omp_set_num_threads", 1, PACKAGE="mpMap2")
		rf <- estimateRF(cross, gbLimit = 1)
		.Call("omp_set_num_threads", 4, PACKAGE="mpMap2")
		rf2 <- estimateRF(cross, gbLimit = 1)
		expect_identical(rf, rf2)

		f2Pedigree <- f2Pedigree(100)
		map <- sim.map(len = 100, n.mar = 400, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)

		.Call("omp_set_num_threads", 1, PACKAGE="mpMap2")
		rf <- estimateRF(cross, gbLimit = 1)
		.Call("omp_set_num_threads", 4, PACKAGE="mpMap2")
		rf2 <- estimateRF(cross, gbLimit = 1)
		expect_identical(rf, rf2)

	})
