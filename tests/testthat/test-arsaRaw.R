context("Test multithreading of arsaRaw")
if(class(try(.Call("omp_set_num_threads", 1, PACKAGE="mpMap2"), silent=TRUE)) != "try-error")
{
	test_that("Test that correct ordering is generated for an F2 population", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:101))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")

		.Call("omp_set_num_threads", 1, PACKAGE="mpMap2")
		orderedSingleThreaded <- orderCross(grouped)
		.Call("omp_set_num_threads", 16, PACKAGE="mpMap2")
		orderedMultiThreaded <- orderCross(grouped)
		
		correlationMultiThreaded <- cor(match(names(map[[1]]), markers(orderedMultiThreaded)), 1:101)
		correlationSingleThreaded <- cor(match(names(map[[1]]), markers(orderedSingleThreaded)), 1:101)
		expect_equal(abs(correlationMultiThreaded), 1, tolerance = 1e-1)
		expect_equal(abs(correlationSingleThreaded), 1, tolerance = 1e-1)
		.Call("omp_set_num_threads", 1, PACKAGE="mpMap2")
	})
}
