context("Test ordering function")
mpMap2::omp_set_num_threads(1)
test_that("Test that not having rf data generates an error",
	{
		f2Pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:101))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		imputed <- impute(grouped)

		grouped@rf <- NULL
		expect_that(orderCross(grouped), throws_error("did not contain recombination fraction"))
		imputed@rf <- NULL
		orderCross(imputed)
	})
test_that("Test that correct ordering is generated for an F2 population", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:101))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		ordered <- orderCross(grouped)
		correlated <- cor(match(markers(ordered), names(map[[1]])), 1:101)
		#This isn't due to numerical accuracy. The ordering is NOT perfect, so the correlation is not 1 or -1. 
		expect_equal(abs(correlated), 1, tolerance = 1e-1)
	})
test_that("Test that identical orderings are generate with and without imputedTheta slot", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:101))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		set.seed(1)
		ordered <- orderCross(grouped)

		imputed <- impute(grouped)
		imputed@rf <- NULL
		set.seed(1)
		imputed <- orderCross(imputed)
		expect_identical(markers(imputed), markers(ordered))
	})
test_that("Test that correct ordering is generated for an F2 population with two chromosomes", 
	{
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = rep(100, 2), n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:202))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 2, method = "average", clusterBy = "theta")
		ordered <- orderCross(grouped)

		correlationChromosome1 <- cor(match(names(map[[1]]), markers(ordered)), 1:101)
		correlationChromosome2 <- cor(match(names(map[[2]]), markers(ordered)), 1:101)
		expect_equal(abs(correlationChromosome1), 1, tolerance = 1e-1)
		expect_equal(abs(correlationChromosome2), 1, tolerance = 1e-1)
	})
test_that("Test that we can order a single marker group",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = 100, n.mar = 1, anchor.tel=FALSE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, verbose = FALSE)
		grouped <- new("mpcrossLG", rf, rf = rf@rf, lg = new ("lg", groups = c("D1M1" = 1L), allGroups = 1L))
		imputed <- impute(grouped)
		expect_error(ordered <- orderCross(imputed), NA)
	})
