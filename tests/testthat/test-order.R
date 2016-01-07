context("Test ordering function")
test_that("Test that correct ordering is generated for an F2 population", 
	{
		f2Pedigree <- f2Pedigree(10000)
		map <- sim.map(len = 100, n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:101))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		ordered <- orderCross(grouped)
		correlated <- cor(match(markers(ordered), names(map[[1]])), 1:101)
		#This isn't due to numerical accuracy. The ordering is NOT perfect, so the correlation is not 1 or -1. 
		expect_equal(abs(correlated), 1, tolerance = 1e-3)
	})
test_that("Test that correct ordering is generated for an F2 population with two chromosomes", 
	{
		f2Pedigree <- f2Pedigree(10000)
		map <- sim.map(len = rep(100, 2), n.mar = 101, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross <- subset(cross, markers = sample(1:202))
		rf <- estimateRF(cross)
		grouped <- formGroups(rf, groups = 2, method = "average", clusterBy = "theta")
		ordered <- orderCross(grouped)

		correlationChromosome1 <- cor(match(names(map[[1]]), markers(ordered)), 1:101)
		correlationChromosome2 <- cor(match(names(map[[2]]), markers(ordered)), 1:101)
		expect_equal(abs(correlationChromosome1), 1, tolerance = 1e-3)
		expect_equal(abs(correlationChromosome2), 1, tolerance = 1e-3)
	})
