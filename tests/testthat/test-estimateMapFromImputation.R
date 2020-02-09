context("Test estimateMapFromImputation")
test_that("Check that estimation of gap sizes is approximately correct, for four parent designs", 
	{
		set.seed(1)
		pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 2000, selfingGenerations = 5, intercrossingGenerations = 0)
		pedigree@selfing <- "finite"
		map <- qtl::sim.map(len = rep(100, 1), n.mar = rep(101, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + multiparentSNP(keepHets = TRUE)
		cross <- subset(cross, markers = c(1:30, 70:100))
		capture.output(rf <- estimateRF(cross))
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		estimated.map <- estimateMap(grouped, maxOffset = 10)

		mapped <- new("mpcrossMapped", grouped, map = estimated.map)
		mapped@geneticData[[1]]@pedigree@selfing <- "infinite"
		suppressWarnings(imputed <- imputeFounders(mapped, errorProb = 0.01))
		reestimated <- estimateMapFromImputation(imputed)
		expect_equal(reestimated@map[[1]]["D1M70"] - reestimated@map[[1]]["D1M30"], 40, tolerance = 0.06, check.attributes = FALSE)
	})
test_that("Check that estimation of gap sizes is approximately correct, for eight parent designs", 
	{
		pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 3000, selfingGenerations = 5, intercrossingGenerations = 0)
		pedigree@selfing <- "infinite"
		map <- qtl::sim.map(len = rep(100, 1), n.mar = rep(101, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + multiparentSNP(keepHets = FALSE)
		cross <- subset(cross, markers = c(1:30, 70:100))
		capture.output(rf <- estimateRF(cross))
		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
		estimated.map <- estimateMap(grouped, maxOffset = 10)

		mapped <- new("mpcrossMapped", grouped, map = estimated.map)
		mapped@geneticData[[1]]@pedigree@selfing <- "infinite"
		suppressWarnings(imputed <- imputeFounders(mapped, errorProb = 0.01))
		reestimated <- estimateMapFromImputation(imputed)
		expect_equal(reestimated@map[[1]]["D1M70"] - reestimated@map[[1]]["D1M30"], 40, tolerance = 0.1, check.attributes = FALSE)
	})
#test_that("Check that estimation of gap sizes is approximately correct, for eight parent designs", 
#	{
#		pedigree <- sixteenParentPedigreeRandomFunnels(initialPopulationSize = 2000, selfingGenerations = 5, intercrossingGenerations = 0)
#		pedigree@selfing <- "infinite"
#		map <- qtl::sim.map(len = rep(100, 1), n.mar = rep(101, 1), anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
#		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + multiparentSNP(keepHets = TRUE)
#		cross <- subset(cross, markers = c(1:10, 90:100))
#		rf <- estimateRFSingleDesign(cross, verbose = TRUE)
#		grouped <- formGroups(rf, groups = 1, method = "average", clusterBy = "theta")
#		estimated.map <- estimateMap(grouped, maxOffset = 10)
#
#		mapped <- new("mpcrossMapped", grouped, map = estimated.map)
#		imputed <- imputeFounders(mapped, errorProb = 0.1)
#		reestimated <- estimateMapFromImputation(imputed)
#		expect_equal(reestimated@map[[1]]["D1M90"] - reestimated@map[[1]]["D1M10"], 80, tolerance = 0.02, check.attributes = FALSE)
#	})
