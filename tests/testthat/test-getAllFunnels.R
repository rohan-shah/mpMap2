context("getAllFunnels")

test_that("Checking that code to get all funnels works, for F2",
	{
		map <- qtl::sim.map(len = c(10), n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(10)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_identical(getAllFunnels(cross, standardised = FALSE), cbind(rep(1L, 10), rep(2L, 10)))
		expect_identical(getAllFunnels(cross, standardised = TRUE), cbind(rep(1L, 10), rep(2L, 10)))
	})
test_that("Checking that code to get all funnels works, for 4-way",
	{
		map <- qtl::sim.map(len = c(10), n.mar = 2, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, intercrossingGenerations = 0)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)

		expect_identical(getAllFunnels(cross, standardised = FALSE), cbind(rep(1L, 10), rep(2L, 10), rep(3L, 10), rep(4L, 10)))
		expect_identical(getAllFunnels(cross@geneticData[[1]], standardised = FALSE), cbind(rep(1L, 10), rep(2L, 10), rep(3L, 10), rep(4L, 10)))
		
		expect_identical(getAllFunnels(cross, standardised = TRUE), cbind(rep(1L, 10), rep(2L, 10), rep(3L, 10), rep(4L, 10)))
		expect_identical(getAllFunnels(cross@geneticData[[1]], standardised = TRUE), cbind(rep(1L, 10), rep(2L, 10), rep(3L, 10), rep(4L, 10)))

		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 10, selfingGenerations = 1, intercrossingGenerations = 1)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		naMatrix <- rep(as.integer(NA), 4*10)
		dim(naMatrix) <- c(10, 4)
		
		expect_identical(getAllFunnels(cross, standardised = FALSE), naMatrix)
		expect_identical(getAllFunnels(cross@geneticData[[1]], standardised = FALSE), naMatrix)
		expect_identical(getAllFunnels(cross, standardised = TRUE), naMatrix)
		expect_identical(getAllFunnels(cross@geneticData[[1]], standardised = TRUE), naMatrix)

		pedigree <- fourParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0)
		
		#simulated pedigrees are always standardised, so mix things up a bit. 
		newMother <- pedigree@mother
		newMother[seq(5, 10, by = 2)] <- pedigree@mother[seq(6, 10, by = 2)]
		newMother[seq(6, 10, by = 2)] <- pedigree@mother[seq(5, 10, by = 2)]
		pedigree@mother <- newMother
		
		newFather <- pedigree@father
		newFather[seq(5, 10, by = 2)] <- pedigree@father[seq(6, 10, by = 2)]
		newFather[seq(6, 10, by = 2)] <- pedigree@father[seq(5, 10, by = 2)]
		pedigree@father <- newFather

		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		expect_false(identical(getAllFunnels(cross, standardised = FALSE), getAllFunnels(cross, standardised = TRUE)))

		funnelsStandardised <- unique(apply(getAllFunnels(cross, standardised = TRUE), 1, function(x) sum(x*(4^(0:3)))))
		funnelsNotStandardised <- unique(apply(getAllFunnels(cross, standardised = FALSE), 1, function(x) sum(x*(4^(0:3)))))
		#Check that we only have three different funnels in the standardised version. 
		expect_identical(length(funnelsStandardised), 3L)
		#We have three funnels in the non-standardised version, but they're different. 
		expect_identical(length(funnelsNotStandardised), 3L)
		expect_identical(length(intersect(funnelsStandardised, funnelsNotStandardised)), 0L)
	})
