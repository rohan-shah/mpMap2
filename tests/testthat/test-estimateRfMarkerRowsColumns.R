context("markerRows and markerColumns arguments of estimateRF")

test_that("Check markerRows and markerColumns",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = c(50, 50), n.mar = 11, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		rfAll <- estimateRF(cross, keepLod = TRUE, keepLkhd=TRUE)
		rfBetweenChromosomes <- estimateRF(cross, keepLod = TRUE, keepLkhd=TRUE, markerRows = 1:11, markerColumns = 12:22)

		expect_identical(rfBetweenChromosomes@rf@theta[12:22,1:11], rfAll@rf@theta[12:22,1:11])
		expect_identical(rfBetweenChromosomes@rf@theta[1:11,12:22], rfAll@rf@theta[1:11,12:22])

		expect_identical(rfBetweenChromosomes@rf@lod[12:22,1:11], rfAll@rf@lod[12:22,1:11])
		expect_identical(rfBetweenChromosomes@rf@lod[1:11,12:22], rfAll@rf@lod[1:11,12:22])

		expect_identical(rfBetweenChromosomes@rf@lkhd[12:22,1:11], rfAll@rf@lkhd[12:22,1:11])
		expect_identical(rfBetweenChromosomes@rf@lkhd[1:11,12:22], rfAll@rf@lkhd[1:11,12:22])
	})

test_that("Check iterative computation using markerRows and markerColumns",
	{
		pedigree <- f2Pedigree(100)
		map <- qtl::sim.map(len = c(50, 50), n.mar = 11, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane) + biparentalDominant()
		rfAll <- estimateRF(cross, keepLod = TRUE, keepLkhd=TRUE)
		rf1 <- estimateRF(cross, keepLod = TRUE, keepLkhd=TRUE, markerRows = 1:11, markerColumns = 1:22)
		suppressWarnings(rf2 <- estimateRF(rf1, keepLod = TRUE, keepLkhd=TRUE, markerRows = 12:22, markerColumns = 1:22))

		expect_identical(rf2@rf@lod, rfAll@rf@lod)
		expect_identical(rf2@rf@lkhd, rfAll@rf@lkhd)
		expect_identical(rf2@rf@theta, rfAll@rf@theta)
	})
