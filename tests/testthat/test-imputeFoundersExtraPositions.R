context("Founder imputation with extra positions")
test_that("Test that the positions are put in the correct order",
	{
		map <- sim.map(len = c(100, 100), n.mar = c(101, 101), anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(10)
		pedigree@selfing <- "finite"

		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		suppressWarnings(result <- imputeFounders(mapped, errorProb = 0, extraPositions = list("2" = c("C2P2" = 60.5, "C2P1" = 50.5), "1" = c("C1P2" = 40.5, "C1P1" = 30.5))))

		names <- colnames(result@geneticData[[1]]@imputed@data)

		index <- match("C1P1", names)
		expect_equal(names[index-1], "D1M31")
		expect_equal(names[index+1], "D1M32")

		index <- match("C1P2", names)
		expect_equal(names[index-1], "D1M41")
		expect_equal(names[index+1], "D1M42")

		index <- match("C2P1", names)
		expect_equal(names[index-1], "D2M51")
		expect_equal(names[index+1], "D2M52")

		index <- match("C2P2", names)
		expect_equal(names[index-1], "D2M61")
		expect_equal(names[index+1], "D2M62")
	})
