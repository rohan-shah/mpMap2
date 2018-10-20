context("Founder imputation with extra positions")
test_that("Test that the positions are put in the correct order",
	{
		map <- qtl::sim.map(len = c(100, 100), n.mar = c(101, 101), anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
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
test_that("Test that imputations at positions look right for F2 populations",
	{
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map[[1]] <- map[[1]][c(1:20, 80:101)]
		pedigree <- f2Pedigree(500)
		pedigree@selfing <- "finite"

		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		suppressWarnings(result <- imputeFounders(mapped, errorProb = 0, extraPositions = list("1" = c("P1" = 30, "P2" = 40, "P3" = 50, "P4" = 60, "P5" = 70))))

		#Different homozygotes
		indices1 <- which(result@geneticData[[1]]@imputed@data[,20] == 1 & result@geneticData[[1]]@imputed@data[,26] == 2)
		indices2 <- which(result@geneticData[[1]]@imputed@data[,20] == 2 & result@geneticData[[1]]@imputed@data[,26] == 1)

		#Check that everything between two homozygotes is imputed as a hetrozygote. 
		expect_equal(as.integer(result@geneticData[[1]]@imputed@data[indices1,c("P1", "P2", "P3", "P4", "P5")]), rep(c(3L, 3L, 3L, 3L, 2L), each = length(indices1)))
		expect_equal(as.integer(result@geneticData[[1]]@imputed@data[indices2,c("P1", "P2", "P3", "P4", "P5")]), rep(c(3L, 3L, 3L, 3L, 1L), each = length(indices2)))

		#Check that everything between two identical homozygotes is imputed the same as the end-points
		indices <- which(result@geneticData[[1]]@imputed@data[,20] == 1 & result@geneticData[[1]]@imputed@data[,26] == 1)
		expect_equal(as.integer(result@geneticData[[1]]@imputed@data[indices,c("P1", "P2", "P3", "P4", "P5")]), rep(1L, length(indices)*5))

		indices <- which(result@geneticData[[1]]@imputed@data[,20] == 2 & result@geneticData[[1]]@imputed@data[,26] == 2)
		expect_equal(as.integer(result@geneticData[[1]]@imputed@data[indices,c("P1", "P2", "P3", "P4", "P5")]), rep(2L, length(indices)*5))
	})
test_that("Test that imputations at positions look right for 4-way population using a single funnel",
	{
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		map[[1]] <- map[[1]][c(1:20, 80:101)]
		pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 1000, selfingGenerations = 1, intercrossingGenerations = 0)
		pedigree@selfing <- "finite"

		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		suppressWarnings(result <- imputeFounders(mapped, errorProb = 0, extraPositions = list("1" = c("P1" = 30, "P2" = 40, "P3" = 50, "P4" = 60, "P5" = 70))))

		#Certain different homozygotes should be seperated by hetrozygotes - E.g. (1, 3), (1, 4), (2, 3), (2, 4). 
		indices <- which((result@geneticData[[1]]@imputed@data[,20] == 1 | result@geneticData[[1]]@imputed@data[,20] == 2) & (result@geneticData[[1]]@imputed@data[,26] == 3 | result@geneticData[[1]]@imputed@data[,26] == 4))
		indices <- c(indices, which((result@geneticData[[1]]@imputed@data[,20] == 3 | result@geneticData[[1]]@imputed@data[,20] == 4) & (result@geneticData[[1]]@imputed@data[,26] == 1 | result@geneticData[[1]]@imputed@data[,26] == 2)))
		expect_equal(as.logical(apply(result@geneticData[[1]]@imputed@data[indices,c("P1", "P2", "P3", "P4", "P5")], 1, function(x) any(x >= 5))), rep(TRUE, length(indices)))

		#But some hetrozygotes are impossible in this pedigree. In that case the differente homozygotes are NOT seperated by a run of hets. 
		indices <- which((result@geneticData[[1]]@imputed@data[,20] == 1 & result@geneticData[[1]]@imputed@data[,26] == 2) | (result@geneticData[[1]]@imputed@data[,20] == 2 & result@geneticData[[1]]@imputed@data[,26] == 1))
		indices <- c(indices, which((result@geneticData[[1]]@imputed@data[,20] == 3 & result@geneticData[[1]]@imputed@data[,26] == 4) | (result@geneticData[[1]]@imputed@data[,20] == 4 & result@geneticData[[1]]@imputed@data[,26] == 3)))
		#Check that everything between is imputed as a hetrozygote
		expect_equal(as.logical(result@geneticData[[1]]@imputed@data[indices,c("P1", "P2", "P3", "P4", "P5")] <= 4), rep(TRUE, length(indices)*5))

		#The same homozygote
		indices <- which(result@geneticData[[1]]@imputed@data[,20] <= 4 & result@geneticData[[1]]@imputed@data[,26] <= 4 & result@geneticData[[1]]@imputed@data[,20] == result@geneticData[[1]]@imputed@data[,26])
		#Check that everything between two identical homozygotes is imputed as the same homozygote
		expectedValue <- as.integer(do.call(rbind, lapply(as.list(result@geneticData[[1]]@imputed@data[indices,20]), function(x) rep(x, 5L))))
		expect_equal(as.integer(result@geneticData[[1]]@imputed@data[indices,c("P1", "P2", "P3", "P4", "P5")]), expectedValue)
	})
