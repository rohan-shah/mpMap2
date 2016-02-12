context("Test imputation")
pedigree <- f2Pedigree(10000)
map <- sim.map(len = rep(100, 3), n.mar = 11, anchor.tel = TRUE, include.x = FALSE, eq.spacing=TRUE)
cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)

test_that("Check that impute warns when re-imputing",
	{
		grouped <- formGroups(rf, groups = 3, clusterBy = "theta", method = "average")
		imputed <- impute(grouped)
		expect_warning(imputed <- impute(imputed), "Existing imputation data wil be overwritten")
	})
test_that("Check that impute generates an error if there is no recombination fraction data", 
	{
		grouped <- formGroups(rf, groups = 3, clusterBy = "theta", method = "average")
		grouped@rf <- NULL
		expect_that(impute(grouped), throws_error("did not contain recombination fraction"))
	})
test_that("Check that impute just makes copies of parts of the diagonal, if there are no NA values",
	{
		grouped <- formGroups(rf, groups = 3, clusterBy = "theta", method = "average")
		imputed <- impute(grouped)
		for(group in 1:3)
		{
			expect_identical(imputed@lg@imputedTheta[[group]], subset(imputed@rf@theta, markers = which(imputed@lg@groups == group)))
		}
	})
test_that("Check that impute copies the most closely linked data",
	{
		copiedMap <- map
		#Ensure that markers 1 and 2 are the same
		copiedMap[[1]][1:2] <- 0:1
		cross <- simulateMPCross(map=copiedMap, pedigree=pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		grouped <- formGroups(rf, groups = 1, clusterBy = "theta", method = "average")

		groupedWithNA <- grouped
		thetaAsMatrix <- as(grouped@rf@theta, "matrix")
		thetaAsMatrix[(-1):(-3),1] <- thetaAsMatrix[1,(-1):(-3)] <- NA
		groupedWithNA@rf@theta <- as(thetaAsMatrix, "rawSymmetricMatrix")
		imputed <- impute(groupedWithNA)
		expect_identical(imputed@lg@imputedTheta[[1]][4:11, 1, drop=TRUE], imputed@lg@imputedTheta[[1]][4:11, 2, drop=TRUE])
		expect_identical(imputed@lg@imputedTheta[[1]][1, 4:11, drop=TRUE], imputed@lg@imputedTheta[[1]][2, 4:11, drop=TRUE])
	})
test_that("Check that imputed marker data is discarded if the linkage groups are recalculated",
	{
		#Ensure that markers 1 and 2 are the same
		copiedMap <- map
		copiedMap[[1]][1:2] <- 0:1
		cross <- simulateMPCross(map=copiedMap, pedigree=pedigree, mapFunction = haldane)
		rf <- estimateRF(cross, keepLod = TRUE, keepLkhd = TRUE)
		grouped <- formGroups(rf, groups = 3, clusterBy = "theta", method = "average")
		imputed <- impute(grouped)

		grouped <- formGroups(imputed, groups = 2, clusterBy = "theta", method = "average")
		expect_identical(length(grouped@lg@allGroups), 2L)
		expect_null(grouped@lg@imputedTheta)
	})
rm(pedigree, cross, rf, map)
