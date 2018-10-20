context("lg validation")

pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
pedigree@selfing <- "finite"
map <- qtl::sim.map(len = rep(100, 2), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
rf <- estimateRF(cross)
grouped <- formGroups(rf, groups=2, clusterBy = "theta")
lg <- grouped@lg

test_that("Simulated lg passes validation",
	{
		expect_identical(validObject(lg, complete=TRUE), TRUE)
	})
test_that("Slot groups cannot contain NA values",
	{
		copied <- lg
		copied@groups[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot allGroups cannot contain NA values",
	{
		copied <- lg
		copied@allGroups[1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot allGroups must contain nonnegative values",
	{
		copied <- lg
		copied@groups[1] <- -1L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot groups must contain values that are also in allGroups",
	{
		copied <- lg
		copied@groups[1] <- 99999L
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot imputedTheta must have the correct length",
	{
		copied <- lg
		copied@imputedTheta <- list()
		expect_that(validObject(copied, complete=TRUE), throws_error("had the wrong length"))
		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 1)
		expect_that(validObject(copied, complete=TRUE), throws_error("had the wrong length"))
		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 2)
		expect_that(validObject(copied, complete=TRUE), throws_error("list of rawSymmetricMatrix objects"))
	})
test_that("Slot imputedTheta must have objects of the correct class",
	{
		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 2)
		expect_that(validObject(copied, complete=TRUE), throws_error("list of rawSymmetricMatrix objects"))
		names(copied@imputedTheta) <- c("1", "2")
		copied@imputedTheta[["1"]] <- copied@imputedTheta[["2"]] <- lg
		expect_that(validObject(copied, complete=TRUE), throws_error("list of rawSymmetricMatrix objects"))
	})
test_that("Slot imputedTheta must have objects with the right number of markers",
	{
		nMarkers1 <- sum(lg@groups == 1)
		nMarkers2 <- sum(lg@groups == 2)
		marker1Names <- names(which(lg@groups == 1))
		marker2Names <- names(which(lg@groups == 2))
		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 2)
		names(copied@imputedTheta) <- c("1", "2")
		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer((nMarkers1-1)*nMarkers1/2)), markers = marker1Names[1:(nMarkers1-1)], levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer((nMarkers2-1)*nMarkers2/2)), markers = marker2Names[1:(nMarkers2-1)], levels = rf@rf@theta@levels)
		expect_that(validObject(copied, complete=TRUE), throws_error("imputedTheta contained objects with the wrong length"))

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer((nMarkers1-1)*nMarkers1/2)), markers = marker1Names[1:(nMarkers1-1)], levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names[1:nMarkers2], levels = rf@rf@theta@levels)
		expect_that(validObject(copied, complete=TRUE), throws_error("imputedTheta contained objects with the wrong length"))

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names[1:nMarkers1], levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer((nMarkers2-1)*nMarkers2/2)), markers = marker2Names[1:(nMarkers2-1)], levels = rf@rf@theta@levels)
		expect_that(validObject(copied, complete=TRUE), throws_error("imputedTheta contained objects with the wrong length"))

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = rf@rf@theta@levels)
		expect_error(validObject(copied, complete=TRUE), NA)
	})
test_that("Slot imputedTheta must have objects with the correct marker names",
	{
		marker1Names <- names(which(lg@groups == 1))
		marker2Names <- names(which(lg@groups == 2))
		nMarkers1 <- sum(lg@groups == 1)
		nMarkers2 <- sum(lg@groups == 2)

		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 2)
		names(copied@imputedTheta) <- c("1", "2")
		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = rf@rf@theta@levels)
		copied@imputedTheta[["1"]]@markers[1] <- "invalidMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error("object@imputedTheta were inconsistent with those in slot object@groups"))
	})
test_that("Slot imputedTheta must have objects with the correct levels",
	{
		marker1Names <- names(which(lg@groups == 1))
		marker2Names <- names(which(lg@groups == 2))
		nMarkers1 <- sum(lg@groups == 1)
		nMarkers2 <- sum(lg@groups == 2)

		copied <- lg
		copied@imputedTheta <- vector(mode = "list", length = 2)
		names(copied@imputedTheta) <- c("1", "2")
		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = 0)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = rf@rf@theta@levels)
		expect_that(validObject(copied, complete=TRUE), throws_error("levels must be the same"))

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = 0)
		expect_that(validObject(copied, complete=TRUE), throws_error("levels must be the same"))

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = rf@rf@theta@levels)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = rf@rf@theta@levels)
		expect_error(validObject(copied, complete=TRUE), NA)

		copied@imputedTheta[["1"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers1*(nMarkers1+1)/2)), markers = marker1Names, levels = 0)
		copied@imputedTheta[["2"]] <- new("rawSymmetricMatrix", data = as.raw(integer(nMarkers2*(nMarkers2+1)/2)), markers = marker2Names, levels = 0)
		expect_error(validObject(copied, complete=TRUE), NA)

	})
test_that("Slot imputedTheta must have the correct names",
	{
		imputed <- impute(grouped)
		copied <- imputed@lg
	
		names(copied@imputedTheta) <- NULL
		expect_that(validObject(copied), throws_error())

		names(copied@imputedTheta)[1] <- "a"
		expect_that(validObject(copied), throws_error())
	})
rm(pedigree, map, cross, rf, grouped, lg)
