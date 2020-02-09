context("Test subset function")

test_that("Checking subset on object of class mpcross by markers, with a single dataset",
	{
		#Test function for an object with a single dataset
		testSingle <- function(cross, subsetted, marker)
		{
			expect_identical(validObject(subsetted, complete=TRUE), TRUE)
			expect_identical(subsetted@geneticData[[1]]@finals, cross@geneticData[[1]]@finals[,marker,drop=F])
			expect_identical(subsetted@geneticData[[1]]@founders, cross@geneticData[[1]]@founders[,marker,drop=F])
			expect_identical(subsetted@geneticData[[1]]@pedigree, cross@geneticData[[1]]@pedigree)
			expect_identical(length(subsetted@geneticData[[1]]@hetData), 1L)
			expect_identical(cross@geneticData[[1]]@hetData[[marker]], subsetted@geneticData[[1]]@hetData[[1]])
		}
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		
		#Checks on marker indices
		expect_that(subset(cross, markers = -1), throws_error())
		expect_that(subset(cross, markers = (-3):0), throws_error())
		expect_that(subset(cross, markers = 1:12), throws_error())
		expect_error(subset(cross, markers = 1:11), NA)

		#Test without dominant markers
		subsetted <- subset(cross, markers = 2)
		testSingle(cross, subsetted, 2)

		#Tests with dominant markers
		crossDominant <- cross + biparentalDominant()
		subsetted <- subset(crossDominant, markers = 2)
		testSingle(crossDominant, subsetted, 2)
		subsetted <- subset(crossDominant, markers = 3)
		testSingle(crossDominant, subsetted, 3)
		subsetted <- subset(crossDominant, markers = 4)
		testSingle(crossDominant, subsetted, 4)
	})
test_that("Subset refuses to duplicate markers and lines",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(1000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		expect_that(subset(cross, markers = rep(1:nMarkers(cross), each = 2)), throws_error("Duplicates detected"))
		expect_that(subset(cross, lines = rep(rownames(cross@geneticData[[1]]@finals), each = 2)), throws_error("Duplicates detected"))
	})
test_that("Subset changes the pedigree from detailedPedigree to pedigree when subsetting by lines",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(1000)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		expect_that(cross@geneticData[[1]]@pedigree, is_a("detailedPedigree"))

		subsetted <- subset(cross, markers = 1:nMarkers(cross))
		expect_that(subsetted@geneticData[[1]]@pedigree, is_a("detailedPedigree"))

		subsetted <- subset(cross, lines = 1:2)
		expect_true(!is(subsetted@geneticData[[1]]@pedigree, "detailedPedigree"))
	})
test_that("Checking subset on object of class mpcross by markers, with two datasets",
	{
		map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		f2Pedigree <- f2Pedigree(1000)

		#Test function for an object with a pair of datasets
		testPair <- function(cross, subsetted, marker)
		{
			expect_identical(length(cross@geneticData), 2L)
			expect_identical(length(subsetted@geneticData), 2L)

			expect_identical(validObject(subsetted, complete=TRUE), TRUE)
			expect_identical(subsetted@geneticData[[1]]@finals, cross@geneticData[[1]]@finals[,marker,drop=F])
			expect_identical(subsetted@geneticData[[2]]@finals, cross@geneticData[[2]]@finals[,marker,drop=F])

			expect_identical(subsetted@geneticData[[1]]@founders, cross@geneticData[[1]]@founders[,marker,drop=F])
			expect_identical(subsetted@geneticData[[2]]@founders, cross@geneticData[[2]]@founders[,marker,drop=F])

			expect_identical(subsetted@geneticData[[1]]@pedigree, cross@geneticData[[1]]@pedigree)
			expect_identical(subsetted@geneticData[[2]]@pedigree, cross@geneticData[[2]]@pedigree)

			expect_identical(cross@geneticData[[1]]@hetData[[marker]], subsetted@geneticData[[1]]@hetData[[1]])
			expect_identical(cross@geneticData[[2]]@hetData[[marker]], subsetted@geneticData[[2]]@hetData[[1]])
		}
		cross1 <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		cross2 <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		pedigreeSubset <- cross1@geneticData[[1]]@pedigree@lineNames %in% rownames(cross1@geneticData[[1]]@finals)
		cross1@geneticData[[1]]@pedigree@lineNames[pedigreeSubset] <- paste0(cross1@geneticData[[1]]@pedigree@lineNames[pedigreeSubset], ",2")
		rownames(cross1@geneticData[[1]]@finals) <- paste0(rownames(cross1@geneticData[[1]]@finals), ",2")

		#Test codominant
		cross <- cross1 + cross2
		subsetted <- subset(cross, markers = 2)
		testPair(cross, subsetted, 2)

		#Test dominant
		crossDominant <- cross
		crossDominant@geneticData[[1]] <- crossDominant@geneticData[[1]] + biparentalDominant()
		crossDominant@geneticData[[2]] <- crossDominant@geneticData[[2]] + biparentalDominant()
		subsetted <- subset(crossDominant, markers = 2)
		testPair(crossDominant, subsetted, 2)
		subsetted <- subset(crossDominant, markers = 3)
		testPair(crossDominant, subsetted, 3)
		subsetted <- subset(crossDominant, markers = 4)
		testPair(crossDominant, subsetted, 4)
	})
test_that("Checking subset by lines when imputation data is present",
	{
		map <- qtl::sim.map(len = 100, n.mar = 101, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
		pedigree <- f2Pedigree(500)
		cross <- simulateMPCross(map=map, pedigree=pedigree, mapFunction = haldane)
		mapped <- new("mpcrossMapped", cross, map = map)
		suppressWarnings(result <- imputeFounders(mapped, errorProb = 0))
		subsetted <- subset(result, lines = sample(lineNames(result)))
		expect_error(validObject(subsetted, complete = TRUE), NA)
	})
