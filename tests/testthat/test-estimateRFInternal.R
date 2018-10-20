context("mpMap2:::estimateRFInternal validation tests")
test_that("Arguments must have correct types",
	{
		#recombValues must be numeric
		expect_that(mpMap2:::estimateRFInternal(object = NULL, recombValues = list(), lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("recombinationFractions must be"))

		#recombValues must contain 0 and 0.5
		expect_that(mpMap2:::estimateRFInternal(object = NULL, recombValues = 0.0, lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("recombinationFractions did not contain"))
		expect_that(mpMap2:::estimateRFInternal(object = NULL, recombValues = 0.5, lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("recombinationFractions did not contain"))

		#Input object must be an S4
		expect_that(mpMap2:::estimateRFInternal(object = NULL, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("object must be an S4"))

		#Input markerRows must be an integer vector
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = 100, n.mar = 4000, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = list(), markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("markerRows must be an integer"))

		#Input markerColumns must be an integer
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = 1:2, markerColumns = list(), keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("markerColumns must be an integer"))

		#Input lineWeights must be a list
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = 1:2, markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), throws_error("lineWeights must be a list"))

		#gbLimit must be a single number
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = 1:2, verbose = NULL), throws_error("gbLimit must be a single numeric"))
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = c(1,2), verbose = NULL), throws_error("gbLimit must be a single numeric"))

		#lineWeights must have the right length
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = -1, verbose = NULL), throws_error("lineWeights had the wrong number of entries"))

		#keepLod must be a boolean
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = list(), keepLkhd = NULL, gbLimit = -1, verbose = NULL), throws_error("keepLod must be a boolean"))

#keepLkhd must be a boolean
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = list(), gbLimit = -1, verbose = NULL), throws_error("keepLkhd must be a boolean"))

		#verbose must be a boolean
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list()), throws_error("verbose must be a"))

		#There must be at least one design
		copiedCross <- cross
		copiedCross@geneticData <- new("geneticDataList")
		expect_that(mpMap2:::estimateRFInternal(object = copiedCross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("There must be at least one design"))

		#Must be at least one markerRow value
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = as.integer(c()), markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("markerRows must have at least one entry"))

		#Must be at least one markerCloumn value
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = as.integer(c()), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("markerColumns must have at least one entry"))

		#Entry of lineWeights has the wrong lengthi
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:10, markerColumns = 1:10, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("lineWeights had the wrong length"))

		#Cannot have negative values for markerRows
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = -(1:2), markerColumns = -(1:2), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("Invalid values for input markerRows"))

		#Cannot have negative values for markerColumns
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = 1:2, markerColumns = -(1:2), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("Invalid value for input markerColumns"))

		#Number of values to estimate cannot be zero
		expect_that(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = 1000, markerColumns = 1, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), throws_error("give a region that is contained in the lower triangular part"))
	})
