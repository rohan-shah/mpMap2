context("mpMap2:::estimateRFInternal validation tests")
test_that("Arguments must have correct types",
	{
		#recombValues must be numeric
		expect_error(mpMap2:::estimateRFInternal(object = NULL, recombValues = list(), lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "recombinationFractions must be", class = "Rcpp::not_compatible")

		#recombValues must contain 0 and 0.5
		expect_error(mpMap2:::estimateRFInternal(object = NULL, recombValues = 0.0, lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "recombinationFractions did not contain", class = "std::runtime_error")
		expect_error(mpMap2:::estimateRFInternal(object = NULL, recombValues = 0.5, lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "recombinationFractions did not contain", class = "std::runtime_error")

		#Input object must be an S4
		expect_error(mpMap2:::estimateRFInternal(object = NULL, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = NULL, markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "object must be an S4", class = "Rcpp::not_compatible")

		#Input markerRows must be an integer vector
		f2Pedigree <- f2Pedigree(1000)
		map <- qtl::sim.map(len = 100, n.mar = 4000, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
		cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = list(), markerColumns = NULL, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "markerRows must be an integer", class="Rcpp::not_compatible")

		#Input markerColumns must be an integer
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = NULL, markerRows = 1:2, markerColumns = list(), keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "markerColumns must be an integer", class="Rcpp::not_compatible")

		#Input lineWeights must be a list
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = 1:2, markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = NULL, verbose = NULL), "lineWeights must be a list", class="Rcpp::not_compatible")

		#gbLimit must be a single number
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = 1:2, verbose = NULL), "gbLimit must be a single numeric", class="Rcpp::not_compatible")
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = c(1,2), verbose = NULL), "gbLimit must be a single numeric", class="Rcpp::not_compatible")

		#lineWeights must have the right length
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = NULL, keepLkhd = NULL, gbLimit = -1, verbose = NULL), "lineWeights had the wrong number of entries", class="std::runtime_error")

		#keepLod must be a boolean
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = list(), keepLkhd = NULL, gbLimit = -1, verbose = NULL), "keepLod must be a boolean", class="std::runtime_error")

#keepLkhd must be a boolean
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = list(), gbLimit = -1, verbose = NULL), "keepLkhd must be a boolean", class="std::runtime_error")

		#verbose must be a boolean
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list()), "verbose must be a", class="std::runtime_error")

		#There must be at least one design
		copiedCross <- cross
		copiedCross@geneticData <- new("geneticDataList")
		expect_error(mpMap2:::estimateRFInternal(object = copiedCross, recombValues = c(0.0, 0.5), lineWeights = list(), markerRows = 1:2, markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "There must be at least one design", class="std::runtime_error")

		#Must be at least one markerRow value
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = as.integer(c()), markerColumns = 1:2, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "markerRows must have at least one entry", class="std::runtime_error")

		#Must be at least one markerCloumn value
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:2, markerColumns = as.integer(c()), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "markerColumns must have at least one entry", class="std::runtime_error")

		#Entry of lineWeights has the wrong lengthi
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(1), markerRows = 1:10, markerColumns = 1:10, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "lineWeights had the wrong length", class="std::runtime_error")

		#Cannot have negative values for markerRows
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = -(1:2), markerColumns = -(1:2), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "Invalid values for input markerRows", class="std::runtime_error")

		#Cannot have negative values for markerColumns
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = 1:2, markerColumns = -(1:2), keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "Invalid value for input markerColumns", class="std::runtime_error")

		#Number of values to estimate cannot be zero
		expect_error(mpMap2:::estimateRFInternal(object = cross, recombValues = c(0.0, 0.5), lineWeights = list(rep(1, 1000)), markerRows = 1000, markerColumns = 1, keepLod = TRUE, keepLkhd = TRUE, gbLimit = -1, verbose = list(verbose=FALSE, progressStyle = 1L)), "give a region that is contained in the lower triangular part", class="std::runtime_error")
	})
