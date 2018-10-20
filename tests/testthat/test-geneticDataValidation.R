context("geneticData validation")

#Set up a basic mpcross object
pedigree <- twoParentPedigree(initialPopulationSize=100, selfingGenerations=0, nSeeds=3, intercrossingGenerations=1)
map <- qtl::sim.map(len = rep(100, 1), n.mar = 12, anchor.tel = T, include.x=FALSE, sex.sp=FALSE, eq.spacing=T)
cross <- simulateMPCross(map = map, pedigree=pedigree, mapFunction = haldaneToRf, seed=1)
geneticData <- cross@geneticData[[1]]

test_that("Simulated geneticData object passes validation", 
	{
		expect_identical(validObject(geneticData), TRUE)
	})
test_that("geneticData allows numeric matrices", 
	{
		copied <- geneticData
		storage.mode(copied@founders) <- "numeric"
		storage.mode(copied@finals) <- "numeric"
		for(i in 1:length(copied@hetData)) storage.mode(copied@hetData[[i]]) <- "numeric"

		expect_identical(validObject(copied), TRUE)
	})
test_that("Slot founders must be numeric", 
	{
		copied <- geneticData
		storage.mode(copied@founders) <- "logical"
		expect_that(validObject(copied), throws_error())
	})
test_that("Slot founders cannot contain only NA", 
	{
		copied <- geneticData
		copied@founders[] <- NA
		expect_that(validObject(copied), throws_error())
	})
test_that("Slot founders cannot contain non-integer values", 
	{
		#floating point value in founders given an error
		copied <- geneticData
		copied@founders[1,1] <- 0.5
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slot finals must be numeric", 
	{
		copied <- geneticData
		storage.mode(copied@finals) <- "logical"
		expect_that(validObject(copied), throws_error())
	})
test_that("Slot finals cannot contain only NA", 
	{
		copied <- geneticData
		copied@finals[] <- NA
		expect_that(validObject(copied), throws_error())
	})
test_that("Slot finals cannot contain non-integer values",
	{
		#floating point value in finals gives an error
		copied <- geneticData
		copied@finals[1,1] <- 0.5
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
#Apparently setting the dimension wrong like this throws an error before the validObject call anyway
test_that("Slot founders must be two-dimensional",
	{
		copied <- geneticData
		expect_that(dim(copied@founders) <- c(2,3,4), throws_error())
	})
test_that("Slot finals must be two-dimensional",
	{
		copied <- geneticData
		expect_that(dim(copied@finals) <- c(2,3,4), throws_error())
	})
test_that("Slots founders must have dimnames",
	{
		copied <- geneticData
		dimnames(copied@founders) <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		expect_that(dimnames(copied@founders)[[1]] <- NULL, throws_error())
	})
test_that("Slots finals must have dimnames",
	{
		copied <- geneticData
		dimnames(copied@finals) <- NULL
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		expect_that(dimnames(copied@finals)[[1]] <- NULL, throws_error())
	})
test_that("Dimnames of founders cannot contain NA",
	{
		copied <- geneticData
		dimnames(copied@founders)[[1]][1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		dimnames(copied@founders)[[2]][1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Dimnames of finals cannot contain NA",
	{
		copied <- geneticData
		dimnames(copied@finals)[[1]][1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		dimnames(copied@finals)[[2]][1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Row and column names of slot founders must be unique",
	{
		copied <- geneticData
		colnames(copied@founders)[1] <- colnames(copied@founders)[2]
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		rownames(copied@founders)[1] <- rownames(copied@founders)[2]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Row and column names of slot finals must be unique",
	{
		copied <- geneticData
		colnames(copied@finals)[1] <- colnames(copied@finals)[2]
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		rownames(copied@finals)[1] <- rownames(copied@finals)[2]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Slots founders, finals and hetData must agree on the number of markers",
	{
		#Add marker to founders
		copied <- geneticData
		copied@founders <- cbind(copied@founders, newMarker = copied@founders[,1])
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Add marker to finals
		copied <- geneticData
		copied@finals <- cbind(copied@finals, newMarker = copied@finals[,1])
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Add marker to hetData
		copied <- geneticData
		copied@hetData[["newMarker"]] <- copied@hetData[[1]]
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Add marker to finals and founders
		copied <- geneticData
		copied@founders <- cbind(copied@founders, newMarker = copied@founders[,1])
		copied@finals <- cbind(copied@finals, newMarker = copied@finals[,1])
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Add marker to finals and hetData
		copied <- geneticData
		copied@hetData[["newMarker"]] <- copied@hetData[[1]]
		copied@finals <- cbind(copied@finals, newMarker = copied@finals[,1])
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Add marker to founders and hetData
		copied <- geneticData
		copied@founders <- cbind(copied@founders, newMarker = copied@founders[,1])
		copied@hetData[["newMarker"]] <- copied@hetData[[1]]
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Adding marker to all three is fine
		copied <- geneticData
		copied@founders <- cbind(copied@founders, newMarker = copied@founders[,1])
		copied@finals <- cbind(copied@finals, newMarker = copied@finals[,1])
		copied@hetData[["newMarker"]] <- copied@hetData[[1]]
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
test_that("Slots founders, finals and hetData must agree on the marker names",
	{
		copied <- geneticData
		colnames(copied@founders)[1] <- "newMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		colnames(copied@finals)[1] <- "newMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		names(copied@hetData)[1] <- "newMarker"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("All founders must be in the pedigree",
	{
		copied <- geneticData
		rownames(copied@founders)[1] <- "newFounder"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("All finals must be in the pedigree",
	{
		copied <- geneticData
		rownames(copied@finals)[1] <- "newFinal"
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("If slot pedigree has class detailedPedigree, then the founders must match the initials",
	{
		copied <- geneticData
		if(!inherits(copied@pedigree, "detailedPedigree")) stop("Slot pedigree should inherit from detailedPedigree")
		rownames(copied@founders)[1] <- copied@pedigree@lineNames[20]
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("If slot pedigree has class detailedPedigree, then the finals must match the observed",
	{
		copied <- geneticData
		if(!inherits(copied@pedigree, "detailedPedigree")) stop("Slot pedigree should inherit from detailedPedigree")
		rownames(copied@finals)[1] <- min(which(!copied@pedigree@observed))
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
#Now the stuff in the C code checked by alleleDataErrors
test_that("NA founders are allowed, so long as the hetData entry has 0 rows, and the finals are all NA",
	{
		copied <- geneticData
		copied@founders[1,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		copied@founders[1,1] <- NA
		copied@finals[,1] <- NA
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		copied@founders[1,1] <- NA
		copied@hetData[[1]] <- matrix(0, 0, 3)
		expect_that(validObject(copied, complete=TRUE), throws_error())

		copied <- geneticData
		copied@founders[1,1] <- NA
		copied@hetData[[1]] <- matrix(0, 0, 3)
		copied@finals[,1] <- NA
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
test_that("Every marker allele must be present in hetData", 
	{
		copied <- geneticData
		copied@hetData[[1]] <- matrix(0, 0, 3)
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Heterozygotes cannot be coded as homozygotes", 
	{
		copied <- geneticData
		copied@hetData[[1]][2:3,3] <- 1
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Values in first two columns of hetData matrix must be founder alleles", 
	{
		#Completely invalid values
		copied <- geneticData
		copied@hetData[[1]][1,1] <- max(copied@founders[,1])+1
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#More completely invalid values
		copied <- geneticData
		copied@hetData[[1]][1,2] <- max(copied@founders[,1])+1
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Case where both founders have a 1 allele, but the hetData contains 2 alleles
		copied <- geneticData
		copied@founders[,1] <- 1
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Values in finals must be valid according to hetData", 
	{
		#Remove heterozygotes from hetData
		copied <- geneticData
		indices <- copied@hetData[[1]][,1] == copied@hetData[[1]][,2]
		copied@hetData[[1]] <- copied@hetData[[1]][indices,]
		expect_that(validObject(copied, complete=TRUE), throws_error())

		#Completely invalid value
		copied <- geneticData
		copied@finals[1,1] <- 20
		expect_that(validObject(copied, complete=TRUE), throws_error())
	})
test_that("Arbitrary encodings are allowed for the founder alleles",
	{
		copied <- geneticData
		copied@founders[1,1] <- 100
		copied@founders[2,1] <- 201
		copied@finals[copied@finals[,1] == 1,1] <- 100
		copied@finals[copied@finals[,1] == 2,1] <- 201
		copied@hetData[[1]][copied@hetData[[1]] == 1] <- 100
		copied@hetData[[1]][copied@hetData[[1]] == 2] <- 201
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
test_that("Monomorphic markers are allowed",
	{
		copied <- geneticData
		copied@founders[,1] <- 1
		copied@finals[,1] <- 1
		copied@hetData[[1]] <- rbind(c(1,1,1))
		expect_identical(validObject(copied, complete=TRUE), TRUE)
	})
rm(pedigree, map, cross, geneticData)
