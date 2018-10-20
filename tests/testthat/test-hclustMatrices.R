context("Test hclust matrix calls to construct distance matrices")
functionNames <- c("hclustThetaMatrix", "hclustCombinedMatrix", "hclustLodMatrix")
requiresLod <- c("hclustCombinedMatrix", "hclustLodMatrix")
test_that("hclustMatrices throws appropriate exceptions",
{
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(10)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	#Throws, wrong class
	for(functionName in functionNames)
	{
		expect_that(.Call(functionName, cross, as.list(1:11), PACKAGE="mpMap2"), throws_error())
	}
	#Throws, wrong type for second argument
	for(functionName in functionNames)
	{
		expect_that(.Call(functionName, cross, 1:11, PACKAGE="mpMap2"), throws_error())
	}

	rf <- estimateRF(cross)

	#Throws, not all markers present
	for(functionName in functionName)
	{
		expect_that(.Call(functionName, rf, as.list(1:10), PACKAGE="mpMap2"), throws_error())
	}

	#Throws, markers are duplicated
	for(functionName in functionName)
	{
		expect_that(.Call(functionName, rf, as.list(c(1, 1:11)), PACKAGE="mpMap2"), throws_error())
	}
	for(functionName in functionName)
	{
		expect_that(.Call(functionName, rf, list(1:5, 5:11), PACKAGE="mpMap2"), throws_error())
	}
	#Throws because Lod was not calculated
	for(functionName in requiresLod)
	{
		expect_that(.Call(functionName, rf, as.list(1:11), PACKAGE="mpMap2"), throws_error())
	}


	rf <- estimateRF(cross, keepLod = TRUE)
	#Doesn't throw
	for(functionName in functionName)
	{
		.Call(functionName, rf, as.list(1:11), PACKAGE="mpMap2")
	}

})
test_that("hclustThetaMatrix call gives expected results",
{
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(10)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	
	#This takes the average of the first two theta values
	newTheta <- .Call("hclustThetaMatrix", rf, c(list(1:2), 3:11), PACKAGE="mpMap2")
	oldTheta <- rf@rf@theta[1:11,1:11]
	averagedPart <- apply(oldTheta[1:2,-(1:2)], 2, mean)
	expectedNewTheta <- rbind(c(mean(oldTheta[1:2,1:2]), averagedPart), cbind(averagedPart, oldTheta[3:11, 3:11]))
	expectedNewTheta <- as.dist(expectedNewTheta)
	expectedNewTheta <- unclass(expectedNewTheta)
	attributes(expectedNewTheta) <- NULL
	expect_equal(newTheta, expectedNewTheta)

	#This takes the average of the first and last values
	newTheta <- .Call("hclustThetaMatrix", rf, c(list(c(1,11)), 2:10), PACKAGE="mpMap2")
	oldTheta <- rf@rf@theta[1:11,1:11]
	averagedPart <- apply(oldTheta[c(1,11),-c(1,11)], 2, mean)
	expectedNewTheta <- rbind(c(mean(oldTheta[c(1,11), c(1,11)]), averagedPart), cbind(averagedPart, oldTheta[2:10, 2:10]))
	expectedNewTheta <- as.dist(expectedNewTheta)
	expectedNewTheta <- unclass(expectedNewTheta)
	attributes(expectedNewTheta) <- NULL
	expect_equal(newTheta, expectedNewTheta)

	#This takes the average of everything
	newTheta <- .Call("hclustThetaMatrix", rf, c(list(1:11)), PACKAGE="mpMap2")
	expectedNewTheta <- mean(rf@rf@theta[1:11,1:11])
	expectedNewTheta <- as.dist(expectedNewTheta)
	expectedNewTheta <- unclass(expectedNewTheta)
	attributes(expectedNewTheta) <- NULL
	expect_equal(newTheta, expectedNewTheta)

	#This takes the average of the first ten
	newTheta <- .Call("hclustThetaMatrix", rf, c(list(1:10), 11), PACKAGE="mpMap2")
	expectedNewTheta <- rbind(c(mean(rf@rf@theta[1:10,1:10]), mean(rf@rf@theta[1:10,11])), c(mean(rf@rf@theta[1:10,11]), rf@rf@theta[11,11]))
	expectedNewTheta <- as.dist(expectedNewTheta)
	expectedNewTheta <- unclass(expectedNewTheta)
	attributes(expectedNewTheta) <- NULL
	expect_equal(newTheta, expectedNewTheta)
	
	#This reverses the order
	newTheta <- .Call("hclustThetaMatrix", rf, as.list(11:1), PACKAGE="mpMap2")
	expectedNewTheta <- rf@rf@theta[11:1,11:1]
	expectedNewTheta <- as.dist(expectedNewTheta)
	expectedNewTheta <- unclass(expectedNewTheta)
	attributes(expectedNewTheta) <- NULL
	expect_equal(newTheta, expectedNewTheta)
	boolean <- lower.tri(rf@rf@theta[1:11,1:11], diag=FALSE)
	expect_equal(newTheta, rf@rf@theta[1:11,1:11][11:1,11:1][boolean])
})
test_that("hclustLodMatrix call gives expected results",
{
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(10)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross, keepLod = TRUE)
	maxLod <- max(rf@rf@lod)
	
	#This takes the average of the first two lod values
	newLod <- .Call("hclustLodMatrix", rf, c(list(1:2), 3:11), PACKAGE="mpMap2")
	oldLod <- maxLod - as(rf@rf@lod, "matrix")
	averagedPart <- apply(oldLod[1:2,-(1:2)], 2, mean)
	expectedNewLod <- rbind(c(mean(oldLod[1:2,1:2]), averagedPart), cbind(averagedPart, oldLod[3:11, 3:11]))
	expectedNewLod <- as.dist(expectedNewLod)
	expectedNewLod <- unclass(expectedNewLod)
	attributes(expectedNewLod) <- NULL
	expect_equal(newLod, expectedNewLod)

	#This takes the average of the first and last values
	newLod <- .Call("hclustLodMatrix", rf, c(list(c(1,11)), 2:10), PACKAGE="mpMap2")
	averagedPart <- apply(oldLod[c(1,11),-c(1,11)], 2, mean)
	expectedNewLod <- rbind(c(mean(oldLod[c(1,11), c(1,11)]), averagedPart), cbind(averagedPart, oldLod[2:10, 2:10]))
	expectedNewLod <- as.dist(expectedNewLod)
	expectedNewLod <- unclass(expectedNewLod)
	attributes(expectedNewLod) <- NULL
	expect_equal(newLod, expectedNewLod)

	#This takes the average of everything
	newLod <- .Call("hclustLodMatrix", rf, c(list(1:11)), PACKAGE="mpMap2")
	expectedNewLod <- mean(oldLod)
	expectedNewLod <- as.dist(expectedNewLod)
	expectedNewLod <- unclass(expectedNewLod)
	attributes(expectedNewLod) <- NULL
	expect_equal(newLod, expectedNewLod)

	#This takes the average of the first ten
	newLod <- .Call("hclustLodMatrix", rf, c(list(1:10), 11), PACKAGE="mpMap2")
	expectedNewLod <- rbind(c(mean(oldLod[1:10,1:10]), mean(oldLod[1:10,11])), c(mean(oldLod[1:10,11]), oldLod[11,11]))
	expectedNewLod <- as.dist(expectedNewLod)
	expectedNewLod <- unclass(expectedNewLod)
	attributes(expectedNewLod) <- NULL
	expect_equal(newLod, expectedNewLod)
	
	#This reverses the order
	newLod <- .Call("hclustLodMatrix", rf, as.list(11:1), PACKAGE="mpMap2")
	expectedNewLod <- oldLod[11:1,11:1]
	expectedNewLod <- as.dist(expectedNewLod)
	expectedNewLod <- unclass(expectedNewLod)
	attributes(expectedNewLod) <- NULL
	expect_equal(newLod, expectedNewLod)
	boolean <- lower.tri(rf@rf@theta[1:11,1:11], diag=FALSE)
	expect_equal(newLod, oldLod[11:1,11:1][boolean])
})
test_that("hclustCombinedMatrix call gives expected results",
{
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(10)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross, keepLod = TRUE)
	
	minimumDistance <- min(diff(rf@rf@theta@levels))
	theta <- rf@rf@theta[1:11,1:11]
	lod <- as(rf@rf@lod, "matrix")
	maxLod <- max(rf@rf@lod)
	combined <- theta + ((maxLod - lod)/max(lod))*minimumDistance
	
	#This takes the average of the first two combined values
	newMatrix <- .Call("hclustCombinedMatrix", rf, c(list(1:2), 3:11), PACKAGE="mpMap2")
	averagedPart <- apply(combined[1:2,-(1:2)], 2, mean)
	expectedNewMatrix <- rbind(c(mean(combined[1:2,1:2]), averagedPart), cbind(averagedPart, combined[3:11, 3:11]))
	expectedNewMatrix <- as.dist(expectedNewMatrix)
	expectedNewMatrix <- unclass(expectedNewMatrix)
	attributes(expectedNewMatrix) <- NULL
	expect_equal(newMatrix, expectedNewMatrix)

	#This takes the average of the first and last values
	newMatrix <- .Call("hclustCombinedMatrix", rf, c(list(c(1,11)), 2:10), PACKAGE="mpMap2")
	averagedPart <- apply(combined[c(1,11),-c(1,11)], 2, mean)
	expectedNewMatrix <- rbind(c(mean(combined[c(1,11), c(1,11)]), averagedPart), cbind(averagedPart, combined[2:10, 2:10]))
	expectedNewMatrix <- as.dist(expectedNewMatrix)
	expectedNewMatrix <- unclass(expectedNewMatrix)
	attributes(expectedNewMatrix) <- NULL
	expect_equal(newMatrix, expectedNewMatrix)

	#This takes the average of everything
	newMatrix <- .Call("hclustCombinedMatrix", rf, c(list(1:11)), PACKAGE="mpMap2")
	expectedNewMatrix <- mean(combined)
	expectedNewMatrix <- as.dist(expectedNewMatrix)
	expectedNewMatrix <- unclass(expectedNewMatrix)
	attributes(expectedNewMatrix) <- NULL
	expect_equal(newMatrix, expectedNewMatrix)

	#This takes the average of the first ten
	newMatrix <- .Call("hclustCombinedMatrix", rf, c(list(1:10), 11), PACKAGE="mpMap2")
	expectedNewMatrix <- rbind(c(mean(combined[1:10,1:10]), mean(combined[1:10,11])), c(mean(combined[1:10,11]), combined[11,11]))
	expectedNewMatrix <- as.dist(expectedNewMatrix)
	expectedNewMatrix <- unclass(expectedNewMatrix)
	attributes(expectedNewMatrix) <- NULL
	expect_equal(newMatrix, expectedNewMatrix)
	
	#This reverses the order
	newMatrix <- .Call("hclustCombinedMatrix", rf, as.list(11:1), PACKAGE="mpMap2")
	expectedNewMatrix <- combined[11:1,11:1]
	expectedNewMatrix <- as.dist(expectedNewMatrix)
	expectedNewMatrix <- unclass(expectedNewMatrix)
	attributes(expectedNewMatrix) <- NULL
	expect_equal(newMatrix, expectedNewMatrix)
	boolean <- lower.tri(combined, diag=FALSE)
	expect_equal(newMatrix, combined[11:1,11:1][boolean])
})

rm(functionNames)
