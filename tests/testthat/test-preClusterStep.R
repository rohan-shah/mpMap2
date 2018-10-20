context("Test pre cluster call to aggregate markers")

test_that("preClusterStep throws appropriate exceptions",
{
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(500)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	#Throws, wrong class
	expect_that(.Call("preClusterStep", cross, PACKAGE="mpMap2"), throws_error())
	rf <- estimateRF(cross)
	#Doesn't throw
	preCluster <- .Call("preClusterStep", rf, PACKAGE="mpMap2")

	rf@rf@theta@levels <- c(0.1, 0.5)
	#Throws, 0, not contained
	expect_that(preCluster <- .Call("preClusterStep", rf, PACKAGE="mpMap2"), throws_error())
})
test_that("preClusterStep call gives expected results",
{
	#We only need the rawSymmetricMatrix object, but it's easiest to construct it this way. 
	map <- qtl::sim.map(len = 100, n.mar = 11, anchor.tel=TRUE, include.x=FALSE, eq.spacing=TRUE)
	f2Pedigree <- f2Pedigree(1)
	cross <- simulateMPCross(map=map, pedigree=f2Pedigree, mapFunction = haldane)
	rf <- estimateRF(cross)
	rf@rf@theta@data[] <- as.raw(0xff)
	rf@rf@theta@data[1:3] <- as.raw(0)
	expect_identical(.Call("preClusterStep", rf, PACKAGE="mpMap2"), c(as.list(3:11), list(1:2)))

	rf@rf@theta@data[] <- as.raw(0xff)
	rf@rf@theta@data[c(6,9, 10)] <- as.raw(0)
	expect_identical(.Call("preClusterStep", rf, PACKAGE="mpMap2"), c(as.list(c(1:2, 5:11)), list(3:4)))

	rf@rf@theta@data[] <- as.raw(0xff)
	rf@rf@theta@data[c(1:3, 6, 9, 10)] <- as.raw(0)
	expect_identical(.Call("preClusterStep", rf, PACKAGE="mpMap2"), c(as.list(5:11), list(1:2, 3:4)))

	rf@rf@theta@data[1:10] <- as.raw(0)
	expect_identical(.Call("preClusterStep", rf, PACKAGE="mpMap2"), c(as.list(5:11), list(1:4)))

	for(i in c(4:5, 7:8))
	{
		copied <- rf
		copied@rf@theta@data[i] <- as.raw(0xff)
		expect_identical(.Call("preClusterStep", copied, PACKAGE="mpMap2"), c(as.list(5:11), list(1:2, 3:4)))
	}

})
