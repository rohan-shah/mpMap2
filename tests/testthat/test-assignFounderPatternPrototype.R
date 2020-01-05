context("Test assignFounderPatternPrototype")
test_that("Check that works with a biparental marker, 8-way design",
{
	map <- qtl::sim.map(len = c(100, 100), n.mar = c(20, 20), anchor.tel = TRUE, include.x = FALSE, eq.spacing=TRUE)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)

	newFounders <- rep(c(rep(0L, 4), rep(1L, 4)), times = nMarkers(cross))
	dim(newFounders) <- dim(founders(cross))
	dimnames(newFounders) <- dimnames(founders(cross))

	changed <- cross + mpMap2:::assignFounderPatternPrototype(newFounders)
	changed_cpp <- cross + mpMap2:::assignFounderPattern(newFounders)
	expect_identical(changed, changed_cpp)

	translationVector <- vector(length = 36, mode = "integer")
	translationVector[1:4] <- 0L
	translationVector[5:8] <- 1L
	translationVector[9:36] <- 2L

	manuallyChanged <- cross

	manuallyChangedFounders <- translationVector[manuallyChanged@geneticData[[1]]@founders]
	dim(manuallyChangedFounders) <- dim(founders(manuallyChanged))
	dimnames(manuallyChangedFounders) <- dimnames(founders(manuallyChanged))
	manuallyChanged@geneticData[[1]]@founders <- manuallyChangedFounders

	manuallyChangedFinals <- translationVector[manuallyChanged@geneticData[[1]]@finals]
	dim(manuallyChangedFinals) <- dim(finals(manuallyChanged))
	dimnames(manuallyChangedFinals) <- dimnames(finals(manuallyChanged))
	manuallyChanged@geneticData[[1]]@finals <- manuallyChangedFinals

	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)
})
test_that("Check that works with random SNP markers, 8-way design",
{
	map <- qtl::sim.map(len = c(100, 100), n.mar = c(20, 20), anchor.tel = TRUE, include.x = FALSE, eq.spacing=TRUE)
	pedigree <- eightParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)

	#First with hets
	changed <- cross + multiparentSNP(keepHets = TRUE)
	manuallyChanged <- cross + mpMap2:::assignFounderPatternPrototype(changed@geneticData[[1]]@founders)
	manuallyChanged_cpp <- cross + mpMap2:::assignFounderPattern(changed@geneticData[[1]]@founders)
	expect_identical(manuallyChanged, manuallyChanged_cpp)

	#We need to order the hetData rows for comparisons
	newOrder <- manuallyChanged@geneticData[[1]]@hetData[[1]]
	expect_identical(newOrder[order(newOrder[,3], newOrder[,2], newOrder[,1]),], rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L)))
	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)

	#Without hets
	changed <- cross + multiparentSNP(keepHets = FALSE)
	manuallyChanged <- cross + mpMap2:::assignFounderPatternPrototype(changed@geneticData[[1]]@founders) + removeHets()
	newOrder <- manuallyChanged@geneticData[[1]]@hetData[[1]]
	expect_identical(newOrder[order(newOrder[, 3]), ], rbind(c(0L, 0L, 0L), c(1L, 1L, 1L)))
	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)
})
test_that("Check that works with a biparental marker, 4-way design",
{
	map <- qtl::sim.map(len = c(100, 100), n.mar = c(20, 20), anchor.tel = TRUE, include.x = FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)

	newFounders <- rep(c(rep(0L, 2), rep(1L, 2)), times = nMarkers(cross))
	dim(newFounders) <- dim(founders(cross))
	dimnames(newFounders) <- dimnames(founders(cross))
	changed <- cross + mpMap2:::assignFounderPatternPrototype(newFounders)
	changed_cpp <- cross + mpMap2:::assignFounderPattern(newFounders)
	expect_identical(changed, changed_cpp)

	translationVector <- vector(length = 10, mode = "integer")
	translationVector[1:2] <- 0L
	translationVector[3:4] <- 1L
	translationVector[5:10] <- 2L

	manuallyChanged <- cross

	manuallyChangedFounders <- translationVector[manuallyChanged@geneticData[[1]]@founders]
	dim(manuallyChangedFounders) <- dim(founders(manuallyChanged))
	dimnames(manuallyChangedFounders) <- dimnames(founders(manuallyChanged))
	manuallyChanged@geneticData[[1]]@founders <- manuallyChangedFounders

	manuallyChangedFinals <- translationVector[manuallyChanged@geneticData[[1]]@finals]
	dim(manuallyChangedFinals) <- dim(finals(manuallyChanged))
	dimnames(manuallyChangedFinals) <- dimnames(finals(manuallyChanged))
	manuallyChanged@geneticData[[1]]@finals <- manuallyChangedFinals

	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)
})
test_that("Check that works with random SNP markers, 4-way design",
{
	map <- qtl::sim.map(len = c(100, 100), n.mar = c(20, 20), anchor.tel = TRUE, include.x = FALSE, eq.spacing=TRUE)
	pedigree <- fourParentPedigreeSingleFunnel(initialPopulationSize = 100, selfingGenerations = 1, intercrossingGenerations = 0, nSeeds = 1)
	cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)

	#First with hets
	changed <- cross + multiparentSNP(keepHets = TRUE)
	manuallyChanged <- cross + mpMap2:::assignFounderPatternPrototype(changed@geneticData[[1]]@founders)
	manuallyChanged_cpp <- cross + mpMap2:::assignFounderPattern(changed@geneticData[[1]]@founders)
	expect_identical(manuallyChanged, manuallyChanged_cpp)

	#We need to order the hetData rows for comparisons
	newOrder <- manuallyChanged@geneticData[[1]]@hetData[[1]]
	expect_identical(newOrder[order(newOrder[,3], newOrder[,2], newOrder[,1]),], rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L)))
	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 0L, 2L), c(0L, 1L, 2L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)

	#Without hets
	changed <- cross + multiparentSNP(keepHets = FALSE)
	manuallyChanged <- cross + mpMap2:::assignFounderPatternPrototype(changed@geneticData[[1]]@founders)
	manuallyChanged_cpp <- cross + mpMap2:::assignFounderPattern(changed@geneticData[[1]]@founders)
	expect_identical(manuallyChanged, manuallyChanged_cpp)
	
	manuallyChanged <- manuallyChanged + removeHets()
	manuallyChanged_cpp <- manuallyChanged_cpp + removeHets()
	expect_identical(manuallyChanged, manuallyChanged_cpp)

	newOrder <- manuallyChanged@geneticData[[1]]@hetData[[1]]
	expect_identical(newOrder[order(newOrder[,3], newOrder[,2], newOrder[,1]),], rbind(c(0L, 0L, 0L), c(1L, 1L, 1L)))
	for(i in 1:length(manuallyChanged@geneticData[[1]]@hetData))
	{
		manuallyChanged@geneticData[[1]]@hetData[[i]] <- rbind(c(0L, 0L, 0L), c(1L, 1L, 1L))
		#We need to order the hetData rows for comparisons
		changedPart <- changed@geneticData[[1]]@hetData[[i]]
		changed@geneticData[[1]]@hetData[[i]] <- changedPart[order(changedPart[,3], changedPart[, 2], changedPart[,1]), ]
	}
	expect_identical(changed, manuallyChanged)
})
