library(mpMap2)
library(purrr)
if(!exists("previousMap"))
{
	previousMapFile <- "/OSM/CBR/AG_MAGICLM/source/data/8way_mapV2/investigate2B2DInfinite.RData"
	load(previousMapFile)
	previousMap <- combined
	rm(combined)
}

#Extract start of chromosomes 1A, 1B and 1D.
toOutput <- subset(previousMap, chromosomes = c("1A", "1B", "1D"))
subMap <- lapply(toOutput@map, function(x) keep(x, function(y) y < 35))
toOutput <- subset(toOutput, markers = unlist(lapply(subMap, names)))

class(subMap) <- "map"
toOutput <- mpcrossMapped(toOutput, map = subMap)

originalImputationMap <- previousMap@geneticData[[1]]@imputed@map
subImputationMap <- lapply(originalImputationMap, function(x) keep(x, function(y) y < 35))
toOutput@geneticData[[1]]@imputed <- subset(previousMap@geneticData[[1]]@imputed, positions = unlist(lapply(subImputationMap, names)))

eightParentSubsetMap <- redact(toOutput)
save(eightParentSubsetMap, file = "../../data/eightParentSubsetMap.RData")
