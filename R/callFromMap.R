#' @export
#' @title Call markers based on an existing map
#' @description This function uses an existing genetic map to call genetic markers, including markers polymorphic on multiple chromosomes. 
#' @param rawData Raw data for a genetic marker.
#' @param thresholdChromosomes The test-statistic threshold for declaring a marker to be polymorphic on a chromosome. 
#' @param thresholdAlleleClusters The p-value threshold for declaring two underlying founder alleles to have different marker alleles. Multiple possible values should be input. 
#' @param maxChromosomes The maximum number of chromosomes that a marker can be polymorphic on
#' @param existingImputations An object of class mpcrossMapped from the mpMap2 package, containing data about imputed underlying genotypes. 
#' @param tDistributionPValue Paramater controlling the size of each detected cluster, ranging from 0 to 1. Small values result in small clusters, and large values result in large clusters.
#' 
#' @details
#' This function uses an existing genetic map to call a genetic marker. There are a number of advantages to this approach
#' \describe{
#' 	\item{1.}{It can correctly call markers which are polymorphic on multiple chromosomes, therefore converting one marker into two.}
#'	\item{2.}{It avoids incorrectly calling markers polymorphic on multiple chromosomes. Incorrect calling can lead to supurious genetic interactions.}
#' 	\item{3.}{It can call markers that initially appear to be monomorphic in the population. }
#' 	\item{4.}{It can call additional marker alleles for markers that would otherwise be ignored. }
#' }
#' 
#' Once a genetic map has been constructed, it should be used to impute underlying founder genotypes at an equally spaced grid of points using function \code{\link[mpMap2]{imputeFounders}}. 
#' The steps in the algorithm are as follows:
#' \describe{
#' 	\item{1.}{Determine which chromosomes the marker is associated to, and where on those chromosomes. This is determined using function \code{\link[mpMap2]{addExtraMarkerFromRawCall}}, which is itself based on a manova model. 
#'	  The marker is assumed associated to chromosomes for which the test statistic is greater than \code{thresholdChromosomes}. 
#' 	  An appropriate value for \code{thresholdChromosomes} can be determined by looking at the results of \code{\link[mpMap2]{addExtraMarkerFromRawCall}}, for a number of different markers.}
#' 	\item{2.}{Determine the distribution of marker alleles, at all the associated genetic locations. 
#'    This is done by taking the founders to be the vertices of a graph, and connecting founders which seem to part of the same marker allele. 
#' 	  The resulting graph should be a union of disjoint complete graphs (cliques). }
#' 	\item{3.}{We now have a preliminary assignment of marker alleles to lines, where the assignment may be of 1, 2, 3 or more \emph{different} marker alleles, depending on how many chromosomes the marker is associated with.
#' 	  For example, if the marker is associated with two chromosomes, then there will be two marker alleles for each line. 
#' 	  For each unique combination of marker alleles, we take the lines which have that assignment of marker alleles, and fit a skew-t distribution. }
#' 	\item{4.}{For each fitted distribution, determine a confidence region using p-value \code{tDistributionPValue}. }
#' 	\item{5.}{Use these confidence regions to construct marker calls at each associated location.}
#' }
callFromMap <- function(rawData, thresholdChromosomes = 100, thresholdAlleleClusters = c(1e-10, 1e-20, 1e-30, 1e-40), maxChromosomes = 2, existingImputations, tDistributionPValue = 0.6, useOnlyExtraImputationPoints = TRUE)
{
	rawResult <- mpMap2::addExtraMarkerFromRawCall(mpcrossMapped = existingImputations, newMarker = rawData, useOnlyExtraImputationPoints = useOnlyExtraImputationPoints)
	chromosomes <- names(existingImputations@map)
	chromosomeScores <- sapply(chromosomes, function(x) max(rawResult@data[names(rawResult@map[[x]])]))
	chromosomeScores <- sort(chromosomeScores, decreasing=TRUE)
	if(sum(chromosomeScores > thresholdChromosomes) > maxChromosomes) 
	{
		result <- list(called = FALSE)
		attr(result, "failures") <- "Detected on too many chromosomes"
		return(result)
	}
	if(sum(chromosomeScores > thresholdChromosomes) == 0)
	{
		result <- list(called = FALSE)
		attr(result, "failures") <- "Detected on zero chromosomes"
		return(result)
	}

	bestChromosomes <- names(chromosomeScores[chromosomeScores > thresholdChromosomes])
	bestPositionsChromosomes <- sapply(bestChromosomes, function(x) names(which.max(rawResult@data[names(rawResult@map[[x]])])))

	thresholdAlleleClusters <- sort(thresholdAlleleClusters, decreasing = FALSE)
	internalResults <- vector(mode = "character", length = length(thresholdAlleleClusters))
	for(index in 1:length(thresholdAlleleClusters))
	{
		thresholdAlleleCluster <- thresholdAlleleClusters[index]
		result <- callFromMapInternal(bestPositionsChromosomes = bestPositionsChromosomes, rawData = rawData, thresholdAlleleCluster = thresholdAlleleCluster, existingImputations = existingImputations, tDistributionPValue = tDistributionPValue)
		if(!is.null(result) && !is.character(result)) return(result)
		if(is.character(result)) internalResults[index] <- result
	}
	result <- list(called = FALSE)
	attr(result, "failures") <- internalResults
	return(result)
}
callFromMapInternal <- function(bestPositionsChromosomes, rawData, thresholdAlleleCluster, existingImputations, tDistributionPValue)
{
	nFounders <- nFounders(existingImputations)
	dataPerPosition <- clusters <- list()
	pValuesMatrices <- list()
	for(position in bestPositionsChromosomes)
	{
		dataPerPosition[[position]] <- factor(existingImputations@geneticData[[1]]@imputed@data[,position])
		results <- matrix(0, nrow = nFounders, ncol = nFounders)
		for(i in 1:(nFounders-1))
		{
			for(j in (i+1):nFounders)
			{
				contrastRow <- rep(0, nFounders)
				contrastRow[i] <- 1
				contrastRow[j] <- -1
				currentModelData <- dataPerPosition[[position]]
				contrasts(currentModelData) <- cbind(condInterest = contrastRow)
				model <- lm(rawData ~ currentModelData)
				testResult <- car::linearHypothesis(model, "currentModelDatacondInterest")
				
				#Taken from print.linearHypothesis.mlm
				SSPE.qr <- qr(testResult$SSPE)
				eigs <- Re(eigen(qr.coef(SSPE.qr, testResult$SSPH), symmetric = FALSE)$values)
				pillai <- car:::Pillai(eigs, testResult$df, testResult$df.residual)
				ok <- (pillai[2] >= 0) & (pillai[3] > 0) & (pillai[4] > 0)
				ok <- !is.na(ok) & ok
				if(!ok)
				{
					return("Problem computing Pillai-Bartlett test statistic")
				}
				pvalue <- pf(pillai[2], pillai[3], pillai[4], lower.tail = FALSE)

				results[j, i] <- results[i, j] <- pvalue
			}
		}
		pValuesMatrices[[position]] <- results
		adjacencyMatrix <- matrix(0, nrow = nFounders, ncol = nFounders)
		adjacencyMatrix[results > thresholdAlleleCluster] <- 1
		diag(adjacencyMatrix) <- 1
		maxCliques <- igraph::max_cliques(igraph::graph_from_adjacency_matrix(adjacencyMatrix))
		#The max cliques should actually partition the graph, so check that
		maxCliqueVector <- do.call(c, maxCliques)
		if(length(maxCliqueVector) != nFounders || length(unique(maxCliqueVector)) != nFounders) return("Overlapping cliques")
		#If marker is monomorphic, return NULL
		if(length(maxCliques) == 1) return("monoallelic")
		clusters[[position]] <- maxCliques
	}
	classifyPosition <- function(position)
	{
		mapping <- vector(mode = "integer", length = nFounders)
		for(clusterCounter in 1:length(clusters[[position]])) mapping[clusters[[position]][[clusterCounter]]] <- clusterCounter
		return(mapping[as.integer(dataPerPosition[[position]])])
	}
	groups <- sapply(bestPositionsChromosomes, classifyPosition, simplify = FALSE)
	combinedGroups <- do.call(interaction, groups)
	nCombinedGroups <- length(levels(combinedGroups))
	levels(combinedGroups) <- as.character(1:nCombinedGroups)

	overallClusterAssignments <- vector(mode = "integer", length = nrow(rawData))
	names(overallClusterAssignments) <- rownames(rawData)
	overallClusterAssignments[] <- NA
	clusterAssignmentsPerPosition <- list()
	insideNumberOfGroups <- vector(mode = "integer", length = nrow(rawData))
	clusterBoundaries <- list()
	for(group in 1:nCombinedGroups)
	{
		groupData <- rawData[combinedGroups == group,]
		groupData <- as.data.frame(groupData)
		colnames(groupData) <- c("x", "y")
		done <- FALSE
		try(
			{
				setTimeLimit(120, 120, transient = TRUE)
				model <- selm(cbind(x, y) ~ 1, family = "ST", data = groupData, fixed.param = list(alpha = 0))
				distribution <- extractSECdistr(model, compNames = c("x", "y"))
				plotRange <- cbind(range(groupData[,1]), range(groupData[, 2]))
				pdf(NULL)
					plotResults <- plot(distribution, col = 2, probs = tDistributionPValue, landmarks = "", range = plotRange)
				dev.off()
				done <- TRUE
			}, silent = TRUE
		)
		if(!done)
		{
			try(
				{
					setTimeLimit(120, 120, transient = TRUE)
					model <- selm(cbind(x, y) ~ 1, family = "ST", data = groupData)
					distribution <- extractSECdistr(model, compNames = c("x", "y"))
					plotRange <- cbind(range(groupData[,1]), range(groupData[, 2]))
					pdf(NULL)
						plotResults <- plot(distribution, col = 2, probs = tDistributionPValue, landmarks = "", range = plotRange)
					dev.off()
					done <- TRUE
				}, silent = TRUE
			)
		}
		if(!done) 
		{
			setTimeLimit()
			return("Error fitting t distribution")
		}

		clusterBoundaries[[group]] <- plotResults$plot$contourLines[[1]]
		data <- cbind(plotResults$plot$contourLines[[1]]$x, plotResults$plot$contourLines[[1]]$y)
		contourValue <- mean(sn::dmst(data, dp = plotResults$object@dp))

		isInside <- sn::dmst(rawData, dp = plotResults$object@dp) > contourValue

		#Not *exactly* sure why this occassionally comes out as NA. Numerical overflow for the occasional point?
		isInside[is.na(isInside)] <- FALSE

		overallClusterAssignments[isInside] <- group
		insideNumberOfGroups[isInside] <- insideNumberOfGroups[isInside] + 1
	}
	#If it's in multiple groups, then mark the cluster assignment as NA.
	overallClusterAssignments[insideNumberOfGroups > 1] <- NA
	overallClusterAssignments <- factor(overallClusterAssignments, levels = 1:nCombinedGroups)
	classificationsPerPosition <- list()
	for(position in bestPositionsChromosomes)
	{
		mapping <- vector(mode = "integer", length = nFounders)
		for(clusterCounter in 1:length(clusters[[position]])) mapping[clusters[[position]][[clusterCounter]]] <- clusterCounter
		currentPositionClassificationTable <- mapping[dataPerPosition[[position]]]
		conversionTable <- table(overallClusterAssignments, currentPositionClassificationTable)
		#Switch how we assign clusters to marker alleles. Potentially the commented out version is better?
		for(k in 1:ncol(conversionTable)) conversionTable[,k] <- conversionTable[,k] / sum(conversionTable[,k])
		
		result <- vector(mode = "integer", length = nrow(rawData))
		for(group in 1:nCombinedGroups)
		{
			result[overallClusterAssignments == group] <- which.max(conversionTable[group,])
		}
		names(result) <- rownames(rawData)
		classificationsPerPosition[[position]] <- list(finals = result, founders = mapping)
	}
	setTimeLimit()
	return(list(overallAssignment = as.integer(overallClusterAssignments), classificationsPerPosition = classificationsPerPosition, pValuesMatrices = pValuesMatrices, preliminaryGroups = combinedGroups, clusterBoundaries = clusterBoundaries, called = TRUE))
}
