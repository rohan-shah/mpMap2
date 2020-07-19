cmake_minimum_required(VERSION 3.1)
project(mpMap2)

option(USE_OPENMP "Use OpenMP" OFF)
option(USE_BOOST "Use Boost" OFF)
option(INTERNAL_CHECKS "Run internal consistency checks" OFF)

if(INTERNAL_CHECKS)
	add_definitions(-DINTERNAL_CHECKS)
endif()

if(USE_BOOST)
	#Boost stuff. Only used for the boost graph library. 
	set(Boost_USE_STATIC_LIBS       OFF)
	set(Boost_USE_MULTITHREADED     ON)
	set(Boost_USE_STATIC_RUNTIME    OFF)
	#Find boost packages
	find_package(Boost 1.50.0 REQUIRED COMPONENTS graph regex)
	add_definitions(-DBOOST_NO_AUTO_PTR)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/src)
if(APPLE)
	find_library(R_LIBRARY R REQUIRED)
else()
	find_package(R REQUIRED)
endif()
find_package(Rcpp REQUIRED)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
add_subdirectory(src)

add_custom_target(copyPackage ALL)	
set(RFILES biparentalDominant.R combineGenotypes.R detailedPedigree-class.R estimateRF.R expand.R f2Pedigree.R formGroups.R fourParentPedigreeRandomFunnels.R fourParentPedigreeSingleFunnel.R fullHetData.R geneticData-class.R hetData-class.R lg-class.R map-class.R mapFunctions.R markers.R mpcross-class.R mpcross.R multiparentSNP.R multiparentSNPPrototype.R nFounders.R nLines.R nMarkers.R pedigree-class.R pedigree.R pedigreeGraph-class.R pedigreeGraph.R pedigreeToGraph.R print.R Rcpp_exceptions.R removeHets.R rf-class.R rilPedigree.R roxygen.R show.R simulateMPCross.R subset.R twoParentPedigree.R validation.R rawSymmetricMatrix.R orderCross.R eightWayPedigreeRandomFunnels.R impute.R sixteenParentPedigreeRandomFunnels.R eightWayPedigreeSingleFunnel.R imputeFounders.R estimateMap.R jitterMap.R founders.R finals.R hetData.R fixedNumberOfFounderAlleles.R compressedProbabilities.R backcrossPedigree.R eightWayPedigreeImproperFunnels.R reorderPedigree.R testDistortion.R lineNames.R selfing.R as.mpInterval.R computeGenotypeProbabilities.R purdyToPedigree.R getAllFunnels.R simulatePhenotypes.R listCodingErrors.R getChromosomes.R getPositions.R changeMarkerPosition.R toMpMap.R stripPedigree.R getIntercrossingAndSelfingGenerations.R transposeProbabilities.R plotProbabilities.R reverseChromosomes.R plot.R plotMosaic.R exportMapToPretzl.R createSNPTemplate.R extraImputationPoints.R addExtraMarkers.R assignFounderPatternPrototype.R assignFounderPattern.R combineKeepRF.R imputationGenerics.R estimateRFSingleDesign.R expandedProbabilities.R estimateMapFromImputation.R canSkipValidity.R additionOperators.R computeAllEpistaticChiSquared.R generateIntervalMidPoints.R generateGridPositions.R addExtraMarkerFromRawCall.R singleLocusProbabilities.R redact.R callFromMap.R probabilityData.R)
#Copy package to binary directory. This works differently on windows and linux
if(WIN32)
	if("${CMAKE_GENERATOR}" STREQUAL "NMake Makefiles")
		#Copy DESCRIPTION
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION")
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/NAMESPACE" "${CMAKE_CURRENT_BINARY_DIR}/NAMESPACE")
		#Copy R files
		foreach(RFILE IN LISTS RFILES)
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/R/${RFILE}")
		endforeach()
		#Run roxygen2
		add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e \"library(roxygen2);roxygen2::roxygenise(clean=TRUE)\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
	else()
		message(STATUS "R_COMMAND = ${R_COMMAND}")
		foreach(BUILD_TYPE IN LISTS CMAKE_CONFIGURATION_TYPES)
			#Copy DESCRIPTION
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/DESCRIPTION")
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/NAMESPACE" "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/NAMESPACE")
			#Copy R files
			foreach(RFILE IN LISTS RFILES)
				add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/R/${RFILE}")
			endforeach()
			#Run roxygen2
			add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e" "\"library(roxygen2);roxygen2::roxygenise(clean=TRUE)\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}")
		endforeach()
	endif()
else()
	#Copy DESCRIPTION
	add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION")
	add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/NAMESPACE" "${CMAKE_CURRENT_BINARY_DIR}/NAMESPACE")
	#Copy R files
	foreach(RFILE IN LISTS RFILES)
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/R/${RFILE}")
	endforeach()
	#Run roxygen2
	add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e" "\"library(roxygen2);roxygen2::roxygenise(clean=TRUE)\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()
if(WIN32)
	install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL --no-multiarch \$\{CMAKE_INSTALL_CONFIG_NAME\} WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")")
else()
	install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL . WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")")
endif()
