cmake_minimum_required(VERSION 3.1)
project(mpMap2)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/src)
find_package(R REQUIRED)
find_package(Rcpp REQUIRED)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
add_subdirectory(src)

add_custom_target(copyPackage ALL)	
set(HEADERS alleleDataErrors.h combineGenotypes.h estimateRFCheckFunnels.h estimateRFSpecificDesign.h generateGenotypes.h intercrossingAndSelfingGenerations.h orderFunnel.h recodeHetsAsNA.h checkHets.h crc32.h estimateRF.h funnelsToUniqueValues.h getFunnel.h markerPatternsToUniqueValues.h recodeFoundersFinalsHets.h sortPedigreeLineNames.h unitTypes.hpp fourParentPedigreeRandomFunnels.h matrixChunks.h rawSymmetricMatrix.h dspMatrix.h)
set(RFILES biparentalDominant.R combineGenotypes.R detailedPedigree-class.R estimateRF.R expand.R f2Pedigree.R formGroups.R fourParentPedigreeRandomFunnels.R fourParentPedigreeSingleFunnel.R fullHetData.R geneticData-class.R hetData-class.R lg-class.R map-class.R mapFunctions.R markers.R mpcross-class.R mpcross.R multiparentSNP.R nFounders.R nLines.R nMarkers.R pedigree-class.R pedigree.R pedigreeGraph-class.R pedigreeGraph.R pedigreeToGraph.R print.R Rcpp_exceptions.R removeHets.R rf-class.R rilPedigree.R roxygen.R show.R simulateMPCross.R subset.R twoParentPedigree.R validation.R rawSymmetricMatrix.R)
#Copy package to binary directory. This works differently on windows and linux
if(WIN32)
	if("${CMAKE_GENERATOR}" STREQUAL "NMake Makefiles")
		#Copy DESCRIPTION
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION")
		#Copy include files
		foreach(HEADER IN LISTS HEADERS)
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/src/${HEADER}" "${CMAKE_CURRENT_BINARY_DIR}/inst/include/${HEADER}")
		endforeach()
		#Copy R files
		foreach(RFILE IN LISTS RFILES)
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/R/${RFILE}")
		endforeach()
		#Run roxygen2
		add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e \"library(roxygen2); roxygen2::roxygenise()\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
	else()
		foreach(BUILD_TYPE IN LISTS CMAKE_CONFIGURATION_TYPES)
			#Copy DESCRIPTION
			add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/DESCRIPTION")
			#Copy include files
			foreach(HEADER IN LISTS HEADERS)
				add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} ARGS "-E copy_if_different \"${CMAKE_CURRENT_SOURCE_DIR}/src/${HEADER}\" \"${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/inst/include/${HEADER}\"")
			endforeach()
			#Copy R files
			foreach(RFILE IN LISTS RFILES)
				add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}/R/${RFILE}")
			endforeach()
			#Run roxygen2
			add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e \"library(roxygen2); roxygen2::roxygenise()\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_TYPE}")
		endforeach()
	endif()
else()
	#Copy DESCRIPTION
	add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/DESCRIPTION.in" "${CMAKE_CURRENT_BINARY_DIR}/DESCRIPTION")
	#Copy include files
	foreach(HEADER IN LISTS HEADERS)
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/src/${HEADER}" "${CMAKE_CURRENT_BINARY_DIR}/inst/include/${HEADER}")
	endforeach()
	#Copy R files
	foreach(RFILE IN LISTS RFILES)
		add_custom_command(TARGET copyPackage COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/R/${RFILE}" "${CMAKE_CURRENT_BINARY_DIR}/R/${RFILE}")
	endforeach()
	#Run roxygen2
	add_custom_command(TARGET copyPackage COMMAND "${R_COMMAND}" ARGS "-e" "\"library(roxygen2);roxygen2::roxygenise()\"" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()
if(WIN32)
	install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL --no-multiarch \$\{CMAKE_INSTALL_CONFIG_NAME\} WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")")
else()
	install(CODE "execute_process(COMMAND \"${R_COMMAND}\" CMD INSTALL . WORKING_DIRECTORY \"${PROJECT_BINARY_DIR}\")")
endif()
