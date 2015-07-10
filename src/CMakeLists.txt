cmake_minimum_required(VERSION 3.1)
project(mpMap2)

#Construct R.lib if we're using Visual studio on Windows
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
	#R does not appear to tell us whether a debug build was selected, so just build a release build every time 
	message(STATUS "Checking that 64-bit R was selected")
	execute_process(COMMAND ${R_COMMAND} --no-restore --no-save -e "q()" OUTPUT_VARIABLE R_OUTPUT)
	string(FIND "${R_OUTPUT}" "64-bit" R_OUTPUT_64_BIT_INDEX)
	string(FIND "${R_OUTPUT}" "32-bit" R_OUTPUT_32_BIT_INDEX)
	if(NOT ${R_OUTPUT_32_BIT_INDEX} EQUAL -1)
		message(FATAL_ERROR "Checking that 64-bit R was selected - 32-bit R found")
	elseif(${R_OUTPUT_64_BIT_INDEX} EQUAL -1)
	 	message(STATUS ${R_OUTPUT})
		message(FATAL_ERROR "Checking that 64-bit R was selected - Unable to identify as either 64 or 32 bit")
	endif()
	message(STATUS "Checking that 64-bit R was selected - passed")

	message(STATUS "Checking that 64-bit compiler was selected")
	execute_process(COMMAND ${CMAKE_CXX_COMPILER} ERROR_VARIABLE COMPILER_OUTPUT_STRING OUTPUT_QUIET)
	string(FIND "${COMPILER_OUTPUT_STRING}" x86 COMPILER_OUTPUT_STRING_X86_INDEX)
	if(NOT ${COMPILER_OUTPUT_STRING_X86_INDEX} EQUAL -1)
		message(FATAL_ERROR "Checking that 64-bit compiler was selected - FAILED\nPlease ensure that you select the x64 developer tools prompt")
	endif()
	string(FIND "${COMPILER_OUTPUT_STRING}" "x64" COMPILER_OUTPUT_STRING_X64_INDEX)
	if(${COMPILER_OUTPUT_STRING_X64_INDEX} EQUAL -1)
		message(FATAL_ERROR "Checking that 64-bit compiler was selected - FAILED\nPlease ensure that you select the x64 developer tools prompt")
	endif()
	message(STATUS "Checking that 64-bit compiler was selected - passed")

	#If we're using MSVC, work out where lib.exe is
	get_filename_component(COMPILER_PATH "${CMAKE_CXX_COMPILER}" DIRECTORY CACHE)
	message(STATUS "Searching for lib.exe")
	find_program(LIB_EXECUTABLE lib PATHS ${COMPILER_PATH})
	if(${LIB_EXECUTABLE} STREQUAL LIB_EXECUTABLE-NOTFOUND)
		message(FATAL_ERROR "Searching for lib.exe - not found")
	else()
		message(STATUS "Searching for lib.exe - found")
	endif()

	#If we're using MSVC, work out where dumpbin.exe is
	message(STATUS "Searching for dumpbin.exe")
	find_program(DUMPBIN_EXECUTABLE dumpbin PATHS ${COMPILER_PATH})
	if(${DUMPBIN_EXECUTABLE} STREQUAL DUMPBIN_EXECUTABLE-NOTFOUND)
		message(FATAL_ERROR "Searching for dumpbin.exe - not found")
	else()
		message(STATUS "Searching for dumpbin.exe - found")
	endif()

	message(STATUS "Generating import library for R.dll")
	file(WRITE "${PROJECT_BINARY_DIR}/R.def" "LIBRARY R\n")
	file(APPEND "${PROJECT_BINARY_DIR}/R.def" "EXPORTS\n")
	execute_process(COMMAND "${DUMPBIN_EXECUTABLE}" /exports "${R_HOME}/bin/x64/R.dll" ERROR_QUIET OUTPUT_VARIABLE DUMPBIN_OUTPUT)
	#Strip off the junk at the start and end
	string(REGEX MATCHALL "[ |\t]+[0-9]+[ |\t]+[0-9A-F]+[ |\t]+0[0-9A-F]+[ |\t]+[^\n|\r]+\n" R_EXPORTS_LIST_UNFILTERED "${DUMPBIN_OUTPUT}")
	foreach(symbolline IN LISTS R_EXPORTS_LIST_UNFILTERED)
		string(REGEX REPLACE "[ |\t]+[0-9]+[ |\t]+[0-9A-F]+[ |\t]+0[0-9A-F]+[ |\t]+([^\n|\r]+)\n" "\t\\1" symbol "${symbolline}")
		file(APPEND "${PROJECT_BINARY_DIR}/R.def" "${symbol}\n")
	endforeach(symbolline)
	execute_process(COMMAND "${LIB_EXECUTABLE}" /MACHINE:X64 /nodefaultlib "/def:${PROJECT_BINARY_DIR}/R.def" "/out:${PROJECT_BINARY_DIR}/R.lib" ERROR_QUIET OUTPUT_QUIET)
	file(REMOVE "${PROJECT_BINARY_DIR}/R.def")
	if(NOT EXISTS "${PROJECT_BINARY_DIR}/R.lib")
		message(FATAL_ERROR "Generating import library for R.dll - Error R.lib not found")
	endif()
	message(STATUS "Generating import library for R.dll - done")

	add_definitions(-Dsnprintf=_snprintf)

	#Put a Makefile.win in the source directory
	find_program(FULL_MAKE_PROGRAM nmake)
	set(FULL_MAKE_PROGRAM "\"${FULL_MAKE_PROGRAM}\" /f ")
	configure_file(${PROJECT_SOURCE_DIR}/Makefile.in ${PROJECT_BINARY_DIR}/Makefile.win @ONLY)
endif()

#Make sure everything installed to the right place
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")

#Now add the shared libarry target
set(SourceFiles alleleDataErrors.cpp checkHets.cpp combineGenotypes.cpp crc32.cpp estimateRF.cpp estimateRFCheckFunnels.cpp estimateRFSpecificDesign.cpp funnelsToUniqueValues.cpp generateGenotypes.cpp getFunnel.cpp intercrossingAndSelfingGenerations.cpp markerPatternsToUniqueValues.cpp orderFunnel.cpp recodeFoundersFinalsHets.cpp register.cpp replaceHetsWithNA.cpp convertGeneticData.cpp)
set(HeaderFiles matrices.hpp constructLookupTable.hpp probabilities.hpp probabilities2.hpp convertGeneticData.h)
add_library(mpMap2 SHARED ${SourceFiles} ${HeaderFiles} ${mpMap2_MOC_SOURCES})
target_link_libraries(mpMap2 PRIVATE Rcpp)
target_include_directories(mpMap2 PUBLIC ${PROJECT_SOURCE_DIR}/../inst/include ${PROJECT_SOURCE_DIR})
target_compile_definitions(mpMap2 PRIVATE CUSTOM_STATIC_RCPP)

#Ensure that it's installed to the right place
if(WIN32)
	target_link_libraries(mpMap2 PRIVATE ${PROJECT_BINARY_DIR}/R.lib)
endif()
