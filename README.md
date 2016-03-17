## mpMap2: Map construction for multi-parent crosses

**Authors**: Rohan Shah, Emma Huang

R/mpMap2 is a software package for constructing genetic linkage maps from a family of experimental designs known as **multi-parent crosses**. These types of designs begin with 2^n inbred founder lines, and may incorporate random intermating and inbreeding by selfing. These designs have found recent application in Arabadopsis, barley, rice, maize, tomatoes and wheat. However there do not appear to be appropriate computational tools to perform map construction using data generated from current experimental populations, which may include thousands of lines and over 100,000 markers.

From a statistical point of view R/mpMap2 is extremely flexible. It can handle 2-, 4-, 8- and 16- parent designs, and allows for the explicit modelling of hetrozygotes in every design. It will also impute underlying hetrozygote genotypes, allowing the identification of residual hetrozygosity for future experiments. 

From a computational point of view, R/mpMap2 aims to allow the user to analyse these types of large populations using the minimum possible computational resources. It outputs diagnostic information about the amount of memory used and the progress of operations, and can be instructed to use a limited amount of working memory for some calculations. Unlike the previous version of this package, R/mpMap2 uses only simple OpenMP multi-threading, making it simpler to compile and run.

R/mpMap2 also has an associated visualisation tool R/mpMapInteractive2, which uses the Qt graphics framework. This tool can be used to visually inspect and alter the data during the mapping process.

This package builds on R/mpMap by Emma Huang. Compared to that previous version, the following significant improvements have been made:
1. R/mpMap2 now uses the S4 object system, and enforces stricter checks on the data. 
2. R/mpMap2 supports designs with finite generations of selfing, so that the genotyped lines do not have to be inbred. 
3. As well as the 4-parent and 8-parent designs, R/mpMap2 now supports biparental designs and 16-parent designs.
4. R/mpMap2 outputs diagnostics about the amount of memory required, and can output the progress of most operations. 
5. When the number of markers and lines is large R/mpMap2 will use significantly less memory and run computations significantly faster. This allows the analysis of data-sets which were previously infeasible. 
6. R/mpMap2 allows the use of multiple data-sets to construct a single genetic map.
7. R/mpMap2 uses only simple OpenMP parallelisation and is therefore simpler to compile and run than R/mpMap, which made use of GPUs and MPI. 

For further details see the package vignette. 

##Package compilation

This package can be compiled in two ways; either using the standard package compilation commands (E.g. R CMD INSTALL) or by using the included CMake build files. The CMake build files allow compilation using Visual Studio on Windows.

###Compilation on Windows using CMake and Visual Studio

1. Download and compile the customized version of Rcpp from github repository rohan-shah/Rcpp. See the associated Readme file for details on compiling that code. 
2. Run the cmake gui. 
3. Set Rcpp_DIR to the Rcpp binaries directory. 
4. Set R_COMMAND to <R_HOME>/bin/x64/R.exe. Ensure that you choose the 64-bit version. 
5. Enter the source directory and the binaries directory (E.g. <mpMap2Root>/build for Visual Studio 64-bit output, or <mpMap2Root>/release for NMake Makefiles output)
6. If the output is going to be NMake Makefiles, set CMAKE_BUILD_TYPE appropriately (E.g. as either Release or Debug)
7. Hit Configure and when prompted choose a Visual Studio 64-bit output, or NMake Makefiles.
8. When configuring succeeds, hit generate. 

The configuration scripts generate an import library for R.dll. This means that the scripts must be able to run cl.exe and lib.exe. If this step fails, check that cl.exe and lib.exe can run. If not, you may need to set up the correct environment for the compiler (by running a script such as vcvarsx86_amd64.bat) before running cmake. 

The package can now be compiled by either running nmake in the binaries directory (NMake Makefiles) or opening mpMap2.sln in the binaries directory. Once the package is compiled a properly formed R package (including NAMESPACE, DESCRIPTION, .R files and C code) will have been constructed in the binaries directory. If Visual Studio output was selected, the package directory will be <mpMap2Binaries>/<buildType> (E.g. mpMap2/build/Release for a release build). If NMake Makefiles output was selected, the package will be <mpMap2Binaries> (E.g mpMap2/build). 

The package can be built using the INSTALL target, which is run using "nmake install" if NMake Makefiles output was selected. Alternatively, after the package is built you can run the standard R CMD INSTALL command. 
