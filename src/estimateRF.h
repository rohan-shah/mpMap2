#ifndef _RFHAPS_H
#define _RFHAPS_H
#include <Rcpp.h>
struct funnelType
{
	int val[16];
};
/** Estimate pairwise recombination fractions
  * 
  * Estimate pairwise recombination fractions using numerical maximum likelihood
  * @param object The mpcross object to use
  * @param recombinationFractions The recombination fractions to test. Must be given in increasing order. 
  * @param markerRows The rows of markers (as a pair of indices) which we wish to consider
  * @param markerColumns The columns of markers (as a pair of indices) which we wish to consider
  * @param lineWeights The line weights, in case we wish to correct for some kind of distortion
  * @param keepLod Boolean telling whether or not to return the likelihood ratio statistic for testing the estimated value being different from 0.
  * @param gbLimit The number of gigabytes to use for the results matrix. A value of negative 1 indicates no limit.
  * @param keepLkhd Boolean telling whether or not to return the maximum likelihood value
  * @param verbose Boolean telling whether or not to output diagnostic and progress information
  * @return A list returning the specified data. In the case of theta, the values are returned as a raw vector. Each entry is an index into the possible recombination fractions. This saves us a factor of 8 in terms of memory usage. The raw vector is indexed column-major, but only contains the values for the upper triangular part of the matrix. 
 **/
SEXP estimateRF(SEXP object, SEXP recombinationFractions, SEXP markerRows, SEXP markerColumns, SEXP lineWeights, SEXP keepLod, SEXP keepLkhd, SEXP gbLimit, SEXP verbose);
#endif
