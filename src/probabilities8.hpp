#ifndef PROBABILITIES_8_HEADER_GUARD
#define PROBABILITIES_8_HEADER_GUARD
#include "probabilities.hpp"
#include "probabilities2.hpp"
#include "probabilities4.hpp"
template<> struct probabilityData<8>
{
public:
	/*See Karl Bromans paper on intermediate generations. This mask converts allele encodings (0 - 2) into indices into 
	the array of 4 different probabilities. In terms of indices, 
	Zero = homozygote, One = other homozygote, Two = hetrozygote
	In terms of values, see Table one of the paper. Zero = first equation of table, one = second equation, etc. Note that
	we combine equations 4 and 5 into a single state. 
	*/
	static const int intermediateProbabilitiesMask[64][64];
	/*This mask takes in the two alleles at a *single* location and returns a value encoding that genotype. */
	static const int intermediateAllelesMask[8][8];
	static const int infiniteMask[8][8];
private:
	probabilityData(){}
};
const int probabilityData<8>::intermediateProbabilitiesMask[][64] = 
{{0, 1, 2, 2, 3, 3, 3, 3, 1, 4, 5, 5, 6, 6, 6, 6, 2, 5, 7, 8, 9, 9, 9,
   9, 2, 5, 8, 7, 9, 9, 9, 9, 3, 6, 9, 9, 10, 11, 12, 12, 3, 6, 9, 9, 
  11, 10, 12, 12, 3, 6, 9, 9, 12, 12, 10, 11, 3, 6, 9, 9, 12, 12, 11, 
  10}, {1, 13, 14, 14, 15, 15, 15, 15, 16, 1, 17, 17, 18, 18, 18, 18, 
  17, 14, 8, 19, 20, 20, 20, 20, 17, 14, 19, 8, 20, 20, 20, 20, 18, 
  15, 20, 20, 11, 21, 22, 22, 18, 15, 20, 20, 21, 11, 22, 22, 18, 15, 
  20, 20, 22, 22, 11, 21, 18, 15, 20, 20, 22, 22, 21, 11}, {2, 14, 23,
   24, 25, 25, 25, 25, 17, 5, 24, 26, 27, 27, 27, 27, 28, 29, 2, 17, 
  30, 30, 30, 30, 29, 31, 14, 5, 32, 32, 32, 32, 30, 32, 25, 27, 12, 
  22, 33, 33, 30, 32, 25, 27, 22, 12, 33, 33, 30, 32, 25, 27, 33, 33, 
  12, 22, 30, 32, 25, 27, 33, 33, 22, 12}, {2, 14, 24, 23, 25, 25, 25,
   25, 17, 5, 26, 24, 27, 27, 27, 27, 29, 31, 5, 14, 32, 32, 32, 32, 
  28, 29, 17, 2, 30, 30, 30, 30, 30, 32, 27, 25, 12, 22, 33, 33, 30, 
  32, 27, 25, 22, 12, 33, 33, 30, 32, 27, 25, 33, 33, 12, 22, 30, 32, 
  27, 25, 33, 33, 22, 12}, {3, 15, 25, 25, 34, 35, 36, 36, 18, 6, 27, 
  27, 35, 37, 38, 38, 30, 32, 9, 20, 36, 38, 39, 39, 30, 32, 20, 9, 
  36, 38, 39, 39, 40, 41, 42, 42, 3, 18, 30, 30, 41, 43, 44, 44, 15, 
  6, 32, 32, 42, 44, 45, 45, 25, 27, 9, 20, 42, 44, 45, 45, 25, 27, 
  20, 9}, {3, 15, 25, 25, 35, 34, 36, 36, 18, 6, 27, 27, 37, 35, 38, 
  38, 30, 32, 9, 20, 38, 36, 39, 39, 30, 32, 20, 9, 38, 36, 39, 39, 
  41, 43, 44, 44, 6, 15, 32, 32, 40, 41, 42, 42, 18, 3, 30, 30, 42, 
  44, 45, 45, 27, 25, 9, 20, 42, 44, 45, 45, 27, 25, 20, 9}, {3, 15, 
  25, 25, 36, 36, 34, 35, 18, 6, 27, 27, 38, 38, 35, 37, 30, 32, 9, 
  20, 39, 39, 36, 38, 30, 32, 20, 9, 39, 39, 36, 38, 42, 44, 45, 45, 
  9, 20, 25, 27, 42, 44, 45, 45, 20, 9, 25, 27, 40, 41, 42, 42, 30, 
  30, 3, 18, 41, 43, 44, 44, 32, 32, 15, 6}, {3, 15, 25, 25, 36, 36, 
  35, 34, 18, 6, 27, 27, 38, 38, 37, 35, 30, 32, 9, 20, 39, 39, 38, 
  36, 30, 32, 20, 9, 39, 39, 38, 36, 42, 44, 45, 45, 9, 20, 27, 25, 
  42, 44, 45, 45, 20, 9, 27, 25, 41, 43, 44, 44, 32, 32, 6, 15, 40, 
  41, 42, 42, 30, 30, 18, 3}, {1, 16, 17, 17, 18, 18, 18, 18, 13, 1, 
  14, 14, 15, 15, 15, 15, 14, 17, 8, 19, 20, 20, 20, 20, 14, 17, 19, 
  8, 20, 20, 20, 20, 15, 18, 20, 20, 11, 21, 22, 22, 15, 18, 20, 20, 
  21, 11, 22, 22, 15, 18, 20, 20, 22, 22, 11, 21, 15, 18, 20, 20, 22, 
  22, 21, 11}, {4, 1, 5, 5, 6, 6, 6, 6, 1, 0, 2, 2, 3, 3, 3, 3, 5, 2, 
  7, 8, 9, 9, 9, 9, 5, 2, 8, 7, 9, 9, 9, 9, 6, 3, 9, 9, 10, 11, 12, 
  12, 6, 3, 9, 9, 11, 10, 12, 12, 6, 3, 9, 9, 12, 12, 10, 11, 6, 3, 9,
   9, 12, 12, 11, 10}, {5, 17, 24, 26, 27, 27, 27, 27, 14, 2, 23, 24, 
  25, 25, 25, 25, 29, 28, 2, 17, 30, 30, 30, 30, 31, 29, 14, 5, 32, 
  32, 32, 32, 32, 30, 25, 27, 12, 22, 33, 33, 32, 30, 25, 27, 22, 12, 
  33, 33, 32, 30, 25, 27, 33, 33, 12, 22, 32, 30, 25, 27, 33, 33, 22, 
  12}, {5, 17, 26, 24, 27, 27, 27, 27, 14, 2, 24, 23, 25, 25, 25, 25, 
  31, 29, 5, 14, 32, 32, 32, 32, 29, 28, 17, 2, 30, 30, 30, 30, 32, 
  30, 27, 25, 12, 22, 33, 33, 32, 30, 27, 25, 22, 12, 33, 33, 32, 30, 
  27, 25, 33, 33, 12, 22, 32, 30, 27, 25, 33, 33, 22, 12}, {6, 18, 27,
   27, 35, 37, 38, 38, 15, 3, 25, 25, 34, 35, 36, 36, 32, 30, 9, 20, 
  36, 38, 39, 39, 32, 30, 20, 9, 36, 38, 39, 39, 41, 40, 42, 42, 3, 
  18, 30, 30, 43, 41, 44, 44, 15, 6, 32, 32, 44, 42, 45, 45, 25, 27, 
  9, 20, 44, 42, 45, 45, 25, 27, 20, 9}, {6, 18, 27, 27, 37, 35, 38, 
  38, 15, 3, 25, 25, 35, 34, 36, 36, 32, 30, 9, 20, 38, 36, 39, 39, 
  32, 30, 20, 9, 38, 36, 39, 39, 43, 41, 44, 44, 6, 15, 32, 32, 41, 
  40, 42, 42, 18, 3, 30, 30, 44, 42, 45, 45, 27, 25, 9, 20, 44, 42, 
  45, 45, 27, 25, 20, 9}, {6, 18, 27, 27, 38, 38, 35, 37, 15, 3, 25, 
  25, 36, 36, 34, 35, 32, 30, 9, 20, 39, 39, 36, 38, 32, 30, 20, 9, 
  39, 39, 36, 38, 44, 42, 45, 45, 9, 20, 25, 27, 44, 42, 45, 45, 20, 
  9, 25, 27, 41, 40, 42, 42, 30, 30, 3, 18, 43, 41, 44, 44, 32, 32, 
  15, 6}, {6, 18, 27, 27, 38, 38, 37, 35, 15, 3, 25, 25, 36, 36, 35, 
  34, 32, 30, 9, 20, 39, 39, 38, 36, 32, 30, 20, 9, 39, 39, 38, 36, 
  44, 42, 45, 45, 9, 20, 27, 25, 44, 42, 45, 45, 20, 9, 27, 25, 43, 
  41, 44, 44, 32, 32, 6, 15, 41, 40, 42, 42, 30, 30, 18, 3}, {2, 17, 
  28, 29, 30, 30, 30, 30, 14, 5, 29, 31, 32, 32, 32, 32, 23, 24, 2, 
  14, 25, 25, 25, 25, 24, 26, 17, 5, 27, 27, 27, 27, 25, 27, 30, 32, 
  12, 22, 33, 33, 25, 27, 30, 32, 22, 12, 33, 33, 25, 27, 30, 32, 33, 
  33, 12, 22, 25, 27, 30, 32, 33, 33, 22, 12}, {5, 14, 29, 31, 32, 32,
   32, 32, 17, 2, 28, 29, 30, 30, 30, 30, 24, 23, 2, 14, 25, 25, 25, 
  25, 26, 24, 17, 5, 27, 27, 27, 27, 27, 25, 30, 32, 12, 22, 33, 33, 
  27, 25, 30, 32, 22, 12, 33, 33, 27, 25, 30, 32, 33, 33, 12, 22, 27, 
  25, 30, 32, 33, 33, 22, 12}, {7, 8, 2, 5, 9, 9, 9, 9, 8, 7, 2, 5, 9,
   9, 9, 9, 2, 2, 0, 1, 3, 3, 3, 3, 5, 5, 1, 4, 6, 6, 6, 6, 9, 9, 3, 
  6, 10, 11, 12, 12, 9, 9, 3, 6, 11, 10, 12, 12, 9, 9, 3, 6, 12, 12, 
  10, 11, 9, 9, 3, 6, 12, 12, 11, 10}, {8, 19, 17, 14, 20, 20, 20, 20,
   19, 8, 17, 14, 20, 20, 20, 20, 14, 14, 1, 13, 15, 15, 15, 15, 17, 
  17, 16, 1, 18, 18, 18, 18, 20, 20, 18, 15, 11, 21, 22, 22, 20, 20, 
  18, 15, 21, 11, 22, 22, 20, 20, 18, 15, 22, 22, 11, 21, 20, 20, 18, 
  15, 22, 22, 21, 11}, {9, 20, 30, 32, 36, 38, 39, 39, 20, 9, 30, 32, 
  36, 38, 39, 39, 25, 25, 3, 15, 34, 35, 36, 36, 27, 27, 18, 6, 35, 
  37, 38, 38, 42, 42, 40, 41, 3, 18, 30, 30, 44, 44, 41, 43, 15, 6, 
  32, 32, 45, 45, 42, 44, 25, 27, 9, 20, 45, 45, 42, 44, 25, 27, 20, 
  9}, {9, 20, 30, 32, 38, 36, 39, 39, 20, 9, 30, 32, 38, 36, 39, 39, 
  25, 25, 3, 15, 35, 34, 36, 36, 27, 27, 18, 6, 37, 35, 38, 38, 44, 
  44, 41, 43, 6, 15, 32, 32, 42, 42, 40, 41, 18, 3, 30, 30, 45, 45, 
  42, 44, 27, 25, 9, 20, 45, 45, 42, 44, 27, 25, 20, 9}, {9, 20, 30, 
  32, 39, 39, 36, 38, 20, 9, 30, 32, 39, 39, 36, 38, 25, 25, 3, 15, 
  36, 36, 34, 35, 27, 27, 18, 6, 38, 38, 35, 37, 45, 45, 42, 44, 9, 
  20, 25, 27, 45, 45, 42, 44, 20, 9, 25, 27, 42, 42, 40, 41, 30, 30, 
  3, 18, 44, 44, 41, 43, 32, 32, 15, 6}, {9, 20, 30, 32, 39, 39, 38, 
  36, 20, 9, 30, 32, 39, 39, 38, 36, 25, 25, 3, 15, 36, 36, 35, 34, 
  27, 27, 18, 6, 38, 38, 37, 35, 45, 45, 42, 44, 9, 20, 27, 25, 45, 
  45, 42, 44, 20, 9, 27, 25, 44, 44, 41, 43, 32, 32, 6, 15, 42, 42, 
  40, 41, 30, 30, 18, 3}, {2, 17, 29, 28, 30, 30, 30, 30, 14, 5, 31, 
  29, 32, 32, 32, 32, 24, 26, 5, 17, 27, 27, 27, 27, 23, 24, 14, 2, 
  25, 25, 25, 25, 25, 27, 32, 30, 12, 22, 33, 33, 25, 27, 32, 30, 22, 
  12, 33, 33, 25, 27, 32, 30, 33, 33, 12, 22, 25, 27, 32, 30, 33, 33, 
  22, 12}, {5, 14, 31, 29, 32, 32, 32, 32, 17, 2, 29, 28, 30, 30, 30, 
  30, 26, 24, 5, 17, 27, 27, 27, 27, 24, 23, 14, 2, 25, 25, 25, 25, 
  27, 25, 32, 30, 12, 22, 33, 33, 27, 25, 32, 30, 22, 12, 33, 33, 27, 
  25, 32, 30, 33, 33, 12, 22, 27, 25, 32, 30, 33, 33, 22, 12}, {8, 19,
   14, 17, 20, 20, 20, 20, 19, 8, 14, 17, 20, 20, 20, 20, 17, 17, 1, 
  16, 18, 18, 18, 18, 14, 14, 13, 1, 15, 15, 15, 15, 20, 20, 15, 18, 
  11, 21, 22, 22, 20, 20, 15, 18, 21, 11, 22, 22, 20, 20, 15, 18, 22, 
  22, 11, 21, 20, 20, 15, 18, 22, 22, 21, 11}, {7, 8, 5, 2, 9, 9, 9, 
  9, 8, 7, 5, 2, 9, 9, 9, 9, 5, 5, 4, 1, 6, 6, 6, 6, 2, 2, 1, 0, 3, 3,
   3, 3, 9, 9, 6, 3, 10, 11, 12, 12, 9, 9, 6, 3, 11, 10, 12, 12, 9, 9,
   6, 3, 12, 12, 10, 11, 9, 9, 6, 3, 12, 12, 11, 10}, {9, 20, 32, 30, 
  36, 38, 39, 39, 20, 9, 32, 30, 36, 38, 39, 39, 27, 27, 6, 18, 35, 
  37, 38, 38, 25, 25, 15, 3, 34, 35, 36, 36, 42, 42, 41, 40, 3, 18, 
  30, 30, 44, 44, 43, 41, 15, 6, 32, 32, 45, 45, 44, 42, 25, 27, 9, 
  20, 45, 45, 44, 42, 25, 27, 20, 9}, {9, 20, 32, 30, 38, 36, 39, 39, 
  20, 9, 32, 30, 38, 36, 39, 39, 27, 27, 6, 18, 37, 35, 38, 38, 25, 
  25, 15, 3, 35, 34, 36, 36, 44, 44, 43, 41, 6, 15, 32, 32, 42, 42, 
  41, 40, 18, 3, 30, 30, 45, 45, 44, 42, 27, 25, 9, 20, 45, 45, 44, 
  42, 27, 25, 20, 9}, {9, 20, 32, 30, 39, 39, 36, 38, 20, 9, 32, 30, 
  39, 39, 36, 38, 27, 27, 6, 18, 38, 38, 35, 37, 25, 25, 15, 3, 36, 
  36, 34, 35, 45, 45, 44, 42, 9, 20, 25, 27, 45, 45, 44, 42, 20, 9, 
  25, 27, 42, 42, 41, 40, 30, 30, 3, 18, 44, 44, 43, 41, 32, 32, 15, 
  6}, {9, 20, 32, 30, 39, 39, 38, 36, 20, 9, 32, 30, 39, 39, 38, 36, 
  27, 27, 6, 18, 38, 38, 37, 35, 25, 25, 15, 3, 36, 36, 35, 34, 45, 
  45, 44, 42, 9, 20, 27, 25, 45, 45, 44, 42, 20, 9, 27, 25, 44, 44, 
  43, 41, 32, 32, 6, 15, 42, 42, 41, 40, 30, 30, 18, 3}, {3, 18, 30, 
  30, 40, 41, 42, 42, 15, 6, 32, 32, 41, 43, 44, 44, 25, 27, 9, 20, 
  42, 44, 45, 45, 25, 27, 20, 9, 42, 44, 45, 45, 34, 35, 36, 36, 3, 
  15, 25, 25, 35, 37, 38, 38, 18, 6, 27, 27, 36, 38, 39, 39, 30, 32, 
  9, 20, 36, 38, 39, 39, 30, 32, 20, 9}, {6, 15, 32, 32, 41, 43, 44, 
  44, 18, 3, 30, 30, 40, 41, 42, 42, 27, 25, 9, 20, 42, 44, 45, 45, 
  27, 25, 20, 9, 42, 44, 45, 45, 35, 34, 36, 36, 3, 15, 25, 25, 37, 
  35, 38, 38, 18, 6, 27, 27, 38, 36, 39, 39, 30, 32, 9, 20, 38, 36, 
  39, 39, 30, 32, 20, 9}, {9, 20, 25, 27, 42, 44, 45, 45, 20, 9, 25, 
  27, 42, 44, 45, 45, 30, 30, 3, 18, 40, 41, 42, 42, 32, 32, 15, 6, 
  41, 43, 44, 44, 36, 36, 34, 35, 3, 15, 25, 25, 38, 38, 35, 37, 18, 
  6, 27, 27, 39, 39, 36, 38, 30, 32, 9, 20, 39, 39, 36, 38, 30, 32, 
  20, 9}, {9, 20, 27, 25, 42, 44, 45, 45, 20, 9, 27, 25, 42, 44, 45, 
  45, 32, 32, 6, 15, 41, 43, 44, 44, 30, 30, 18, 3, 40, 41, 42, 42, 
  36, 36, 35, 34, 3, 15, 25, 25, 38, 38, 37, 35, 18, 6, 27, 27, 39, 
  39, 38, 36, 30, 32, 9, 20, 39, 39, 38, 36, 30, 32, 20, 9}, {10, 11, 
  12, 12, 3, 6, 9, 9, 11, 10, 12, 12, 3, 6, 9, 9, 12, 12, 10, 11, 3, 
  6, 9, 9, 12, 12, 11, 10, 3, 6, 9, 9, 3, 3, 3, 3, 0, 1, 2, 2, 6, 6, 
  6, 6, 1, 4, 5, 5, 9, 9, 9, 9, 2, 5, 7, 8, 9, 9, 9, 9, 2, 5, 8, 
  7}, {11, 21, 22, 22, 18, 15, 20, 20, 21, 11, 22, 22, 18, 15, 20, 20,
   22, 22, 11, 21, 18, 15, 20, 20, 22, 22, 21, 11, 18, 15, 20, 20, 15,
   15, 15, 15, 1, 13, 14, 14, 18, 18, 18, 18, 16, 1, 17, 17, 20, 20, 
  20, 20, 17, 14, 8, 19, 20, 20, 20, 20, 17, 14, 19, 8}, {12, 22, 33, 
  33, 30, 32, 25, 27, 22, 12, 33, 33, 30, 32, 25, 27, 33, 33, 12, 22, 
  30, 32, 25, 27, 33, 33, 22, 12, 30, 32, 25, 27, 25, 25, 25, 25, 2, 
  14, 23, 24, 27, 27, 27, 27, 17, 5, 24, 26, 30, 30, 30, 30, 28, 29, 
  2, 17, 32, 32, 32, 32, 29, 31, 14, 5}, {12, 22, 33, 33, 30, 32, 27, 
  25, 22, 12, 33, 33, 30, 32, 27, 25, 33, 33, 12, 22, 30, 32, 27, 25, 
  33, 33, 22, 12, 30, 32, 27, 25, 25, 25, 25, 25, 2, 14, 24, 23, 27, 
  27, 27, 27, 17, 5, 26, 24, 32, 32, 32, 32, 29, 31, 5, 14, 30, 30, 
  30, 30, 28, 29, 17, 2}, {3, 18, 30, 30, 41, 40, 42, 42, 15, 6, 32, 
  32, 43, 41, 44, 44, 25, 27, 9, 20, 44, 42, 45, 45, 25, 27, 20, 9, 
  44, 42, 45, 45, 35, 37, 38, 38, 6, 18, 27, 27, 34, 35, 36, 36, 15, 
  3, 25, 25, 36, 38, 39, 39, 32, 30, 9, 20, 36, 38, 39, 39, 32, 30, 
  20, 9}, {6, 15, 32, 32, 43, 41, 44, 44, 18, 3, 30, 30, 41, 40, 42, 
  42, 27, 25, 9, 20, 44, 42, 45, 45, 27, 25, 20, 9, 44, 42, 45, 45, 
  37, 35, 38, 38, 6, 18, 27, 27, 35, 34, 36, 36, 15, 3, 25, 25, 38, 
  36, 39, 39, 32, 30, 9, 20, 38, 36, 39, 39, 32, 30, 20, 9}, {9, 20, 
  25, 27, 44, 42, 45, 45, 20, 9, 25, 27, 44, 42, 45, 45, 30, 30, 3, 
  18, 41, 40, 42, 42, 32, 32, 15, 6, 43, 41, 44, 44, 38, 38, 35, 37, 
  6, 18, 27, 27, 36, 36, 34, 35, 15, 3, 25, 25, 39, 39, 36, 38, 32, 
  30, 9, 20, 39, 39, 36, 38, 32, 30, 20, 9}, {9, 20, 27, 25, 44, 42, 
  45, 45, 20, 9, 27, 25, 44, 42, 45, 45, 32, 32, 6, 15, 43, 41, 44, 
  44, 30, 30, 18, 3, 41, 40, 42, 42, 38, 38, 37, 35, 6, 18, 27, 27, 
  36, 36, 35, 34, 15, 3, 25, 25, 39, 39, 38, 36, 32, 30, 9, 20, 39, 
  39, 38, 36, 32, 30, 20, 9}, {11, 21, 22, 22, 15, 18, 20, 20, 21, 11,
   22, 22, 15, 18, 20, 20, 22, 22, 11, 21, 15, 18, 20, 20, 22, 22, 21,
   11, 15, 18, 20, 20, 18, 18, 18, 18, 1, 16, 17, 17, 15, 15, 15, 15, 
  13, 1, 14, 14, 20, 20, 20, 20, 14, 17, 8, 19, 20, 20, 20, 20, 14, 
  17, 19, 8}, {10, 11, 12, 12, 6, 3, 9, 9, 11, 10, 12, 12, 6, 3, 9, 9,
   12, 12, 10, 11, 6, 3, 9, 9, 12, 12, 11, 10, 6, 3, 9, 9, 6, 6, 6, 6,
   4, 1, 5, 5, 3, 3, 3, 3, 1, 0, 2, 2, 9, 9, 9, 9, 5, 2, 7, 8, 9, 9, 
  9, 9, 5, 2, 8, 7}, {12, 22, 33, 33, 32, 30, 25, 27, 22, 12, 33, 33, 
  32, 30, 25, 27, 33, 33, 12, 22, 32, 30, 25, 27, 33, 33, 22, 12, 32, 
  30, 25, 27, 27, 27, 27, 27, 5, 17, 24, 26, 25, 25, 25, 25, 14, 2, 
  23, 24, 30, 30, 30, 30, 29, 28, 2, 17, 32, 32, 32, 32, 31, 29, 14, 
  5}, {12, 22, 33, 33, 32, 30, 27, 25, 22, 12, 33, 33, 32, 30, 27, 25,
   33, 33, 12, 22, 32, 30, 27, 25, 33, 33, 22, 12, 32, 30, 27, 25, 27,
   27, 27, 27, 5, 17, 26, 24, 25, 25, 25, 25, 14, 2, 24, 23, 32, 32, 
  32, 32, 31, 29, 5, 14, 30, 30, 30, 30, 29, 28, 17, 2}, {3, 18, 30, 
  30, 42, 42, 40, 41, 15, 6, 32, 32, 44, 44, 41, 43, 25, 27, 9, 20, 
  45, 45, 42, 44, 25, 27, 20, 9, 45, 45, 42, 44, 36, 38, 39, 39, 9, 
  20, 30, 32, 36, 38, 39, 39, 20, 9, 30, 32, 34, 35, 36, 36, 25, 25, 
  3, 15, 35, 37, 38, 38, 27, 27, 18, 6}, {6, 15, 32, 32, 44, 44, 41, 
  43, 18, 3, 30, 30, 42, 42, 40, 41, 27, 25, 9, 20, 45, 45, 42, 44, 
  27, 25, 20, 9, 45, 45, 42, 44, 38, 36, 39, 39, 9, 20, 30, 32, 38, 
  36, 39, 39, 20, 9, 30, 32, 35, 34, 36, 36, 25, 25, 3, 15, 37, 35, 
  38, 38, 27, 27, 18, 6}, {9, 20, 25, 27, 45, 45, 42, 44, 20, 9, 25, 
  27, 45, 45, 42, 44, 30, 30, 3, 18, 42, 42, 40, 41, 32, 32, 15, 6, 
  44, 44, 41, 43, 39, 39, 36, 38, 9, 20, 30, 32, 39, 39, 36, 38, 20, 
  9, 30, 32, 36, 36, 34, 35, 25, 25, 3, 15, 38, 38, 35, 37, 27, 27, 
  18, 6}, {9, 20, 27, 25, 45, 45, 42, 44, 20, 9, 27, 25, 45, 45, 42, 
  44, 32, 32, 6, 15, 44, 44, 41, 43, 30, 30, 18, 3, 42, 42, 40, 41, 
  39, 39, 38, 36, 9, 20, 30, 32, 39, 39, 38, 36, 20, 9, 30, 32, 36, 
  36, 35, 34, 25, 25, 3, 15, 38, 38, 37, 35, 27, 27, 18, 6}, {12, 22, 
  33, 33, 25, 27, 30, 32, 22, 12, 33, 33, 25, 27, 30, 32, 33, 33, 12, 
  22, 25, 27, 30, 32, 33, 33, 22, 12, 25, 27, 30, 32, 30, 30, 30, 30, 
  2, 17, 28, 29, 32, 32, 32, 32, 14, 5, 29, 31, 25, 25, 25, 25, 23, 
  24, 2, 14, 27, 27, 27, 27, 24, 26, 17, 5}, {12, 22, 33, 33, 27, 25, 
  30, 32, 22, 12, 33, 33, 27, 25, 30, 32, 33, 33, 12, 22, 27, 25, 30, 
  32, 33, 33, 22, 12, 27, 25, 30, 32, 32, 32, 32, 32, 5, 14, 29, 31, 
  30, 30, 30, 30, 17, 2, 28, 29, 25, 25, 25, 25, 24, 23, 2, 14, 27, 
  27, 27, 27, 26, 24, 17, 5}, {10, 11, 12, 12, 9, 9, 3, 6, 11, 10, 12,
   12, 9, 9, 3, 6, 12, 12, 10, 11, 9, 9, 3, 6, 12, 12, 11, 10, 9, 9, 
  3, 6, 9, 9, 9, 9, 7, 8, 2, 5, 9, 9, 9, 9, 8, 7, 2, 5, 3, 3, 3, 3, 2,
   2, 0, 1, 6, 6, 6, 6, 5, 5, 1, 4}, {11, 21, 22, 22, 20, 20, 18, 15, 
  21, 11, 22, 22, 20, 20, 18, 15, 22, 22, 11, 21, 20, 20, 18, 15, 22, 
  22, 21, 11, 20, 20, 18, 15, 20, 20, 20, 20, 8, 19, 17, 14, 20, 20, 
  20, 20, 19, 8, 17, 14, 15, 15, 15, 15, 14, 14, 1, 13, 18, 18, 18, 
  18, 17, 17, 16, 1}, {3, 18, 30, 30, 42, 42, 41, 40, 15, 6, 32, 32, 
  44, 44, 43, 41, 25, 27, 9, 20, 45, 45, 44, 42, 25, 27, 20, 9, 45, 
  45, 44, 42, 36, 38, 39, 39, 9, 20, 32, 30, 36, 38, 39, 39, 20, 9, 
  32, 30, 35, 37, 38, 38, 27, 27, 6, 18, 34, 35, 36, 36, 25, 25, 15, 
  3}, {6, 15, 32, 32, 44, 44, 43, 41, 18, 3, 30, 30, 42, 42, 41, 40, 
  27, 25, 9, 20, 45, 45, 44, 42, 27, 25, 20, 9, 45, 45, 44, 42, 38, 
  36, 39, 39, 9, 20, 32, 30, 38, 36, 39, 39, 20, 9, 32, 30, 37, 35, 
  38, 38, 27, 27, 6, 18, 35, 34, 36, 36, 25, 25, 15, 3}, {9, 20, 25, 
  27, 45, 45, 44, 42, 20, 9, 25, 27, 45, 45, 44, 42, 30, 30, 3, 18, 
  42, 42, 41, 40, 32, 32, 15, 6, 44, 44, 43, 41, 39, 39, 36, 38, 9, 
  20, 32, 30, 39, 39, 36, 38, 20, 9, 32, 30, 38, 38, 35, 37, 27, 27, 
  6, 18, 36, 36, 34, 35, 25, 25, 15, 3}, {9, 20, 27, 25, 45, 45, 44, 
  42, 20, 9, 27, 25, 45, 45, 44, 42, 32, 32, 6, 15, 44, 44, 43, 41, 
  30, 30, 18, 3, 42, 42, 41, 40, 39, 39, 38, 36, 9, 20, 32, 30, 39, 
  39, 38, 36, 20, 9, 32, 30, 38, 38, 37, 35, 27, 27, 6, 18, 36, 36, 
  35, 34, 25, 25, 15, 3}, {12, 22, 33, 33, 25, 27, 32, 30, 22, 12, 33,
   33, 25, 27, 32, 30, 33, 33, 12, 22, 25, 27, 32, 30, 33, 33, 22, 12,
   25, 27, 32, 30, 30, 30, 30, 30, 2, 17, 29, 28, 32, 32, 32, 32, 14, 
  5, 31, 29, 27, 27, 27, 27, 24, 26, 5, 17, 25, 25, 25, 25, 23, 24, 
  14, 2}, {12, 22, 33, 33, 27, 25, 32, 30, 22, 12, 33, 33, 27, 25, 32,
   30, 33, 33, 12, 22, 27, 25, 32, 30, 33, 33, 22, 12, 27, 25, 32, 30,
   32, 32, 32, 32, 5, 14, 31, 29, 30, 30, 30, 30, 17, 2, 29, 28, 27, 
  27, 27, 27, 26, 24, 5, 17, 25, 25, 25, 25, 24, 23, 14, 2}, {11, 21, 
  22, 22, 20, 20, 15, 18, 21, 11, 22, 22, 20, 20, 15, 18, 22, 22, 11, 
  21, 20, 20, 15, 18, 22, 22, 21, 11, 20, 20, 15, 18, 20, 20, 20, 20, 
  8, 19, 14, 17, 20, 20, 20, 20, 19, 8, 14, 17, 18, 18, 18, 18, 17, 
  17, 1, 16, 15, 15, 15, 15, 14, 14, 13, 1}, {10, 11, 12, 12, 9, 9, 6,
   3, 11, 10, 12, 12, 9, 9, 6, 3, 12, 12, 10, 11, 9, 9, 6, 3, 12, 12, 
  11, 10, 9, 9, 6, 3, 9, 9, 9, 9, 7, 8, 5, 2, 9, 9, 9, 9, 8, 7, 5, 2, 
  6, 6, 6, 6, 5, 5, 4, 1, 3, 3, 3, 3, 2, 2, 1, 0}};
		/* 0 = AA
		   1 = AB
		   2 = AC
		   3 = AD
		   4 = AE
		   5 = AF
		   6 = AG
		   7 = AH
		   8 = BA
		   9 = BB
		   10 = BC
		   11 = BD
		   12 = BE
		   13 = BF
		   14 = BG
		   15 = BH
		   etc
		   */
const int probabilityData<8>::intermediateAllelesMask[][8] = 
			/*Range[#1, #1 + 7] & /@ Range[0, 63, 8]*/
			{ { 0, 1, 2, 3, 4, 5, 6, 7 }, { 8, 9, 10, 11, 12, 13, 14, 15 }, { 16, 17,
			18, 19, 20, 21, 22, 23 }, { 24, 25, 26, 27, 28, 29, 30, 31 }, { 32, 33,
			34, 35, 36, 37, 38, 39 }, { 40, 41, 42, 43, 44, 45, 46, 47 }, { 48, 49,
			50, 51, 52, 53, 54, 55 }, { 56, 57, 58, 59, 60, 61, 62, 63 } };

const int probabilityData<8>::infiniteMask[][8] = 
		{
			{0, 1, 2, 2, 2, 2, 2, 2},
			{1, 0, 2, 2, 2, 2, 2, 2},
			{2, 2, 0, 1, 2, 2, 2, 2},
			{2, 2, 1, 0, 2, 2, 2, 2},
			{2, 2, 2, 2, 0, 1, 2, 2},
			{2, 2, 2, 2, 1, 0, 2, 2},
			{2, 2, 2, 2, 2, 2, 0, 1},
			{2, 2, 2, 2, 2, 2, 1, 0}
		};
template<> void genotypeProbabilitiesNoIntercross<8, true>(double(&prob)[nDifferentProbs], double r, int, std::size_t)
{
	prob[0] = (1-r)*(1-r)/(8*(1 + 2*r));
	prob[1] = (1-r)*r/(8+16*r);
	prob[2] = r/(16+32*r);
}
template<> void genotypeProbabilitiesNoIntercross<8, false>(double(&prob)[nDifferentProbs], double r, int selfingGenerations, std::size_t nFunnels)
{
	double onePlus2R = 1 + 2 *r;
	double rSquared = r*r;
	double oneMinusRSquared = (1-r)*(1-r);
	double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
	double powD1 = std::pow(1 + 2*(-1 + r)*r, selfingGenerations);
	double oneMinusRPow4 = oneMinusRSquared*oneMinusRSquared;
	double oneMinusRCubed = (1-r)*oneMinusRSquared;
  	double oneMinusR = 1 - r;
  	double pow2 = std::pow(2, selfingGenerations);
    
	prob[0] = (oneMinusRSquared*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(16*onePlus2R*pow2);
	prob[1] = 0;
	prob[2] = 0;
	prob[3] = -(oneMinusRSquared*(-1 + powD1))/(64*pow2);
	prob[4] = (oneMinusR*r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(16*onePlus2R*pow2);
	prob[5] = 0;
	prob[6] = -(oneMinusR*(-1 + powD1)*r)/(64*pow2);
	prob[7] = (r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(32*onePlus2R*pow2);
	prob[8] = 0;
	prob[9] = -((-1 + powD1)*r)/(128*pow2);
	prob[10] = (-2 + onePlus2R*powD1 + powOneMinus2R - 4*r + 4*pow2*r - 2*powOneMinus2R*r)/(64*onePlus2R*pow2);
	prob[11] = 0;
	prob[12] = 0;
	prob[13] = 0;
	prob[14] = 0;
	prob[15] = 0;
	prob[16] = 0;
	prob[17] = 0;
	prob[18] = 0;
	prob[19] = 0;
	prob[20] = 0;
	prob[21] = 0;
	prob[22] = 0;
	prob[23] = 0;
	prob[24] = 0;
	prob[25] = 0;
	prob[26] = 0;
	prob[27] = 0;
	prob[28] = 0;
	prob[29] = 0;
	prob[30] = 0;
	prob[31] = 0;
	prob[32] = 0;
	prob[33] = 0;
	prob[34] = (oneMinusRPow4*(powD1 + powOneMinus2R))/(64*pow2);
	prob[35] = (oneMinusRCubed*(powD1 + powOneMinus2R)*r)/(64*pow2);
	prob[36] = (oneMinusRSquared*(powD1 + powOneMinus2R)*r)/(128*pow2);
	prob[37] = (oneMinusRSquared*(powD1 + powOneMinus2R)*rSquared)/(64*pow2);
	prob[38] = (oneMinusR*(powD1 + powOneMinus2R)*rSquared)/(128*pow2);
	prob[39] = ((powD1 + powOneMinus2R)*rSquared)/(256*pow2);
	prob[40] = -(oneMinusRPow4*(-powD1 + powOneMinus2R))/(64*pow2);
	prob[41] = -(oneMinusRCubed*(-powD1 + powOneMinus2R)*r)/(64*pow2);
	prob[42] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*r)/(128*pow2);
	prob[43] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*rSquared)/(64*pow2);
	prob[44] = -(oneMinusR*(-powD1 + powOneMinus2R)*rSquared)/(128*pow2);
	prob[45] = ((powD1 - powOneMinus2R)*rSquared)/(256*pow2);
#ifndef NDEBUG
	        double sum = 0;
	        for(int i = 0; i < 18; i++) sum += prob[i];
#endif
}
template<> void genotypeProbabilitiesWithIntercross<8, true>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)*(1-r)/8 + (2 * r + 1 - powOneMinusR)/64)/(1 + 2 * r);
	prob[1] = prob[2] = (1 - 8*prob[0])/56;
}
template<> void genotypeProbabilitiesWithIntercross<8, false>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		throw std::runtime_error("The case of a single funnel with intercrossing is not yet handled");
	}
	double oneMinusR = (1 - r);
	double powOneMinusR1 = std::pow(oneMinusR, 1+nAIGenerations);
	double powOneMinusR2 = std::pow(oneMinusR, nAIGenerations - 1);
	double powOneMinusR3 = std::pow(oneMinusR, 2*nAIGenerations);
	double onePlus2R = (1 + 2*r);
	double oneMinusRSquared = oneMinusR*oneMinusR; 
	double oneMinus2R = (1 - 2*r);
	double oneMinus2RSquared = oneMinus2R*oneMinus2R;
	double pow2 = std::pow(2, selfingGenerations);
	double powOneMinus2R = std::pow(oneMinus2R, selfingGenerations);
	double quadratic1 = 7 + 2*r*(-5 + 2*r);
	double toPowD1 = std::pow(1 - 2*oneMinusR*r, selfingGenerations);
	double quadratic1Squared = quadratic1*quadratic1;
	double oneMinusRCubed = oneMinusR*oneMinusR*oneMinusR;
	double quadratic2 = 1.0/64.0 + (-(1.0/64.0) + oneMinusRCubed/8)* powOneMinusR2;
	quadratic2 *= quadratic2;
	double quadratic3 = 1 - 8 *(1.0/64.0 + (-(1.0/64.0) + oneMinusRCubed/8)*powOneMinusR2);
	quadratic3 *= quadratic3;
	double quadratic4 = 7 - oneMinus2R*powOneMinusR2*quadratic1*(1 - r) - 7*r;
	quadratic4 *= quadratic4;

	prob[0] = (-16*oneMinusR*(-(oneMinusR*onePlus2R*(-7 + 4*pow2)) - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1*(-3 + 4*pow2 - 2*oneMinus2R*powOneMinus2R - 6*r)) + onePlus2R*(49*oneMinusRSquared + 18*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(512*oneMinusRSquared*onePlus2R*pow2);
	prob[1] = (56*oneMinusRSquared + 24*oneMinus2R*powOneMinusR1*quadratic1 - (49*oneMinusRSquared + 18*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(896*oneMinusRSquared*pow2);
	prob[2] = (56*oneMinusRSquared + 24*oneMinus2R*powOneMinusR1*quadratic1 - (49*oneMinusRSquared + 18*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(448*oneMinusRSquared*pow2);
	prob[3] = (56*oneMinusRSquared + 24*oneMinus2R*powOneMinusR1*quadratic1 - (49*oneMinusRSquared + 18*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(224*oneMinusRSquared*pow2);
	prob[4] = (112*(7*oneMinusRSquared*onePlus2R*(-7 + 4*pow2) - oneMinus2R*powOneMinusR1*quadratic1*(-3 + 4*pow2 - 2*oneMinus2R*powOneMinus2R - 6*r)) + onePlus2R*(2401*oneMinusRSquared - 126*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(25088*oneMinusRSquared*onePlus2R*pow2);
	prob[5] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[6] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[7] = (112*(7*oneMinusRSquared*onePlus2R*(-7 + 4*pow2) - oneMinus2R*powOneMinusR1*quadratic1*(-3 + 4*pow2 - 2*oneMinus2R*powOneMinus2R - 6*r)) + onePlus2R*(2401*oneMinusRSquared - 126*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(12544*oneMinusRSquared*onePlus2R*pow2);
	prob[8] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(6272*oneMinusRSquared*pow2);
	prob[9] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(784*oneMinusRSquared*pow2);
	prob[10] = (112*(7*oneMinusRSquared*onePlus2R*(-7 + 4*pow2) - oneMinus2R*powOneMinusR1*quadratic1*(-3 + 4*pow2 - 2*oneMinus2R*powOneMinus2R - 6*r)) + onePlus2R*(2401*oneMinusRSquared - 126*oneMinus2R*powOneMinusR1*quadratic1 + oneMinus2RSquared*powOneMinusR3*quadratic1Squared)*toPowD1)/(6272*oneMinusRSquared*onePlus2R*pow2);
	prob[11] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[12] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(56*oneMinusR + (-49*oneMinusR - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[13] = (quadratic3*(-powOneMinus2R + toPowD1) + 3136*quadratic2*(powOneMinus2R + toPowD1))/(784*pow2);
	prob[14] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[15] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[16] = (3136*quadratic2*(-powOneMinus2R + toPowD1) + quadratic3*(powOneMinus2R + toPowD1))/(784*pow2);
	prob[17] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[18] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[19] = (quadratic4*toPowD1)/(12544*oneMinusRSquared*pow2);
	prob[20] = (quadratic4*toPowD1)/(784*oneMinusRSquared*pow2);
	prob[21] = (quadratic4*toPowD1)/(6272*oneMinusRSquared*pow2);
	prob[22] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
	prob[23] = (quadratic3*(-powOneMinus2R + toPowD1) + 3136*quadratic2*(powOneMinus2R + toPowD1))/(392*pow2);
	prob[24] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(6272*oneMinusRSquared*pow2);
	prob[25] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(784*oneMinusRSquared*pow2);
	prob[26] = (quadratic4*toPowD1)/(12544*oneMinusRSquared*pow2);
	prob[27] = (quadratic4*toPowD1)/(784*oneMinusRSquared*pow2);
	prob[28] = (3136*quadratic2*(-powOneMinus2R + toPowD1) + quadratic3*(powOneMinus2R + toPowD1))/(392*pow2);
	prob[29] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(6272*oneMinusRSquared*pow2);
	prob[30] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(784*oneMinusRSquared*pow2);
	prob[31] = (quadratic4*toPowD1)/(12544*oneMinusRSquared*pow2);
	prob[32] = (quadratic4*toPowD1)/(784*oneMinusRSquared*pow2);
	prob[33] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
	prob[34] = (quadratic3*(-powOneMinus2R + toPowD1) + 3136*quadratic2*(powOneMinus2R + toPowD1))/(196*pow2);
	prob[35] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[36] = -((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (-7*oneMinusR - 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[37] = (quadratic4*toPowD1)/(6272*oneMinusRSquared*pow2);
	prob[38] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
	prob[39] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
	prob[40] = (3136*quadratic2*(-powOneMinus2R + toPowD1) + quadratic3*(powOneMinus2R + toPowD1))/(196*pow2);
	prob[41] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(3136*oneMinusRSquared*pow2);
	prob[42] = ((7 - oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*(-4*oneMinus2R*oneMinusR*powOneMinus2R*powOneMinusR2*quadratic1 + (7 + 3*oneMinus2R*oneMinusR*powOneMinusR2*quadratic1 - 7*r)*toPowD1))/(1568*oneMinusRSquared*pow2);
	prob[43] = (quadratic4*toPowD1)/(6272*oneMinusRSquared*pow2);
	prob[44] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
	prob[45] = (quadratic4*toPowD1)/(1568*oneMinusRSquared*pow2);
#ifndef NDEBUG
	double sum = 0;
	for(int i = 0; i < 18; i++) sum += prob[i];
#endif
	//This is because we combined some states (see mathematica code)
	prob[0] /= 8;
	prob[1] /= 32;
	prob[2] /= 64;
	prob[3] /= 128;
	prob[4] /= 8;
	prob[5] /= 64;
	prob[6] /= 128;
	prob[7] /= 16;
	prob[8] /= 32;
	prob[9] /= 256;
	prob[10] /= 32;
	prob[11] /= 64;
	prob[12] /= 128;
	prob[13] /= 8;
	prob[14] /= 64;
	prob[15] /= 128;
	prob[16] /= 8;
	prob[17] /= 64;
	prob[18] /= 128;
	prob[19] /= 16;
	prob[20] /= 256;
	prob[21] /= 32;
	prob[22] /= 128;
	prob[23] /= 16;
	prob[24] /= 32;
	prob[25] /= 256;
	prob[26] /= 16;
	prob[27] /= 256;
	prob[28] /= 16;
	prob[29] /= 32;
	prob[30] /= 256;
	prob[31] /= 16;
	prob[32] /= 256;
	prob[33] /= 128;
	prob[34] /= 32;
	prob[35] /= 64;
	prob[36] /= 128;
	prob[37] /= 32;
	prob[38] /= 128;
	prob[39] /= 128;
	prob[40] /= 32;
	prob[41] /= 64;
	prob[42] /= 128;
	prob[43] /= 32;
	prob[44] /= 128;
	prob[45] /= 128;
}
#endif