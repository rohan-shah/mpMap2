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
		{{1, 2, 3, 3, 4, 4, 4, 4, 2, 5, 6, 6, 7, 7, 7, 7, 3, 6, 8, 9, 10, 10, 
  10, 10, 3, 6, 9, 8, 10, 10, 10, 10, 4, 7, 10, 10, 11, 12, 13, 13, 4,
   7, 10, 10, 12, 11, 13, 13, 4, 7, 10, 10, 13, 13, 11, 12, 4, 7, 10, 
  10, 13, 13, 12, 11}, {2, 14, 15, 15, 16, 16, 16, 16, 17, 2, 18, 18, 
  19, 19, 19, 19, 18, 15, 9, 20, 21, 21, 21, 21, 18, 15, 20, 9, 21, 
  21, 21, 21, 19, 16, 21, 21, 12, 22, 23, 23, 19, 16, 21, 21, 22, 12, 
  23, 23, 19, 16, 21, 21, 23, 23, 12, 22, 19, 16, 21, 21, 23, 23, 22, 
  12}, {3, 15, 24, 25, 26, 26, 26, 26, 18, 6, 25, 27, 28, 28, 28, 28, 
  29, 30, 3, 18, 31, 31, 31, 31, 30, 32, 15, 6, 33, 33, 33, 33, 31, 
  33, 26, 28, 13, 23, 34, 34, 31, 33, 26, 28, 23, 13, 34, 34, 31, 33, 
  26, 28, 34, 34, 13, 23, 31, 33, 26, 28, 34, 34, 23, 13}, {3, 15, 25,
   24, 26, 26, 26, 26, 18, 6, 27, 25, 28, 28, 28, 28, 30, 32, 6, 15, 
  33, 33, 33, 33, 29, 30, 18, 3, 31, 31, 31, 31, 31, 33, 28, 26, 13, 
  23, 34, 34, 31, 33, 28, 26, 23, 13, 34, 34, 31, 33, 28, 26, 34, 34, 
  13, 23, 31, 33, 28, 26, 34, 34, 23, 13}, {4, 16, 26, 26, 35, 36, 37,
   37, 19, 7, 28, 28, 36, 38, 39, 39, 31, 33, 10, 21, 37, 39, 40, 40, 
  31, 33, 21, 10, 37, 39, 40, 40, 41, 42, 43, 43, 4, 19, 31, 31, 42, 
  44, 45, 45, 16, 7, 33, 33, 43, 45, 46, 46, 26, 28, 10, 21, 43, 45, 
  46, 46, 26, 28, 21, 10}, {4, 16, 26, 26, 36, 35, 37, 37, 19, 7, 28, 
  28, 38, 36, 39, 39, 31, 33, 10, 21, 39, 37, 40, 40, 31, 33, 21, 10, 
  39, 37, 40, 40, 42, 44, 45, 45, 7, 16, 33, 33, 41, 42, 43, 43, 19, 
  4, 31, 31, 43, 45, 46, 46, 28, 26, 10, 21, 43, 45, 46, 46, 28, 26, 
  21, 10}, {4, 16, 26, 26, 37, 37, 35, 36, 19, 7, 28, 28, 39, 39, 36, 
  38, 31, 33, 10, 21, 40, 40, 37, 39, 31, 33, 21, 10, 40, 40, 37, 39, 
  43, 45, 46, 46, 10, 21, 26, 28, 43, 45, 46, 46, 21, 10, 26, 28, 41, 
  42, 43, 43, 31, 31, 4, 19, 42, 44, 45, 45, 33, 33, 16, 7}, {4, 16, 
  26, 26, 37, 37, 36, 35, 19, 7, 28, 28, 39, 39, 38, 36, 31, 33, 10, 
  21, 40, 40, 39, 37, 31, 33, 21, 10, 40, 40, 39, 37, 43, 45, 46, 46, 
  10, 21, 28, 26, 43, 45, 46, 46, 21, 10, 28, 26, 42, 44, 45, 45, 33, 
  33, 7, 16, 41, 42, 43, 43, 31, 31, 19, 4}, {2, 17, 18, 18, 19, 19, 
  19, 19, 14, 2, 15, 15, 16, 16, 16, 16, 15, 18, 9, 20, 21, 21, 21, 
  21, 15, 18, 20, 9, 21, 21, 21, 21, 16, 19, 21, 21, 12, 22, 23, 23, 
  16, 19, 21, 21, 22, 12, 23, 23, 16, 19, 21, 21, 23, 23, 12, 22, 16, 
  19, 21, 21, 23, 23, 22, 12}, {5, 2, 6, 6, 7, 7, 7, 7, 2, 1, 3, 3, 4,
   4, 4, 4, 6, 3, 8, 9, 10, 10, 10, 10, 6, 3, 9, 8, 10, 10, 10, 10, 7,
   4, 10, 10, 11, 12, 13, 13, 7, 4, 10, 10, 12, 11, 13, 13, 7, 4, 10, 
  10, 13, 13, 11, 12, 7, 4, 10, 10, 13, 13, 12, 11}, {6, 18, 25, 27, 
  28, 28, 28, 28, 15, 3, 24, 25, 26, 26, 26, 26, 30, 29, 3, 18, 31, 
  31, 31, 31, 32, 30, 15, 6, 33, 33, 33, 33, 33, 31, 26, 28, 13, 23, 
  34, 34, 33, 31, 26, 28, 23, 13, 34, 34, 33, 31, 26, 28, 34, 34, 13, 
  23, 33, 31, 26, 28, 34, 34, 23, 13}, {6, 18, 27, 25, 28, 28, 28, 28,
   15, 3, 25, 24, 26, 26, 26, 26, 32, 30, 6, 15, 33, 33, 33, 33, 30, 
  29, 18, 3, 31, 31, 31, 31, 33, 31, 28, 26, 13, 23, 34, 34, 33, 31, 
  28, 26, 23, 13, 34, 34, 33, 31, 28, 26, 34, 34, 13, 23, 33, 31, 28, 
  26, 34, 34, 23, 13}, {7, 19, 28, 28, 36, 38, 39, 39, 16, 4, 26, 26, 
  35, 36, 37, 37, 33, 31, 10, 21, 37, 39, 40, 40, 33, 31, 21, 10, 37, 
  39, 40, 40, 42, 41, 43, 43, 4, 19, 31, 31, 44, 42, 45, 45, 16, 7, 
  33, 33, 45, 43, 46, 46, 26, 28, 10, 21, 45, 43, 46, 46, 26, 28, 21, 
  10}, {7, 19, 28, 28, 38, 36, 39, 39, 16, 4, 26, 26, 36, 35, 37, 37, 
  33, 31, 10, 21, 39, 37, 40, 40, 33, 31, 21, 10, 39, 37, 40, 40, 44, 
  42, 45, 45, 7, 16, 33, 33, 42, 41, 43, 43, 19, 4, 31, 31, 45, 43, 
  46, 46, 28, 26, 10, 21, 45, 43, 46, 46, 28, 26, 21, 10}, {7, 19, 28,
   28, 39, 39, 36, 38, 16, 4, 26, 26, 37, 37, 35, 36, 33, 31, 10, 21, 
  40, 40, 37, 39, 33, 31, 21, 10, 40, 40, 37, 39, 45, 43, 46, 46, 10, 
  21, 26, 28, 45, 43, 46, 46, 21, 10, 26, 28, 42, 41, 43, 43, 31, 31, 
  4, 19, 44, 42, 45, 45, 33, 33, 16, 7}, {7, 19, 28, 28, 39, 39, 38, 
  36, 16, 4, 26, 26, 37, 37, 36, 35, 33, 31, 10, 21, 40, 40, 39, 37, 
  33, 31, 21, 10, 40, 40, 39, 37, 45, 43, 46, 46, 10, 21, 28, 26, 45, 
  43, 46, 46, 21, 10, 28, 26, 44, 42, 45, 45, 33, 33, 7, 16, 42, 41, 
  43, 43, 31, 31, 19, 4}, {3, 18, 29, 30, 31, 31, 31, 31, 15, 6, 30, 
  32, 33, 33, 33, 33, 24, 25, 3, 15, 26, 26, 26, 26, 25, 27, 18, 6, 
  28, 28, 28, 28, 26, 28, 31, 33, 13, 23, 34, 34, 26, 28, 31, 33, 23, 
  13, 34, 34, 26, 28, 31, 33, 34, 34, 13, 23, 26, 28, 31, 33, 34, 34, 
  23, 13}, {6, 15, 30, 32, 33, 33, 33, 33, 18, 3, 29, 30, 31, 31, 31, 
  31, 25, 24, 3, 15, 26, 26, 26, 26, 27, 25, 18, 6, 28, 28, 28, 28, 
  28, 26, 31, 33, 13, 23, 34, 34, 28, 26, 31, 33, 23, 13, 34, 34, 28, 
  26, 31, 33, 34, 34, 13, 23, 28, 26, 31, 33, 34, 34, 23, 13}, {8, 9, 
  3, 6, 10, 10, 10, 10, 9, 8, 3, 6, 10, 10, 10, 10, 3, 3, 1, 2, 4, 4, 
  4, 4, 6, 6, 2, 5, 7, 7, 7, 7, 10, 10, 4, 7, 11, 12, 13, 13, 10, 10, 
  4, 7, 12, 11, 13, 13, 10, 10, 4, 7, 13, 13, 11, 12, 10, 10, 4, 7, 
  13, 13, 12, 11}, {9, 20, 18, 15, 21, 21, 21, 21, 20, 9, 18, 15, 21, 
  21, 21, 21, 15, 15, 2, 14, 16, 16, 16, 16, 18, 18, 17, 2, 19, 19, 
  19, 19, 21, 21, 19, 16, 12, 22, 23, 23, 21, 21, 19, 16, 22, 12, 23, 
  23, 21, 21, 19, 16, 23, 23, 12, 22, 21, 21, 19, 16, 23, 23, 22, 
  12}, {10, 21, 31, 33, 37, 39, 40, 40, 21, 10, 31, 33, 37, 39, 40, 
  40, 26, 26, 4, 16, 35, 36, 37, 37, 28, 28, 19, 7, 36, 38, 39, 39, 
  43, 43, 41, 42, 4, 19, 31, 31, 45, 45, 42, 44, 16, 7, 33, 33, 46, 
  46, 43, 45, 26, 28, 10, 21, 46, 46, 43, 45, 26, 28, 21, 10}, {10, 
  21, 31, 33, 39, 37, 40, 40, 21, 10, 31, 33, 39, 37, 40, 40, 26, 26, 
  4, 16, 36, 35, 37, 37, 28, 28, 19, 7, 38, 36, 39, 39, 45, 45, 42, 
  44, 7, 16, 33, 33, 43, 43, 41, 42, 19, 4, 31, 31, 46, 46, 43, 45, 
  28, 26, 10, 21, 46, 46, 43, 45, 28, 26, 21, 10}, {10, 21, 31, 33, 
  40, 40, 37, 39, 21, 10, 31, 33, 40, 40, 37, 39, 26, 26, 4, 16, 37, 
  37, 35, 36, 28, 28, 19, 7, 39, 39, 36, 38, 46, 46, 43, 45, 10, 21, 
  26, 28, 46, 46, 43, 45, 21, 10, 26, 28, 43, 43, 41, 42, 31, 31, 4, 
  19, 45, 45, 42, 44, 33, 33, 16, 7}, {10, 21, 31, 33, 40, 40, 39, 37,
   21, 10, 31, 33, 40, 40, 39, 37, 26, 26, 4, 16, 37, 37, 36, 35, 28, 
  28, 19, 7, 39, 39, 38, 36, 46, 46, 43, 45, 10, 21, 28, 26, 46, 46, 
  43, 45, 21, 10, 28, 26, 45, 45, 42, 44, 33, 33, 7, 16, 43, 43, 41, 
  42, 31, 31, 19, 4}, {3, 18, 30, 29, 31, 31, 31, 31, 15, 6, 32, 30, 
  33, 33, 33, 33, 25, 27, 6, 18, 28, 28, 28, 28, 24, 25, 15, 3, 26, 
  26, 26, 26, 26, 28, 33, 31, 13, 23, 34, 34, 26, 28, 33, 31, 23, 13, 
  34, 34, 26, 28, 33, 31, 34, 34, 13, 23, 26, 28, 33, 31, 34, 34, 23, 
  13}, {6, 15, 32, 30, 33, 33, 33, 33, 18, 3, 30, 29, 31, 31, 31, 31, 
  27, 25, 6, 18, 28, 28, 28, 28, 25, 24, 15, 3, 26, 26, 26, 26, 28, 
  26, 33, 31, 13, 23, 34, 34, 28, 26, 33, 31, 23, 13, 34, 34, 28, 26, 
  33, 31, 34, 34, 13, 23, 28, 26, 33, 31, 34, 34, 23, 13}, {9, 20, 15,
   18, 21, 21, 21, 21, 20, 9, 15, 18, 21, 21, 21, 21, 18, 18, 2, 17, 
  19, 19, 19, 19, 15, 15, 14, 2, 16, 16, 16, 16, 21, 21, 16, 19, 12, 
  22, 23, 23, 21, 21, 16, 19, 22, 12, 23, 23, 21, 21, 16, 19, 23, 23, 
  12, 22, 21, 21, 16, 19, 23, 23, 22, 12}, {8, 9, 6, 3, 10, 10, 10, 
  10, 9, 8, 6, 3, 10, 10, 10, 10, 6, 6, 5, 2, 7, 7, 7, 7, 3, 3, 2, 1, 
  4, 4, 4, 4, 10, 10, 7, 4, 11, 12, 13, 13, 10, 10, 7, 4, 12, 11, 13, 
  13, 10, 10, 7, 4, 13, 13, 11, 12, 10, 10, 7, 4, 13, 13, 12, 
  11}, {10, 21, 33, 31, 37, 39, 40, 40, 21, 10, 33, 31, 37, 39, 40, 
  40, 28, 28, 7, 19, 36, 38, 39, 39, 26, 26, 16, 4, 35, 36, 37, 37, 
  43, 43, 42, 41, 4, 19, 31, 31, 45, 45, 44, 42, 16, 7, 33, 33, 46, 
  46, 45, 43, 26, 28, 10, 21, 46, 46, 45, 43, 26, 28, 21, 10}, {10, 
  21, 33, 31, 39, 37, 40, 40, 21, 10, 33, 31, 39, 37, 40, 40, 28, 28, 
  7, 19, 38, 36, 39, 39, 26, 26, 16, 4, 36, 35, 37, 37, 45, 45, 44, 
  42, 7, 16, 33, 33, 43, 43, 42, 41, 19, 4, 31, 31, 46, 46, 45, 43, 
  28, 26, 10, 21, 46, 46, 45, 43, 28, 26, 21, 10}, {10, 21, 33, 31, 
  40, 40, 37, 39, 21, 10, 33, 31, 40, 40, 37, 39, 28, 28, 7, 19, 39, 
  39, 36, 38, 26, 26, 16, 4, 37, 37, 35, 36, 46, 46, 45, 43, 10, 21, 
  26, 28, 46, 46, 45, 43, 21, 10, 26, 28, 43, 43, 42, 41, 31, 31, 4, 
  19, 45, 45, 44, 42, 33, 33, 16, 7}, {10, 21, 33, 31, 40, 40, 39, 37,
   21, 10, 33, 31, 40, 40, 39, 37, 28, 28, 7, 19, 39, 39, 38, 36, 26, 
  26, 16, 4, 37, 37, 36, 35, 46, 46, 45, 43, 10, 21, 28, 26, 46, 46, 
  45, 43, 21, 10, 28, 26, 45, 45, 44, 42, 33, 33, 7, 16, 43, 43, 42, 
  41, 31, 31, 19, 4}, {4, 19, 31, 31, 41, 42, 43, 43, 16, 7, 33, 33, 
  42, 44, 45, 45, 26, 28, 10, 21, 43, 45, 46, 46, 26, 28, 21, 10, 43, 
  45, 46, 46, 35, 36, 37, 37, 4, 16, 26, 26, 36, 38, 39, 39, 19, 7, 
  28, 28, 37, 39, 40, 40, 31, 33, 10, 21, 37, 39, 40, 40, 31, 33, 21, 
  10}, {7, 16, 33, 33, 42, 44, 45, 45, 19, 4, 31, 31, 41, 42, 43, 43, 
  28, 26, 10, 21, 43, 45, 46, 46, 28, 26, 21, 10, 43, 45, 46, 46, 36, 
  35, 37, 37, 4, 16, 26, 26, 38, 36, 39, 39, 19, 7, 28, 28, 39, 37, 
  40, 40, 31, 33, 10, 21, 39, 37, 40, 40, 31, 33, 21, 10}, {10, 21, 
  26, 28, 43, 45, 46, 46, 21, 10, 26, 28, 43, 45, 46, 46, 31, 31, 4, 
  19, 41, 42, 43, 43, 33, 33, 16, 7, 42, 44, 45, 45, 37, 37, 35, 36, 
  4, 16, 26, 26, 39, 39, 36, 38, 19, 7, 28, 28, 40, 40, 37, 39, 31, 
  33, 10, 21, 40, 40, 37, 39, 31, 33, 21, 10}, {10, 21, 28, 26, 43, 
  45, 46, 46, 21, 10, 28, 26, 43, 45, 46, 46, 33, 33, 7, 16, 42, 44, 
  45, 45, 31, 31, 19, 4, 41, 42, 43, 43, 37, 37, 36, 35, 4, 16, 26, 
  26, 39, 39, 38, 36, 19, 7, 28, 28, 40, 40, 39, 37, 31, 33, 10, 21, 
  40, 40, 39, 37, 31, 33, 21, 10}, {11, 12, 13, 13, 4, 7, 10, 10, 12, 
  11, 13, 13, 4, 7, 10, 10, 13, 13, 11, 12, 4, 7, 10, 10, 13, 13, 12, 
  11, 4, 7, 10, 10, 4, 4, 4, 4, 1, 2, 3, 3, 7, 7, 7, 7, 2, 5, 6, 6, 
  10, 10, 10, 10, 3, 6, 8, 9, 10, 10, 10, 10, 3, 6, 9, 8}, {12, 22, 
  23, 23, 19, 16, 21, 21, 22, 12, 23, 23, 19, 16, 21, 21, 23, 23, 12, 
  22, 19, 16, 21, 21, 23, 23, 22, 12, 19, 16, 21, 21, 16, 16, 16, 16, 
  2, 14, 15, 15, 19, 19, 19, 19, 17, 2, 18, 18, 21, 21, 21, 21, 18, 
  15, 9, 20, 21, 21, 21, 21, 18, 15, 20, 9}, {13, 23, 34, 34, 31, 33, 
  26, 28, 23, 13, 34, 34, 31, 33, 26, 28, 34, 34, 13, 23, 31, 33, 26, 
  28, 34, 34, 23, 13, 31, 33, 26, 28, 26, 26, 26, 26, 3, 15, 24, 25, 
  28, 28, 28, 28, 18, 6, 25, 27, 31, 31, 31, 31, 29, 30, 3, 18, 33, 
  33, 33, 33, 30, 32, 15, 6}, {13, 23, 34, 34, 31, 33, 28, 26, 23, 13,
   34, 34, 31, 33, 28, 26, 34, 34, 13, 23, 31, 33, 28, 26, 34, 34, 23,
   13, 31, 33, 28, 26, 26, 26, 26, 26, 3, 15, 25, 24, 28, 28, 28, 28, 
  18, 6, 27, 25, 33, 33, 33, 33, 30, 32, 6, 15, 31, 31, 31, 31, 29, 
  30, 18, 3}, {4, 19, 31, 31, 42, 41, 43, 43, 16, 7, 33, 33, 44, 42, 
  45, 45, 26, 28, 10, 21, 45, 43, 46, 46, 26, 28, 21, 10, 45, 43, 46, 
  46, 36, 38, 39, 39, 7, 19, 28, 28, 35, 36, 37, 37, 16, 4, 26, 26, 
  37, 39, 40, 40, 33, 31, 10, 21, 37, 39, 40, 40, 33, 31, 21, 10}, {7,
   16, 33, 33, 44, 42, 45, 45, 19, 4, 31, 31, 42, 41, 43, 43, 28, 26, 
  10, 21, 45, 43, 46, 46, 28, 26, 21, 10, 45, 43, 46, 46, 38, 36, 39, 
  39, 7, 19, 28, 28, 36, 35, 37, 37, 16, 4, 26, 26, 39, 37, 40, 40, 
  33, 31, 10, 21, 39, 37, 40, 40, 33, 31, 21, 10}, {10, 21, 26, 28, 
  45, 43, 46, 46, 21, 10, 26, 28, 45, 43, 46, 46, 31, 31, 4, 19, 42, 
  41, 43, 43, 33, 33, 16, 7, 44, 42, 45, 45, 39, 39, 36, 38, 7, 19, 
  28, 28, 37, 37, 35, 36, 16, 4, 26, 26, 40, 40, 37, 39, 33, 31, 10, 
  21, 40, 40, 37, 39, 33, 31, 21, 10}, {10, 21, 28, 26, 45, 43, 46, 
  46, 21, 10, 28, 26, 45, 43, 46, 46, 33, 33, 7, 16, 44, 42, 45, 45, 
  31, 31, 19, 4, 42, 41, 43, 43, 39, 39, 38, 36, 7, 19, 28, 28, 37, 
  37, 36, 35, 16, 4, 26, 26, 40, 40, 39, 37, 33, 31, 10, 21, 40, 40, 
  39, 37, 33, 31, 21, 10}, {12, 22, 23, 23, 16, 19, 21, 21, 22, 12, 
  23, 23, 16, 19, 21, 21, 23, 23, 12, 22, 16, 19, 21, 21, 23, 23, 22, 
  12, 16, 19, 21, 21, 19, 19, 19, 19, 2, 17, 18, 18, 16, 16, 16, 16, 
  14, 2, 15, 15, 21, 21, 21, 21, 15, 18, 9, 20, 21, 21, 21, 21, 15, 
  18, 20, 9}, {11, 12, 13, 13, 7, 4, 10, 10, 12, 11, 13, 13, 7, 4, 10,
   10, 13, 13, 11, 12, 7, 4, 10, 10, 13, 13, 12, 11, 7, 4, 10, 10, 7, 
  7, 7, 7, 5, 2, 6, 6, 4, 4, 4, 4, 2, 1, 3, 3, 10, 10, 10, 10, 6, 3, 
  8, 9, 10, 10, 10, 10, 6, 3, 9, 8}, {13, 23, 34, 34, 33, 31, 26, 28, 
  23, 13, 34, 34, 33, 31, 26, 28, 34, 34, 13, 23, 33, 31, 26, 28, 34, 
  34, 23, 13, 33, 31, 26, 28, 28, 28, 28, 28, 6, 18, 25, 27, 26, 26, 
  26, 26, 15, 3, 24, 25, 31, 31, 31, 31, 30, 29, 3, 18, 33, 33, 33, 
  33, 32, 30, 15, 6}, {13, 23, 34, 34, 33, 31, 28, 26, 23, 13, 34, 34,
   33, 31, 28, 26, 34, 34, 13, 23, 33, 31, 28, 26, 34, 34, 23, 13, 33,
   31, 28, 26, 28, 28, 28, 28, 6, 18, 27, 25, 26, 26, 26, 26, 15, 3, 
  25, 24, 33, 33, 33, 33, 32, 30, 6, 15, 31, 31, 31, 31, 30, 29, 18, 
  3}, {4, 19, 31, 31, 43, 43, 41, 42, 16, 7, 33, 33, 45, 45, 42, 44, 
  26, 28, 10, 21, 46, 46, 43, 45, 26, 28, 21, 10, 46, 46, 43, 45, 37, 
  39, 40, 40, 10, 21, 31, 33, 37, 39, 40, 40, 21, 10, 31, 33, 35, 36, 
  37, 37, 26, 26, 4, 16, 36, 38, 39, 39, 28, 28, 19, 7}, {7, 16, 33, 
  33, 45, 45, 42, 44, 19, 4, 31, 31, 43, 43, 41, 42, 28, 26, 10, 21, 
  46, 46, 43, 45, 28, 26, 21, 10, 46, 46, 43, 45, 39, 37, 40, 40, 10, 
  21, 31, 33, 39, 37, 40, 40, 21, 10, 31, 33, 36, 35, 37, 37, 26, 26, 
  4, 16, 38, 36, 39, 39, 28, 28, 19, 7}, {10, 21, 26, 28, 46, 46, 43, 
  45, 21, 10, 26, 28, 46, 46, 43, 45, 31, 31, 4, 19, 43, 43, 41, 42, 
  33, 33, 16, 7, 45, 45, 42, 44, 40, 40, 37, 39, 10, 21, 31, 33, 40, 
  40, 37, 39, 21, 10, 31, 33, 37, 37, 35, 36, 26, 26, 4, 16, 39, 39, 
  36, 38, 28, 28, 19, 7}, {10, 21, 28, 26, 46, 46, 43, 45, 21, 10, 28,
   26, 46, 46, 43, 45, 33, 33, 7, 16, 45, 45, 42, 44, 31, 31, 19, 4, 
  43, 43, 41, 42, 40, 40, 39, 37, 10, 21, 31, 33, 40, 40, 39, 37, 21, 
  10, 31, 33, 37, 37, 36, 35, 26, 26, 4, 16, 39, 39, 38, 36, 28, 28, 
  19, 7}, {13, 23, 34, 34, 26, 28, 31, 33, 23, 13, 34, 34, 26, 28, 31,
   33, 34, 34, 13, 23, 26, 28, 31, 33, 34, 34, 23, 13, 26, 28, 31, 33,
   31, 31, 31, 31, 3, 18, 29, 30, 33, 33, 33, 33, 15, 6, 30, 32, 26, 
  26, 26, 26, 24, 25, 3, 15, 28, 28, 28, 28, 25, 27, 18, 6}, {13, 23, 
  34, 34, 28, 26, 31, 33, 23, 13, 34, 34, 28, 26, 31, 33, 34, 34, 13, 
  23, 28, 26, 31, 33, 34, 34, 23, 13, 28, 26, 31, 33, 33, 33, 33, 33, 
  6, 15, 30, 32, 31, 31, 31, 31, 18, 3, 29, 30, 26, 26, 26, 26, 25, 
  24, 3, 15, 28, 28, 28, 28, 27, 25, 18, 6}, {11, 12, 13, 13, 10, 10, 
  4, 7, 12, 11, 13, 13, 10, 10, 4, 7, 13, 13, 11, 12, 10, 10, 4, 7, 
  13, 13, 12, 11, 10, 10, 4, 7, 10, 10, 10, 10, 8, 9, 3, 6, 10, 10, 
  10, 10, 9, 8, 3, 6, 4, 4, 4, 4, 3, 3, 1, 2, 7, 7, 7, 7, 6, 6, 2, 
  5}, {12, 22, 23, 23, 21, 21, 19, 16, 22, 12, 23, 23, 21, 21, 19, 16,
   23, 23, 12, 22, 21, 21, 19, 16, 23, 23, 22, 12, 21, 21, 19, 16, 21,
   21, 21, 21, 9, 20, 18, 15, 21, 21, 21, 21, 20, 9, 18, 15, 16, 16, 
  16, 16, 15, 15, 2, 14, 19, 19, 19, 19, 18, 18, 17, 2}, {4, 19, 31, 
  31, 43, 43, 42, 41, 16, 7, 33, 33, 45, 45, 44, 42, 26, 28, 10, 21, 
  46, 46, 45, 43, 26, 28, 21, 10, 46, 46, 45, 43, 37, 39, 40, 40, 10, 
  21, 33, 31, 37, 39, 40, 40, 21, 10, 33, 31, 36, 38, 39, 39, 28, 28, 
  7, 19, 35, 36, 37, 37, 26, 26, 16, 4}, {7, 16, 33, 33, 45, 45, 44, 
  42, 19, 4, 31, 31, 43, 43, 42, 41, 28, 26, 10, 21, 46, 46, 45, 43, 
  28, 26, 21, 10, 46, 46, 45, 43, 39, 37, 40, 40, 10, 21, 33, 31, 39, 
  37, 40, 40, 21, 10, 33, 31, 38, 36, 39, 39, 28, 28, 7, 19, 36, 35, 
  37, 37, 26, 26, 16, 4}, {10, 21, 26, 28, 46, 46, 45, 43, 21, 10, 26,
   28, 46, 46, 45, 43, 31, 31, 4, 19, 43, 43, 42, 41, 33, 33, 16, 7, 
  45, 45, 44, 42, 40, 40, 37, 39, 10, 21, 33, 31, 40, 40, 37, 39, 21, 
  10, 33, 31, 39, 39, 36, 38, 28, 28, 7, 19, 37, 37, 35, 36, 26, 26, 
  16, 4}, {10, 21, 28, 26, 46, 46, 45, 43, 21, 10, 28, 26, 46, 46, 45,
   43, 33, 33, 7, 16, 45, 45, 44, 42, 31, 31, 19, 4, 43, 43, 42, 41, 
  40, 40, 39, 37, 10, 21, 33, 31, 40, 40, 39, 37, 21, 10, 33, 31, 39, 
  39, 38, 36, 28, 28, 7, 19, 37, 37, 36, 35, 26, 26, 16, 4}, {13, 23, 
  34, 34, 26, 28, 33, 31, 23, 13, 34, 34, 26, 28, 33, 31, 34, 34, 13, 
  23, 26, 28, 33, 31, 34, 34, 23, 13, 26, 28, 33, 31, 31, 31, 31, 31, 
  3, 18, 30, 29, 33, 33, 33, 33, 15, 6, 32, 30, 28, 28, 28, 28, 25, 
  27, 6, 18, 26, 26, 26, 26, 24, 25, 15, 3}, {13, 23, 34, 34, 28, 26, 
  33, 31, 23, 13, 34, 34, 28, 26, 33, 31, 34, 34, 13, 23, 28, 26, 33, 
  31, 34, 34, 23, 13, 28, 26, 33, 31, 33, 33, 33, 33, 6, 15, 32, 30, 
  31, 31, 31, 31, 18, 3, 30, 29, 28, 28, 28, 28, 27, 25, 6, 18, 26, 
  26, 26, 26, 25, 24, 15, 3}, {12, 22, 23, 23, 21, 21, 16, 19, 22, 12,
   23, 23, 21, 21, 16, 19, 23, 23, 12, 22, 21, 21, 16, 19, 23, 23, 22,
   12, 21, 21, 16, 19, 21, 21, 21, 21, 9, 20, 15, 18, 21, 21, 21, 21, 
  20, 9, 15, 18, 19, 19, 19, 19, 18, 18, 2, 17, 16, 16, 16, 16, 15, 
  15, 14, 2}, {11, 12, 13, 13, 10, 10, 7, 4, 12, 11, 13, 13, 10, 10, 
  7, 4, 13, 13, 11, 12, 10, 10, 7, 4, 13, 13, 12, 11, 10, 10, 7, 4, 
  10, 10, 10, 10, 8, 9, 6, 3, 10, 10, 10, 10, 9, 8, 6, 3, 7, 7, 7, 7, 
  6, 6, 5, 2, 4, 4, 4, 4, 3, 3, 2, 1}};
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
}
template<> void genotypeProbabilitiesWithIntercross<8, true>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)/4 + (2 * r + 1 - powOneMinusR)/16)/(1 + 2 * r);
	prob[1] = prob[2] = (1 - 4*prob[0])/12;
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
