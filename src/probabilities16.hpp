#ifndef PROBABILITIES_16_HEADER_GUARD
#define PROBABILITIES_16_HEADER_GUARD
#include "probabilities.hpp"
#include "probabilities2.hpp"
#include "probabilities4.hpp"
template<> struct probabilityData<16>
{
public:
	/*See Karl Bromans paper on intermediate generations. This mask converts allele encodings (0 - 2) into indices into 
	the array of 4 different probabilities. In terms of indices, 
	Zero = homozygote, One = other homozygote, Two = hetrozygote
	In terms of values, see Table one of the paper. Zero = first equation of table, one = second equation, etc. Note that
	we combine equations 4 and 5 into a single state. 
	*/
	static const int intermediateProbabilitiesMask[256][256];
	/*This mask takes in the two alleles at a *single* location and returns a value encoding that genotype. */
	static const int intermediateAllelesMask[16][16];
	static const int infiniteMask[16][16];
private:
	probabilityData(){}
};
const int probabilityData<16>::intermediateProbabilitiesMask[][256] = 
{{0, 1, 2, 2, 3, 3, 3, 3, 1, 5, 6, 6, 7, 7, 7, 7, 2, 6, 9, 10, 11, 11,
   11, 11, 2, 6, 10, 9, 11, 11, 11, 11, 3, 7, 11, 11, 13, 14, 15, 15, 
  3, 7, 11, 11, 14, 13, 15, 15, 3, 7, 11, 11, 15, 15, 13, 14, 3, 7, 
  11, 11, 15, 15, 14, 13}, {1, 21, 22, 22, 23, 23, 23, 23, 25, 1, 26, 
  26, 27, 27, 27, 27, 26, 22, 10, 29, 30, 30, 30, 30, 26, 22, 29, 10, 
  30, 30, 30, 30, 27, 23, 30, 30, 14, 32, 33, 33, 27, 23, 30, 30, 32, 
  14, 33, 33, 27, 23, 30, 30, 33, 33, 14, 32, 27, 23, 30, 30, 33, 33, 
  32, 14}, {2, 22, 38, 39, 40, 40, 40, 40, 26, 6, 39, 42, 43, 43, 43, 
  43, 45, 46, 2, 26, 47, 47, 47, 47, 46, 49, 22, 6, 50, 50, 50, 50, 
  47, 50, 40, 43, 15, 33, 52, 52, 47, 50, 40, 43, 33, 15, 52, 52, 47, 
  50, 40, 43, 52, 52, 15, 33, 47, 50, 40, 43, 52, 52, 33, 15}, {2, 22,
   39, 38, 40, 40, 40, 40, 26, 6, 42, 39, 43, 43, 43, 43, 46, 49, 6, 
  22, 50, 50, 50, 50, 45, 46, 26, 2, 47, 47, 47, 47, 47, 50, 43, 40, 
  15, 33, 52, 52, 47, 50, 43, 40, 33, 15, 52, 52, 47, 50, 43, 40, 52, 
  52, 15, 33, 47, 50, 43, 40, 52, 52, 33, 15}, {3, 23, 40, 40, 56, 57,
   58, 58, 27, 7, 43, 43, 57, 60, 61, 61, 47, 50, 11, 30, 58, 61, 63, 
  63, 47, 50, 30, 11, 58, 61, 63, 63, 65, 66, 67, 67, 3, 27, 47, 47, 
  66, 69, 70, 70, 23, 7, 50, 50, 67, 70, 72, 72, 40, 43, 11, 30, 67, 
  70, 72, 72, 40, 43, 30, 11}, {3, 23, 40, 40, 57, 56, 58, 58, 27, 7, 
  43, 43, 60, 57, 61, 61, 47, 50, 11, 30, 61, 58, 63, 63, 47, 50, 30, 
  11, 61, 58, 63, 63, 66, 69, 70, 70, 7, 23, 50, 50, 65, 66, 67, 67, 
  27, 3, 47, 47, 67, 70, 72, 72, 43, 40, 11, 30, 67, 70, 72, 72, 43, 
  40, 30, 11}, {3, 23, 40, 40, 58, 58, 56, 57, 27, 7, 43, 43, 61, 61, 
  57, 60, 47, 50, 11, 30, 63, 63, 58, 61, 47, 50, 30, 11, 63, 63, 58, 
  61, 67, 70, 72, 72, 11, 30, 40, 43, 67, 70, 72, 72, 30, 11, 40, 43, 
  65, 66, 67, 67, 47, 47, 3, 27, 66, 69, 70, 70, 50, 50, 23, 7}, {3, 
  23, 40, 40, 58, 58, 57, 56, 27, 7, 43, 43, 61, 61, 60, 57, 47, 50, 
  11, 30, 63, 63, 61, 58, 47, 50, 30, 11, 63, 63, 61, 58, 67, 70, 72, 
  72, 11, 30, 43, 40, 67, 70, 72, 72, 30, 11, 43, 40, 66, 69, 70, 70, 
  50, 50, 7, 23, 65, 66, 67, 67, 47, 47, 27, 3}, {1, 25, 26, 26, 27, 
  27, 27, 27, 21, 1, 22, 22, 23, 23, 23, 23, 22, 26, 10, 29, 30, 30, 
  30, 30, 22, 26, 29, 10, 30, 30, 30, 30, 23, 27, 30, 30, 14, 32, 33, 
  33, 23, 27, 30, 30, 32, 14, 33, 33, 23, 27, 30, 30, 33, 33, 14, 32, 
  23, 27, 30, 30, 33, 33, 32, 14}, {5, 1, 6, 6, 7, 7, 7, 7, 1, 0, 2, 
  2, 3, 3, 3, 3, 6, 2, 9, 10, 11, 11, 11, 11, 6, 2, 10, 9, 11, 11, 11,
   11, 7, 3, 11, 11, 13, 14, 15, 15, 7, 3, 11, 11, 14, 13, 15, 15, 7, 
  3, 11, 11, 15, 15, 13, 14, 7, 3, 11, 11, 15, 15, 14, 13}, {6, 26, 
  39, 42, 43, 43, 43, 43, 22, 2, 38, 39, 40, 40, 40, 40, 46, 45, 2, 
  26, 47, 47, 47, 47, 49, 46, 22, 6, 50, 50, 50, 50, 50, 47, 40, 43, 
  15, 33, 52, 52, 50, 47, 40, 43, 33, 15, 52, 52, 50, 47, 40, 43, 52, 
  52, 15, 33, 50, 47, 40, 43, 52, 52, 33, 15}, {6, 26, 42, 39, 43, 43,
   43, 43, 22, 2, 39, 38, 40, 40, 40, 40, 49, 46, 6, 22, 50, 50, 50, 
  50, 46, 45, 26, 2, 47, 47, 47, 47, 50, 47, 43, 40, 15, 33, 52, 52, 
  50, 47, 43, 40, 33, 15, 52, 52, 50, 47, 43, 40, 52, 52, 15, 33, 50, 
  47, 43, 40, 52, 52, 33, 15}, {7, 27, 43, 43, 57, 60, 61, 61, 23, 3, 
  40, 40, 56, 57, 58, 58, 50, 47, 11, 30, 58, 61, 63, 63, 50, 47, 30, 
  11, 58, 61, 63, 63, 66, 65, 67, 67, 3, 27, 47, 47, 69, 66, 70, 70, 
  23, 7, 50, 50, 70, 67, 72, 72, 40, 43, 11, 30, 70, 67, 72, 72, 40, 
  43, 30, 11}, {7, 27, 43, 43, 60, 57, 61, 61, 23, 3, 40, 40, 57, 56, 
  58, 58, 50, 47, 11, 30, 61, 58, 63, 63, 50, 47, 30, 11, 61, 58, 63, 
  63, 69, 66, 70, 70, 7, 23, 50, 50, 66, 65, 67, 67, 27, 3, 47, 47, 
  70, 67, 72, 72, 43, 40, 11, 30, 70, 67, 72, 72, 43, 40, 30, 11}, {7,
   27, 43, 43, 61, 61, 57, 60, 23, 3, 40, 40, 58, 58, 56, 57, 50, 47, 
  11, 30, 63, 63, 58, 61, 50, 47, 30, 11, 63, 63, 58, 61, 70, 67, 72, 
  72, 11, 30, 40, 43, 70, 67, 72, 72, 30, 11, 40, 43, 66, 65, 67, 67, 
  47, 47, 3, 27, 69, 66, 70, 70, 50, 50, 23, 7}, {7, 27, 43, 43, 61, 
  61, 60, 57, 23, 3, 40, 40, 58, 58, 57, 56, 50, 47, 11, 30, 63, 63, 
  61, 58, 50, 47, 30, 11, 63, 63, 61, 58, 70, 67, 72, 72, 11, 30, 43, 
  40, 70, 67, 72, 72, 30, 11, 43, 40, 69, 66, 70, 70, 50, 50, 7, 23, 
  66, 65, 67, 67, 47, 47, 27, 3}, {2, 26, 45, 46, 47, 47, 47, 47, 22, 
  6, 46, 49, 50, 50, 50, 50, 38, 39, 2, 22, 40, 40, 40, 40, 39, 42, 
  26, 6, 43, 43, 43, 43, 40, 43, 47, 50, 15, 33, 52, 52, 40, 43, 47, 
  50, 33, 15, 52, 52, 40, 43, 47, 50, 52, 52, 15, 33, 40, 43, 47, 50, 
  52, 52, 33, 15}, {6, 22, 46, 49, 50, 50, 50, 50, 26, 2, 45, 46, 47, 
  47, 47, 47, 39, 38, 2, 22, 40, 40, 40, 40, 42, 39, 26, 6, 43, 43, 
  43, 43, 43, 40, 47, 50, 15, 33, 52, 52, 43, 40, 47, 50, 33, 15, 52, 
  52, 43, 40, 47, 50, 52, 52, 15, 33, 43, 40, 47, 50, 52, 52, 33, 
  15}, {9, 10, 2, 6, 11, 11, 11, 11, 10, 9, 2, 6, 11, 11, 11, 11, 2, 
  2, 0, 1, 3, 3, 3, 3, 6, 6, 1, 5, 7, 7, 7, 7, 11, 11, 3, 7, 13, 14, 
  15, 15, 11, 11, 3, 7, 14, 13, 15, 15, 11, 11, 3, 7, 15, 15, 13, 14, 
  11, 11, 3, 7, 15, 15, 14, 13}, {10, 29, 26, 22, 30, 30, 30, 30, 29, 
  10, 26, 22, 30, 30, 30, 30, 22, 22, 1, 21, 23, 23, 23, 23, 26, 26, 
  25, 1, 27, 27, 27, 27, 30, 30, 27, 23, 14, 32, 33, 33, 30, 30, 27, 
  23, 32, 14, 33, 33, 30, 30, 27, 23, 33, 33, 14, 32, 30, 30, 27, 23, 
  33, 33, 32, 14}, {11, 30, 47, 50, 58, 61, 63, 63, 30, 11, 47, 50, 
  58, 61, 63, 63, 40, 40, 3, 23, 56, 57, 58, 58, 43, 43, 27, 7, 57, 
  60, 61, 61, 67, 67, 65, 66, 3, 27, 47, 47, 70, 70, 66, 69, 23, 7, 
  50, 50, 72, 72, 67, 70, 40, 43, 11, 30, 72, 72, 67, 70, 40, 43, 30, 
  11}, {11, 30, 47, 50, 61, 58, 63, 63, 30, 11, 47, 50, 61, 58, 63, 
  63, 40, 40, 3, 23, 57, 56, 58, 58, 43, 43, 27, 7, 60, 57, 61, 61, 
  70, 70, 66, 69, 7, 23, 50, 50, 67, 67, 65, 66, 27, 3, 47, 47, 72, 
  72, 67, 70, 43, 40, 11, 30, 72, 72, 67, 70, 43, 40, 30, 11}, {11, 
  30, 47, 50, 63, 63, 58, 61, 30, 11, 47, 50, 63, 63, 58, 61, 40, 40, 
  3, 23, 58, 58, 56, 57, 43, 43, 27, 7, 61, 61, 57, 60, 72, 72, 67, 
  70, 11, 30, 40, 43, 72, 72, 67, 70, 30, 11, 40, 43, 67, 67, 65, 66, 
  47, 47, 3, 27, 70, 70, 66, 69, 50, 50, 23, 7}, {11, 30, 47, 50, 63, 
  63, 61, 58, 30, 11, 47, 50, 63, 63, 61, 58, 40, 40, 3, 23, 58, 58, 
  57, 56, 43, 43, 27, 7, 61, 61, 60, 57, 72, 72, 67, 70, 11, 30, 43, 
  40, 72, 72, 67, 70, 30, 11, 43, 40, 70, 70, 66, 69, 50, 50, 7, 23, 
  67, 67, 65, 66, 47, 47, 27, 3}, {2, 26, 46, 45, 47, 47, 47, 47, 22, 
  6, 49, 46, 50, 50, 50, 50, 39, 42, 6, 26, 43, 43, 43, 43, 38, 39, 
  22, 2, 40, 40, 40, 40, 40, 43, 50, 47, 15, 33, 52, 52, 40, 43, 50, 
  47, 33, 15, 52, 52, 40, 43, 50, 47, 52, 52, 15, 33, 40, 43, 50, 47, 
  52, 52, 33, 15}, {6, 22, 49, 46, 50, 50, 50, 50, 26, 2, 46, 45, 47, 
  47, 47, 47, 42, 39, 6, 26, 43, 43, 43, 43, 39, 38, 22, 2, 40, 40, 
  40, 40, 43, 40, 50, 47, 15, 33, 52, 52, 43, 40, 50, 47, 33, 15, 52, 
  52, 43, 40, 50, 47, 52, 52, 15, 33, 43, 40, 50, 47, 52, 52, 33, 
  15}, {10, 29, 22, 26, 30, 30, 30, 30, 29, 10, 22, 26, 30, 30, 30, 
  30, 26, 26, 1, 25, 27, 27, 27, 27, 22, 22, 21, 1, 23, 23, 23, 23, 
  30, 30, 23, 27, 14, 32, 33, 33, 30, 30, 23, 27, 32, 14, 33, 33, 30, 
  30, 23, 27, 33, 33, 14, 32, 30, 30, 23, 27, 33, 33, 32, 14}, {9, 10,
   6, 2, 11, 11, 11, 11, 10, 9, 6, 2, 11, 11, 11, 11, 6, 6, 5, 1, 7, 
  7, 7, 7, 2, 2, 1, 0, 3, 3, 3, 3, 11, 11, 7, 3, 13, 14, 15, 15, 11, 
  11, 7, 3, 14, 13, 15, 15, 11, 11, 7, 3, 15, 15, 13, 14, 11, 11, 7, 
  3, 15, 15, 14, 13}, {11, 30, 50, 47, 58, 61, 63, 63, 30, 11, 50, 47,
   58, 61, 63, 63, 43, 43, 7, 27, 57, 60, 61, 61, 40, 40, 23, 3, 56, 
  57, 58, 58, 67, 67, 66, 65, 3, 27, 47, 47, 70, 70, 69, 66, 23, 7, 
  50, 50, 72, 72, 70, 67, 40, 43, 11, 30, 72, 72, 70, 67, 40, 43, 30, 
  11}, {11, 30, 50, 47, 61, 58, 63, 63, 30, 11, 50, 47, 61, 58, 63, 
  63, 43, 43, 7, 27, 60, 57, 61, 61, 40, 40, 23, 3, 57, 56, 58, 58, 
  70, 70, 69, 66, 7, 23, 50, 50, 67, 67, 66, 65, 27, 3, 47, 47, 72, 
  72, 70, 67, 43, 40, 11, 30, 72, 72, 70, 67, 43, 40, 30, 11}, {11, 
  30, 50, 47, 63, 63, 58, 61, 30, 11, 50, 47, 63, 63, 58, 61, 43, 43, 
  7, 27, 61, 61, 57, 60, 40, 40, 23, 3, 58, 58, 56, 57, 72, 72, 70, 
  67, 11, 30, 40, 43, 72, 72, 70, 67, 30, 11, 40, 43, 67, 67, 66, 65, 
  47, 47, 3, 27, 70, 70, 69, 66, 50, 50, 23, 7}, {11, 30, 50, 47, 63, 
  63, 61, 58, 30, 11, 50, 47, 63, 63, 61, 58, 43, 43, 7, 27, 61, 61, 
  60, 57, 40, 40, 23, 3, 58, 58, 57, 56, 72, 72, 70, 67, 11, 30, 43, 
  40, 72, 72, 70, 67, 30, 11, 43, 40, 70, 70, 69, 66, 50, 50, 7, 23, 
  67, 67, 66, 65, 47, 47, 27, 3}, {3, 27, 47, 47, 65, 66, 67, 67, 23, 
  7, 50, 50, 66, 69, 70, 70, 40, 43, 11, 30, 67, 70, 72, 72, 40, 43, 
  30, 11, 67, 70, 72, 72, 56, 57, 58, 58, 3, 23, 40, 40, 57, 60, 61, 
  61, 27, 7, 43, 43, 58, 61, 63, 63, 47, 50, 11, 30, 58, 61, 63, 63, 
  47, 50, 30, 11}, {7, 23, 50, 50, 66, 69, 70, 70, 27, 3, 47, 47, 65, 
  66, 67, 67, 43, 40, 11, 30, 67, 70, 72, 72, 43, 40, 30, 11, 67, 70, 
  72, 72, 57, 56, 58, 58, 3, 23, 40, 40, 60, 57, 61, 61, 27, 7, 43, 
  43, 61, 58, 63, 63, 47, 50, 11, 30, 61, 58, 63, 63, 47, 50, 30, 
  11}, {11, 30, 40, 43, 67, 70, 72, 72, 30, 11, 40, 43, 67, 70, 72, 
  72, 47, 47, 3, 27, 65, 66, 67, 67, 50, 50, 23, 7, 66, 69, 70, 70, 
  58, 58, 56, 57, 3, 23, 40, 40, 61, 61, 57, 60, 27, 7, 43, 43, 63, 
  63, 58, 61, 47, 50, 11, 30, 63, 63, 58, 61, 47, 50, 30, 11}, {11, 
  30, 43, 40, 67, 70, 72, 72, 30, 11, 43, 40, 67, 70, 72, 72, 50, 50, 
  7, 23, 66, 69, 70, 70, 47, 47, 27, 3, 65, 66, 67, 67, 58, 58, 57, 
  56, 3, 23, 40, 40, 61, 61, 60, 57, 27, 7, 43, 43, 63, 63, 61, 58, 
  47, 50, 11, 30, 63, 63, 61, 58, 47, 50, 30, 11}, {13, 14, 15, 15, 3,
   7, 11, 11, 14, 13, 15, 15, 3, 7, 11, 11, 15, 15, 13, 14, 3, 7, 11, 
  11, 15, 15, 14, 13, 3, 7, 11, 11, 3, 3, 3, 3, 0, 1, 2, 2, 7, 7, 7, 
  7, 1, 5, 6, 6, 11, 11, 11, 11, 2, 6, 9, 10, 11, 11, 11, 11, 2, 6, 
  10, 9}, {14, 32, 33, 33, 27, 23, 30, 30, 32, 14, 33, 33, 27, 23, 30,
   30, 33, 33, 14, 32, 27, 23, 30, 30, 33, 33, 32, 14, 27, 23, 30, 30,
   23, 23, 23, 23, 1, 21, 22, 22, 27, 27, 27, 27, 25, 1, 26, 26, 30, 
  30, 30, 30, 26, 22, 10, 29, 30, 30, 30, 30, 26, 22, 29, 10}, {15, 
  33, 52, 52, 47, 50, 40, 43, 33, 15, 52, 52, 47, 50, 40, 43, 52, 52, 
  15, 33, 47, 50, 40, 43, 52, 52, 33, 15, 47, 50, 40, 43, 40, 40, 40, 
  40, 2, 22, 38, 39, 43, 43, 43, 43, 26, 6, 39, 42, 47, 47, 47, 47, 
  45, 46, 2, 26, 50, 50, 50, 50, 46, 49, 22, 6}, {15, 33, 52, 52, 47, 
  50, 43, 40, 33, 15, 52, 52, 47, 50, 43, 40, 52, 52, 15, 33, 47, 50, 
  43, 40, 52, 52, 33, 15, 47, 50, 43, 40, 40, 40, 40, 40, 2, 22, 39, 
  38, 43, 43, 43, 43, 26, 6, 42, 39, 50, 50, 50, 50, 46, 49, 6, 22, 
  47, 47, 47, 47, 45, 46, 26, 2}, {3, 27, 47, 47, 66, 65, 67, 67, 23, 
  7, 50, 50, 69, 66, 70, 70, 40, 43, 11, 30, 70, 67, 72, 72, 40, 43, 
  30, 11, 70, 67, 72, 72, 57, 60, 61, 61, 7, 27, 43, 43, 56, 57, 58, 
  58, 23, 3, 40, 40, 58, 61, 63, 63, 50, 47, 11, 30, 58, 61, 63, 63, 
  50, 47, 30, 11}, {7, 23, 50, 50, 69, 66, 70, 70, 27, 3, 47, 47, 66, 
  65, 67, 67, 43, 40, 11, 30, 70, 67, 72, 72, 43, 40, 30, 11, 70, 67, 
  72, 72, 60, 57, 61, 61, 7, 27, 43, 43, 57, 56, 58, 58, 23, 3, 40, 
  40, 61, 58, 63, 63, 50, 47, 11, 30, 61, 58, 63, 63, 50, 47, 30, 
  11}, {11, 30, 40, 43, 70, 67, 72, 72, 30, 11, 40, 43, 70, 67, 72, 
  72, 47, 47, 3, 27, 66, 65, 67, 67, 50, 50, 23, 7, 69, 66, 70, 70, 
  61, 61, 57, 60, 7, 27, 43, 43, 58, 58, 56, 57, 23, 3, 40, 40, 63, 
  63, 58, 61, 50, 47, 11, 30, 63, 63, 58, 61, 50, 47, 30, 11}, {11, 
  30, 43, 40, 70, 67, 72, 72, 30, 11, 43, 40, 70, 67, 72, 72, 50, 50, 
  7, 23, 69, 66, 70, 70, 47, 47, 27, 3, 66, 65, 67, 67, 61, 61, 60, 
  57, 7, 27, 43, 43, 58, 58, 57, 56, 23, 3, 40, 40, 63, 63, 61, 58, 
  50, 47, 11, 30, 63, 63, 61, 58, 50, 47, 30, 11}, {14, 32, 33, 33, 
  23, 27, 30, 30, 32, 14, 33, 33, 23, 27, 30, 30, 33, 33, 14, 32, 23, 
  27, 30, 30, 33, 33, 32, 14, 23, 27, 30, 30, 27, 27, 27, 27, 1, 25, 
  26, 26, 23, 23, 23, 23, 21, 1, 22, 22, 30, 30, 30, 30, 22, 26, 10, 
  29, 30, 30, 30, 30, 22, 26, 29, 10}, {13, 14, 15, 15, 7, 3, 11, 11, 
  14, 13, 15, 15, 7, 3, 11, 11, 15, 15, 13, 14, 7, 3, 11, 11, 15, 15, 
  14, 13, 7, 3, 11, 11, 7, 7, 7, 7, 5, 1, 6, 6, 3, 3, 3, 3, 1, 0, 2, 
  2, 11, 11, 11, 11, 6, 2, 9, 10, 11, 11, 11, 11, 6, 2, 10, 9}, {15, 
  33, 52, 52, 50, 47, 40, 43, 33, 15, 52, 52, 50, 47, 40, 43, 52, 52, 
  15, 33, 50, 47, 40, 43, 52, 52, 33, 15, 50, 47, 40, 43, 43, 43, 43, 
  43, 6, 26, 39, 42, 40, 40, 40, 40, 22, 2, 38, 39, 47, 47, 47, 47, 
  46, 45, 2, 26, 50, 50, 50, 50, 49, 46, 22, 6}, {15, 33, 52, 52, 50, 
  47, 43, 40, 33, 15, 52, 52, 50, 47, 43, 40, 52, 52, 15, 33, 50, 47, 
  43, 40, 52, 52, 33, 15, 50, 47, 43, 40, 43, 43, 43, 43, 6, 26, 42, 
  39, 40, 40, 40, 40, 22, 2, 39, 38, 50, 50, 50, 50, 49, 46, 6, 22, 
  47, 47, 47, 47, 46, 45, 26, 2}, {3, 27, 47, 47, 67, 67, 65, 66, 23, 
  7, 50, 50, 70, 70, 66, 69, 40, 43, 11, 30, 72, 72, 67, 70, 40, 43, 
  30, 11, 72, 72, 67, 70, 58, 61, 63, 63, 11, 30, 47, 50, 58, 61, 63, 
  63, 30, 11, 47, 50, 56, 57, 58, 58, 40, 40, 3, 23, 57, 60, 61, 61, 
  43, 43, 27, 7}, {7, 23, 50, 50, 70, 70, 66, 69, 27, 3, 47, 47, 67, 
  67, 65, 66, 43, 40, 11, 30, 72, 72, 67, 70, 43, 40, 30, 11, 72, 72, 
  67, 70, 61, 58, 63, 63, 11, 30, 47, 50, 61, 58, 63, 63, 30, 11, 47, 
  50, 57, 56, 58, 58, 40, 40, 3, 23, 60, 57, 61, 61, 43, 43, 27, 
  7}, {11, 30, 40, 43, 72, 72, 67, 70, 30, 11, 40, 43, 72, 72, 67, 70,
   47, 47, 3, 27, 67, 67, 65, 66, 50, 50, 23, 7, 70, 70, 66, 69, 63, 
  63, 58, 61, 11, 30, 47, 50, 63, 63, 58, 61, 30, 11, 47, 50, 58, 58, 
  56, 57, 40, 40, 3, 23, 61, 61, 57, 60, 43, 43, 27, 7}, {11, 30, 43, 
  40, 72, 72, 67, 70, 30, 11, 43, 40, 72, 72, 67, 70, 50, 50, 7, 23, 
  70, 70, 66, 69, 47, 47, 27, 3, 67, 67, 65, 66, 63, 63, 61, 58, 11, 
  30, 47, 50, 63, 63, 61, 58, 30, 11, 47, 50, 58, 58, 57, 56, 40, 40, 
  3, 23, 61, 61, 60, 57, 43, 43, 27, 7}, {15, 33, 52, 52, 40, 43, 47, 
  50, 33, 15, 52, 52, 40, 43, 47, 50, 52, 52, 15, 33, 40, 43, 47, 50, 
  52, 52, 33, 15, 40, 43, 47, 50, 47, 47, 47, 47, 2, 26, 45, 46, 50, 
  50, 50, 50, 22, 6, 46, 49, 40, 40, 40, 40, 38, 39, 2, 22, 43, 43, 
  43, 43, 39, 42, 26, 6}, {15, 33, 52, 52, 43, 40, 47, 50, 33, 15, 52,
   52, 43, 40, 47, 50, 52, 52, 15, 33, 43, 40, 47, 50, 52, 52, 33, 15,
   43, 40, 47, 50, 50, 50, 50, 50, 6, 22, 46, 49, 47, 47, 47, 47, 26, 
  2, 45, 46, 40, 40, 40, 40, 39, 38, 2, 22, 43, 43, 43, 43, 42, 39, 
  26, 6}, {13, 14, 15, 15, 11, 11, 3, 7, 14, 13, 15, 15, 11, 11, 3, 7,
   15, 15, 13, 14, 11, 11, 3, 7, 15, 15, 14, 13, 11, 11, 3, 7, 11, 11,
   11, 11, 9, 10, 2, 6, 11, 11, 11, 11, 10, 9, 2, 6, 3, 3, 3, 3, 2, 2,
   0, 1, 7, 7, 7, 7, 6, 6, 1, 5}, {14, 32, 33, 33, 30, 30, 27, 23, 32,
   14, 33, 33, 30, 30, 27, 23, 33, 33, 14, 32, 30, 30, 27, 23, 33, 33,
   32, 14, 30, 30, 27, 23, 30, 30, 30, 30, 10, 29, 26, 22, 30, 30, 30,
   30, 29, 10, 26, 22, 23, 23, 23, 23, 22, 22, 1, 21, 27, 27, 27, 27, 
  26, 26, 25, 1}, {3, 27, 47, 47, 67, 67, 66, 65, 23, 7, 50, 50, 70, 
  70, 69, 66, 40, 43, 11, 30, 72, 72, 70, 67, 40, 43, 30, 11, 72, 72, 
  70, 67, 58, 61, 63, 63, 11, 30, 50, 47, 58, 61, 63, 63, 30, 11, 50, 
  47, 57, 60, 61, 61, 43, 43, 7, 27, 56, 57, 58, 58, 40, 40, 23, 
  3}, {7, 23, 50, 50, 70, 70, 69, 66, 27, 3, 47, 47, 67, 67, 66, 65, 
  43, 40, 11, 30, 72, 72, 70, 67, 43, 40, 30, 11, 72, 72, 70, 67, 61, 
  58, 63, 63, 11, 30, 50, 47, 61, 58, 63, 63, 30, 11, 50, 47, 60, 57, 
  61, 61, 43, 43, 7, 27, 57, 56, 58, 58, 40, 40, 23, 3}, {11, 30, 40, 
  43, 72, 72, 70, 67, 30, 11, 40, 43, 72, 72, 70, 67, 47, 47, 3, 27, 
  67, 67, 66, 65, 50, 50, 23, 7, 70, 70, 69, 66, 63, 63, 58, 61, 11, 
  30, 50, 47, 63, 63, 58, 61, 30, 11, 50, 47, 61, 61, 57, 60, 43, 43, 
  7, 27, 58, 58, 56, 57, 40, 40, 23, 3}, {11, 30, 43, 40, 72, 72, 70, 
  67, 30, 11, 43, 40, 72, 72, 70, 67, 50, 50, 7, 23, 70, 70, 69, 66, 
  47, 47, 27, 3, 67, 67, 66, 65, 63, 63, 61, 58, 11, 30, 50, 47, 63, 
  63, 61, 58, 30, 11, 50, 47, 61, 61, 60, 57, 43, 43, 7, 27, 58, 58, 
  57, 56, 40, 40, 23, 3}, {15, 33, 52, 52, 40, 43, 50, 47, 33, 15, 52,
   52, 40, 43, 50, 47, 52, 52, 15, 33, 40, 43, 50, 47, 52, 52, 33, 15,
   40, 43, 50, 47, 47, 47, 47, 47, 2, 26, 46, 45, 50, 50, 50, 50, 22, 
  6, 49, 46, 43, 43, 43, 43, 39, 42, 6, 26, 40, 40, 40, 40, 38, 39, 
  22, 2}, {15, 33, 52, 52, 43, 40, 50, 47, 33, 15, 52, 52, 43, 40, 50,
   47, 52, 52, 15, 33, 43, 40, 50, 47, 52, 52, 33, 15, 43, 40, 50, 47,
   50, 50, 50, 50, 6, 22, 49, 46, 47, 47, 47, 47, 26, 2, 46, 45, 43, 
  43, 43, 43, 42, 39, 6, 26, 40, 40, 40, 40, 39, 38, 22, 2}, {14, 32, 
  33, 33, 30, 30, 23, 27, 32, 14, 33, 33, 30, 30, 23, 27, 33, 33, 14, 
  32, 30, 30, 23, 27, 33, 33, 32, 14, 30, 30, 23, 27, 30, 30, 30, 30, 
  10, 29, 22, 26, 30, 30, 30, 30, 29, 10, 22, 26, 27, 27, 27, 27, 26, 
  26, 1, 25, 23, 23, 23, 23, 22, 22, 21, 1}, {13, 14, 15, 15, 11, 11, 
  7, 3, 14, 13, 15, 15, 11, 11, 7, 3, 15, 15, 13, 14, 11, 11, 7, 3, 
  15, 15, 14, 13, 11, 11, 7, 3, 11, 11, 11, 11, 9, 10, 6, 2, 11, 11, 
  11, 11, 10, 9, 6, 2, 7, 7, 7, 7, 6, 6, 5, 1, 3, 3, 3, 3, 2, 2, 1, 
  0}};

const int probabilityData<16>::infiniteMask[][16] = 
		{
			{0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
			{1, 0, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
			{2, 2, 0, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
			{2, 2, 1, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
			{3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3},
			{3, 3, 3, 3, 1, 0, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3},
			{3, 3, 3, 3, 2, 2, 0, 1, 3, 3, 3, 3, 3, 3, 3, 3},
			{3, 3, 3, 3, 2, 2, 1, 0, 3, 3, 3, 3, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 1, 0, 2, 2, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 1, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 0, 3, 3, 3, 3},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 1, 2, 2},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 0, 2, 2},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 1},
      {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1, 0}
		};
template<> void genotypeProbabilitiesNoIntercross<16, true>(double(&prob)[nDifferentProbs], double r, int, std::size_t)
{
	prob[0] = (1-r)*(1-r)*(1-r)/(16*(1 + 2*r));
	prob[1] = (1-r)*(1-r)*r/(16+32*r);
	prob[2] = (1 - r)*r/(32+64*r);
  prob[3] = r/(64*(1 + 2 * r));
}
template<> void genotypeProbabilitiesNoIntercross<16, false>(double(&prob)[nDifferentProbs], double r, int selfingGenerations, std::size_t nFunnels)
{
	double onePlus2R = 1 + 2 *r;
	double rSquared = r*r;
	double oneMinusRSquared = (1-r)*(1-r);
	double powOneMinus2R = std::pow(1 - 2 * r, selfingGenerations);
	double powD1 = std::pow(1 + 2*(-1 + r)*r, selfingGenerations);
	double oneMinusRCubed = (1-r)*oneMinusRSquared;
  double oneMinusRPow4 = oneMinusRSquared*oneMinusRSquared;
  double oneMinusRPow5 = oneMinusRCubed*oneMinusRSquared;
  double oneMinusRPow6 = oneMinusRCubed*oneMinusRCubed;
	double oneMinusR = 1 - r;
	double pow2 = std::pow(2, selfingGenerations);
    
  prob[0] = -(oneMinusRCubed*(2 - 2*pow2 - onePlus2R*powD1 + powOneMinus2R + 4*r - 2*powOneMinus2R*r))/(32*onePlus2R*pow2);
  prob[1] = 0;
  prob[2] = 0;
  prob[3] = 0;
  prob[4] = -(oneMinusRCubed*(-1 + powD1))/(256*pow2);
  prob[5] = (oneMinusRSquared*r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(32*onePlus2R*pow2);
  prob[6] = 0;
  prob[7] = 0;
  prob[8] = -(oneMinusRSquared*(-1 + powD1)*r)/(256*pow2);
  prob[9] = (oneMinusR*r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(64*onePlus2R*pow2);
  prob[10] = 0;
  prob[11] = 0;
  prob[12] = -(oneMinusR*(-1 + powD1)*r)/(512*pow2);
  prob[13] = (r*(-2 + 2*pow2 + onePlus2R*powD1 - powOneMinus2R - 4*r + 2*powOneMinus2R*r))/(128*onePlus2R*pow2);
  prob[14] = 0;
  prob[15] = 0;
  prob[16] = -((-1 + powD1)*r)/(1024*pow2);
  prob[17] = (-2 + onePlus2R*powD1 + powOneMinus2R - 4*r + 4*pow2*r - 2*powOneMinus2R*r)/(256*onePlus2R*pow2);
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
  prob[34] = 0;
  prob[35] = 0;
  prob[36] = 0;
  prob[37] = 0;
  prob[38] = 0;
  prob[39] = 0;
  prob[40] = 0;
  prob[41] = 0;
  prob[42] = 0;
  prob[43] = 0;
  prob[44] = 0;
  prob[45] = 0;
  prob[46] = 0;
  prob[47] = 0;
  prob[48] = 0;
  prob[49] = 0;
  prob[50] = 0;
  prob[51] = 0;
  prob[52] = 0;
  prob[53] = 0;
  prob[54] = 0;
  prob[55] = 0;
  prob[56] = 0;
  prob[57] = 0;
  prob[58] = 0;
  prob[59] = 0;
  prob[60] = 0;
  prob[61] = 0;
  prob[62] = 0;
  prob[63] = 0;
  prob[64] = 0;
  prob[65] = 0;
  prob[66] = 0;
  prob[67] = 0;
  prob[68] = 0;
  prob[69] = 0;
  prob[70] = 0;
  prob[71] = 0;
  prob[72] = 0;
  prob[73] = 0;
  prob[74] = 0;
  prob[75] = (oneMinusRPow6*(powD1 + powOneMinus2R))/(256*pow2);
  prob[76] = (oneMinusRPow5*(powD1 + powOneMinus2R)*r)/(256*pow2);
  prob[77] = (oneMinusRPow4*(powD1 + powOneMinus2R)*r)/(512*pow2);
  prob[78] = (oneMinusRCubed*(powD1 + powOneMinus2R)*r)/(1024*pow2);
  prob[79] = (oneMinusRPow4*(powD1 + powOneMinus2R)*rSquared)/(256*pow2);
  prob[80] = (oneMinusRCubed*(powD1 + powOneMinus2R)*rSquared)/(512*pow2);
  prob[81] = (oneMinusRSquared*(powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
  prob[82] = (oneMinusRSquared*(powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
  prob[83] = (oneMinusR*(powD1 + powOneMinus2R)*rSquared)/(2048*pow2);
  prob[84] = ((powD1 + powOneMinus2R)*rSquared)/(4096*pow2);
  prob[85] = -(oneMinusRPow6*(-powD1 + powOneMinus2R))/(256*pow2);
  prob[86] = -(oneMinusRPow5*(-powD1 + powOneMinus2R)*r)/(256*pow2);
  prob[87] = -(oneMinusRPow4*(-powD1 + powOneMinus2R)*r)/(512*pow2);
  prob[88] = -(oneMinusRCubed*(-powD1 + powOneMinus2R)*r)/(1024*pow2);
  prob[89] = -(oneMinusRPow4*(-powD1 + powOneMinus2R)*rSquared)/(256*pow2);
  prob[90] = -(oneMinusRCubed*(-powD1 + powOneMinus2R)*rSquared)/(512*pow2);
  prob[91] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
  prob[92] = -(oneMinusRSquared*(-powD1 + powOneMinus2R)*rSquared)/(1024*pow2);
  prob[93] = -(oneMinusR*(-powD1 + powOneMinus2R)*rSquared)/(2048*pow2);
  prob[94] = ((powD1 - powOneMinus2R)*rSquared)/(4096*pow2);
#ifndef NDEBUG
  double sum = 0;
  for(int i = 0; i < 95; i++) sum += prob[i];
#endif
}
template<> void genotypeProbabilitiesWithIntercross<16, true>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int, std::size_t nFunnels)
{
	double powOneMinusR = std::pow(1 - r, nAIGenerations-1);
	prob[0] = (powOneMinusR*(1-r)*(1-r)*(1-r)*(1-r)/16 + (2 * r + 1 - powOneMinusR)/(16*16))/(1 + 2 * r);
	prob[1] = prob[2] = prob[3] = (1 - 16*prob[0])/(16*15);
}
template<> void genotypeProbabilitiesWithIntercross<16, false>(double(&prob)[nDifferentProbs], int nAIGenerations, double r, int selfingGenerations, std::size_t nFunnels)
{
	if (nFunnels == 1)
	{
		throw std::runtime_error("The case of a single funnel with intercrossing is not yet handled");
	}
  //See mathematica code for rewritting rules
	double oneMinusR = (1 - r);
  double onePlus2R = 1 + 2*r;
  double oneMinusRSquared = oneMinusR*oneMinusR;
  double twoRMinus3 = 2 * r - 3;
  double twoRMinus3Squared = twoRMinus3*twoRMinus3;
  double powOneMinusR1 = std::pow(oneMinusR, 1+nAIGenerations);
  double powOneMinusR2 = std::pow(oneMinusR, nAIGenerations - 1);
  double powOneMinusR3 = std::pow(oneMinusR, nAIGenerations - 2);
  double oneMinusRPow1 = std::pow(oneMinusR, 2*nAIGenerations);
  double pow2 = std::pow(2, selfingGenerations);
  double twoRMinus1 = 2*r - 1;
  double complexPower1 = std::pow(1 + 2 *(-1 + r)* r, selfingGenerations);
  double quadratic1 = 5 + 4 * (-2 + r) * r;
  double quadratic1Squared = quadratic1*quadratic1;
  double twoRMinus1Squared = twoRMinus1*twoRMinus1;
  double twoRMinus1PowD = std::pow(twoRMinus1, selfingGenerations);
  double rMinus1Pow4 = std::pow(r - 1, 4);
  double complexPart1 = 1.0/256.0 + powOneMinusR2 * (-1.0/256.0 + (1.0/16.0) *rMinus1Pow4);
  double complexPart2 = -7 + std::pow(2.0, 3.0 + selfingGenerations) - 4*twoRMinus1*twoRMinus1PowD - 14*r;
  double complexPart3 = 1.0 - 16.0*complexPart1;
  double complexPart3Squared = complexPart3*complexPart3;
  double complexPart1Squared = complexPart1*complexPart1;
  double complexPart4 = complexPower1 + 8*powOneMinusR2*twoRMinus1PowD;
  double complexPart5 = 105*complexPower1 - 16*(-7*complexPower1 + 8*twoRMinus1PowD)*(-2 + r)*r*(2 + (-2 + r)*r);
  double complexPart6 = complexPart3*(15*complexPart4 + powOneMinusR2*(105*complexPower1 + 16*(-2 + r)*r*(2 + (-2 + r)*r)*(7*complexPower1 + 8*twoRMinus1PowD)));
  double complexPart7 = complexPart3*(complexPart5*powOneMinusR2 + 15*(complexPower1 - 8*powOneMinusR2*twoRMinus1PowD));
  double complexPart8 = 15*(-1 + r) + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3;
  double complexPart8Squared = complexPart8;
  double complexPart9 = complexPart8*(-240*oneMinusR + complexPower1*(225*oneMinusR + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3));
  double complexPart10 = 240*oneMinusRSquared + 112*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 - complexPower1*(225*oneMinusRSquared + 98*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared / twoRMinus3Squared);
  double complexPart11 = 57600*complexPart1Squared*(complexPower1 - twoRMinus1PowD) + complexPart3Squared*(complexPower1 + twoRMinus1PowD);
  double complexPart12 = complexPart3Squared*(complexPower1 - twoRMinus1PowD) + 57600*complexPart1Squared*(complexPower1 + twoRMinus1PowD);
  double complexPart13 = std::pow(0.5 - oneMinusR*r, selfingGenerations);
  double complexPart14 = -15*(-15 + std::pow(2.0, 3 + selfingGenerations))*oneMinusR *onePlus2R + complexPart2*oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3;
  double complexPart15 = 32*((-15 + std::pow(2.0, 3 + selfingGenerations))*oneMinusRSquared*onePlus2R + complexPart2*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3) + complexPower1*onePlus2R*(225*oneMinusRSquared + 98*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared / twoRMinus3Squared);
  double complexPart16 = -480*complexPart14*oneMinusR + complexPower1*onePlus2R*(50625*oneMinusRSquared - 1470*powOneMinusR1*quadratic1*twoRMinus1*twoRMinus3 + oneMinusRPow1*quadratic1Squared*twoRMinus1Squared / twoRMinus3Squared);
  double complexPart17 = 240 - 225*complexPower1 + 1680*powOneMinusR2 + powOneMinusR3*(1792*oneMinusR*(-2 + r)*r*(2 + (-2 + r)*r) - complexPower1*quadratic1*twoRMinus1*twoRMinus3*(98*oneMinusR + oneMinusR*powOneMinusR2*quadratic1*twoRMinus1*twoRMinus3));
  double complexPart18 = std::pow(4 - 8*oneMinusR*r, selfingGenerations);
  double complexPart19 = std::pow(2.0, -3*(1 + selfingGenerations));

	prob[0] = complexPart15/(4096*oneMinusRSquared*onePlus2R*pow2);
  prob[1] = complexPart17/(15360*pow2);
  prob[2] = complexPart10/(7680*oneMinusRSquared*pow2);
  prob[3] = complexPart10/(3840*oneMinusRSquared*pow2);
  prob[4] = complexPart10/(1920*oneMinusRSquared*pow2);
  prob[5] = complexPart16/(921600*oneMinusRSquared*onePlus2R*pow2);
  prob[6] = complexPart9/(115200*oneMinusRSquared*pow2);
  prob[7] = complexPart9/(57600*oneMinusRSquared*pow2);
  prob[8] = complexPart9/(28800*oneMinusRSquared*pow2);
  prob[9] = complexPart16/(460800*oneMinusRSquared*onePlus2R*pow2);
  prob[10] = complexPart9/(230400*oneMinusRSquared*pow2);
  prob[11] = complexPart9/(28800*oneMinusRSquared*pow2);
  prob[12] = complexPart9/(14400*oneMinusRSquared*pow2);
  prob[13] = complexPart16/(230400*oneMinusRSquared*onePlus2R*pow2);
  prob[14] = complexPart9/(115200*oneMinusRSquared*pow2);
  prob[15] = complexPart9/(57600*oneMinusRSquared*pow2);
  prob[16] = complexPart9/(7200*oneMinusRSquared*pow2);
  prob[17] = complexPart16/(115200*oneMinusRSquared*onePlus2R*pow2);
  prob[18] = complexPart9/(57600*oneMinusRSquared*pow2);
  prob[19] = complexPart9/(28800*oneMinusRSquared*pow2);
  prob[20] = complexPart9/(14400*oneMinusRSquared*pow2);
  prob[21] = complexPart12/(7200*pow2);
  prob[22] = complexPart6/(7200*pow2);
  prob[23] = complexPart6/(3600*pow2);
  prob[24] = complexPart6/(1800*pow2);
  prob[25] = complexPart11/(7200*pow2);
  prob[26] = complexPart7/(7200*pow2);
  prob[27] = complexPart7/(3600*pow2);
  prob[28] = complexPart7/(1800*pow2);
  prob[29] = (complexPart18*complexPart19*complexPart3Squared)/225;
  prob[30] = (2*complexPart3Squared*complexPower1)/(225*pow2);
  prob[31] = (4*complexPart3Squared*complexPower1)/(225*pow2);
  prob[32] = (complexPart3Squared*complexPower1)/(900*pow2);
  prob[33] = (complexPart13*complexPart3Squared)/225;
  prob[34] = (8*complexPart3Squared*complexPower1)/(225*pow2);
  prob[35] = (complexPart3Squared*complexPower1)/(450*pow2);
  prob[36] = (2*complexPart3Squared*complexPower1)/(225*pow2);
  prob[37] = (4*complexPart3Squared*complexPower1)/(225*pow2);
  prob[38] = complexPart12/(3600*pow2);
  prob[39] = complexPart6/(14400*pow2);
  prob[40] = complexPart6/(1800*pow2);
  prob[41] = complexPart6/(900*pow2);
  prob[42] = (complexPart8Squared*complexPower1)/(460800*oneMinusRSquared*pow2);
  prob[43] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[44] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
  prob[45] = complexPart11/(3600*pow2);
  prob[46] = complexPart7/(14400*pow2);
  prob[47] = complexPart7/(1800*pow2);
  prob[48] = complexPart7/(900*pow2);
  prob[49] = (complexPart8Squared*complexPower1)/(460800*oneMinusRSquared*pow2);
  prob[50] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[51] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
  prob[52] = (complexPart13*complexPart3Squared)/225;
  prob[53] = (16*complexPart3Squared*complexPower1)/(225*pow2);
  prob[54] = (2*complexPart3Squared*complexPower1)/(225*pow2);
  prob[55] = (8*complexPart3Squared*complexPower1)/(225*pow2);
  prob[56] = complexPart12/(1800*pow2);
  prob[57] = complexPart6/(7200*pow2);
  prob[58] = complexPart6/(3600*pow2);
  prob[59] = complexPart6/(450*pow2);
  prob[60] = (complexPart8Squared*complexPower1)/(230400*oneMinusRSquared*pow2);
  prob[61] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
  prob[62] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
  prob[63] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
  prob[64] = (complexPart8Squared*complexPower1)/(3600*oneMinusRSquared*pow2);
  prob[65] = complexPart11/(1800*pow2);
  prob[66] = complexPart7/(7200*pow2);
  prob[67] = complexPart7/(3600*pow2);
  prob[68] = complexPart7/(450*pow2);
  prob[69] = (complexPart8Squared*complexPower1)/(230400*oneMinusRSquared*pow2);
  prob[70] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
  prob[71] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
  prob[72] = (complexPart8Squared*complexPower1)/(57600*oneMinusRSquared*pow2);
  prob[73] = (complexPart8Squared*complexPower1)/(3600*oneMinusRSquared*pow2);
  prob[74] = (8*complexPart3Squared*complexPower1)/(225*pow2);
  prob[75] = complexPart12/(900*pow2);
  prob[76] = complexPart6/(3600*pow2);
  prob[77] = complexPart6/(1800*pow2);
  prob[78] = complexPart6/(900*pow2);
  prob[79] = (complexPart8Squared*complexPower1)/(115200*oneMinusRSquared*pow2);
  prob[80] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[81] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
  prob[82] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[83] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
  prob[84] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
  prob[85] = complexPart11/(900*pow2);
  prob[86] = complexPart7/(3600*pow2);
  prob[87] = complexPart7/(1800*pow2);
  prob[88] = complexPart7/(900*pow2);
  prob[89] = (complexPart8Squared*complexPower1)/(115200*oneMinusRSquared*pow2);
  prob[90] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[91] = (complexPart8Squared*complexPower1)/(14400*oneMinusRSquared*pow2);
  prob[92] = (complexPart8Squared*complexPower1)/(28800*oneMinusRSquared*pow2);
  prob[93] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
  prob[94] = (complexPart8Squared*complexPower1)/(7200*oneMinusRSquared*pow2);
#ifndef NDEBUG
	double sum = 0;
	for(int i = 0; i < 95; i++) sum += prob[i];
#endif
	//This is because we combined some states (see mathematica code)
  prob[0] /= 16;
  prob[1] /= 64;
  prob[2] /= 128;
  prob[3] /= 256;
  prob[4] /= 512;
  prob[5] /= 16;
  prob[6] /= 128;
  prob[7] /= 256;
  prob[8] /= 512;
  prob[9] /= 32;
  prob[10] /= 64;
  prob[11] /= 512;
  prob[12] /= 1024;
  prob[13] /= 64;
  prob[14] /= 128;
  prob[15] /= 256;
  prob[16] /= 2048;
  prob[17] /= 128;
  prob[18] /= 256;
  prob[19] /= 512;
  prob[20] /= 1024;
  prob[21] /= 16;
  prob[22] /= 128;
  prob[23] /= 256;
  prob[24] /= 512;
  prob[25] /= 16;
  prob[26] /= 128;
  prob[27] /= 256;
  prob[28] /= 512;
  prob[29] /= 32;
  prob[30] /= 512;
  prob[31] /= 1024;
  prob[32] /= 64;
  prob[33] /= 256;
  prob[34] /= 2048;
  prob[35] /= 128;
  prob[36] /= 512;
  prob[37] /= 1024;
  prob[38] /= 32;
  prob[39] /= 64;
  prob[40] /= 512;
  prob[41] /= 1024;
  prob[42] /= 32;
  prob[43] /= 512;
  prob[44] /= 1024;
  prob[45] /= 32;
  prob[46] /= 64;
  prob[47] /= 512;
  prob[48] /= 1024;
  prob[49] /= 32;
  prob[50] /= 512;
  prob[51] /= 1024;
  prob[52] /= 256;
  prob[53] /= 4096;
  prob[54] /= 512;
  prob[55] /= 2048;
  prob[56] /= 64;
  prob[57] /= 128;
  prob[58] /= 256;
  prob[59] /= 2048;
  prob[60] /= 64;
  prob[61] /= 256;
  prob[62] /= 2048;
  prob[63] /= 256;
  prob[64] /= 4096;
  prob[65] /= 64;
  prob[66] /= 128;
  prob[67] /= 256;
  prob[68] /= 2048;
  prob[69] /= 64;
  prob[70] /= 256;
  prob[71] /= 2048;
  prob[72] /= 256;
  prob[73] /= 4096;
  prob[74] /= 2048;
  prob[75] /= 128;
  prob[76] /= 256;
  prob[77] /= 512;
  prob[78] /= 1024;
  prob[79] /= 128;
  prob[80] /= 512;
  prob[81] /= 1024;
  prob[82] /= 512;
  prob[83] /= 2048;
  prob[84] /= 2048;
  prob[85] /= 128;
  prob[86] /= 256;
  prob[87] /= 512;
  prob[88] /= 1024;
  prob[89] /= 128;
  prob[90] /= 512;
  prob[91] /= 1024;
  prob[92] /= 512;
  prob[93] /= 2048;
  prob[94] /= 2048;
}
#endif
