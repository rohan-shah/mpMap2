#ifndef VITERBI_HEADER_GUARD
#define VITERBI_HEADER_GUARD
template<int nFounders, bool infiniteSelfing> struct viterbiAlgorithm;
#include <stdexcept>
class impossibleDataException : public std::runtime_error
{
public:
	impossibleDataException(int marker, int line)
		: std::runtime_error("Data is impossible for the given map"), marker(marker), line(line)
	{}
	int marker, line;
};
#include "viterbiInfiniteSelfing.hpp"
#include "viterbiFiniteSelfing.hpp"
#endif
