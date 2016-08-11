#ifndef IMPOSSIBLE_DATA_EXCEPTION_HEADER_GUARD
#define IMPOSSIBLE_DATA_EXCEPTION_HEADER_GUARD
#include <stdexcept>
class impossibleDataException : public std::runtime_error
{
public:
	impossibleDataException(int marker, int line)
		: std::runtime_error("Data is impossible for the given map"), marker(marker), line(line)
	{}
	int marker, line;
};
#endif
