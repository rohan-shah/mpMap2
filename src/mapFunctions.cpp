#include "mapFunctions.h"
#include <cmath>
double haldaneToRf(double x)
{
	return 0.5 * (1 - std::exp(-2 * x/100));
}

