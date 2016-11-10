#include "haldaneToRf.h"
#include <cmath>
void haldaneToRf(const std::vector<double>& inputs, std::vector<double>& outputs)
{
	outputs.clear();
	for(std::vector<double>::const_iterator i = inputs.begin(); i != inputs.end(); i++)
	{
		outputs.push_back(0.5*(1.0 - std::exp(-2.0 * *i / 100.0)));
	}
}
