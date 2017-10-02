#include "getMinAIGenerations.h"
#include <algorithm>
int getMinAIGenerations(const std::vector<int>* intercrossingGenerations)
{
	int minAIGenerations = 0;
	for(std::vector<int>::const_iterator i = intercrossingGenerations->begin(); i != intercrossingGenerations->end(); i++)
	{
		if(*i != 0)
	        {
			if(minAIGenerations == 0) minAIGenerations = *i;
			else minAIGenerations = std::min(minAIGenerations, *i);
		}
	}
	minAIGenerations = std::max(minAIGenerations, 1);
	return minAIGenerations;
}
