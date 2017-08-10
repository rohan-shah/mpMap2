#include "orderFunnel.h"
#include <algorithm>
#include <string>
#include <cstring>
#include <stdexcept>
//The functions in this file put the funnels into cannonical order
void orderFour(int* start)
{
	int min = std::min(start[0], start[1]), max = std::max(start[0], start[1]);
	start[0] = min; start[1] = max;
	min = std::min(start[2], start[3]), max = std::max(start[2], start[3]);
	start[2] = min; start[3] = max;
}
template<int size> void orderByHalves(int* funnel)
{
	//Make a copy of the original.
	int copied[size];
	memcpy(copied, funnel, sizeof(int) * size);
	while(true)
	{
		int* firstHalfMin = std::min_element(copied, copied+size/2), *secondHalfMin = std::min_element(copied+size/2, copied+size);
		//If we've exhausted all the elements looking for a minimum value which is contained in one set of eight but not the other, then both halves are the same. So just return the original again.
		if(*firstHalfMin == std::numeric_limits<int>::max() && *secondHalfMin == std::numeric_limits<int>::max())
		{
			return;
		}
		//If the minimum value is present in both sides, set it to max() and repeat. 
		else if(*firstHalfMin == *secondHalfMin)
		{
			*secondHalfMin = *firstHalfMin = std::numeric_limits<int>::max();
			continue;
		}
		//Otherwise we have a unique minimum value in one of the sides. 
		int* minVal = std::min_element(copied, copied + size);
		//if the smallest value is in the first eight, do nothing
		if ((minVal - copied) / (size / 2) == 0)
		{
			return;
		}
		else
		{
			memcpy(copied, funnel, sizeof(int) * size);
			//if the smallest value is in the second eight, switch the two sets of eight
			memcpy(funnel, copied + (size/2), sizeof(int) * (size/2));
			memcpy(funnel + (size/2), copied, sizeof(int) * (size/2));
			return;
		}
	}
}
//order a set of four (considered as first pair / second pair) so that the pair containing the smallest value comes first, followed by the other pair. Each pair also sorted as (smaller, larger)
void orderFunnel4(int* funnel)
{
	orderByHalves<4>(funnel);
	orderFour(funnel);
}
void orderFunnel8(int* funnel)
{
	orderByHalves<8>(funnel);
	orderFunnel4(funnel);
	orderFunnel4(funnel+4);
}
void orderFunnel16(int* funnel)
{
	orderByHalves<16>(funnel);
	orderFunnel8(funnel);
	orderFunnel8(funnel + 8);
}
void orderFunnel(int* funnel, int nFounders)
{
	if(nFounders == 2)
	{
		int min = std::min(funnel[0], funnel[1]);
		int max = std::max(funnel[0], funnel[1]);
		funnel[0] = min;
		funnel[1] = max;
	}
	else if(nFounders == 4)
	{
		orderFunnel4(funnel);
	}
	else if (nFounders == 8)
	{
		orderFunnel8(funnel);
	}
	else if (nFounders == 16)
	{
		orderFunnel16(funnel);
	}
	else throw std::runtime_error("Input nFounders must be 2, 4, 8 or 16");
}
