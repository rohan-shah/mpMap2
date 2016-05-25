#include "funnelsToUniqueValues.h"
void funnelsToUniqueValues(std::map<funnelEncoding, funnelID>& funnelTranslation, std::vector<funnelID>& lineFunnelIDs, std::vector<funnelEncoding>& lineFunnelEncodings, std::vector<funnelEncoding>& allFunnelEncodings, std::vector<funnelType>& lineFunnels, std::vector<funnelType>& allFunnels, int nFounders)
{
	//First the line Funnels
	funnelTranslation.clear();
	for(std::vector<funnelType>::iterator i = lineFunnels.begin(); i != lineFunnels.end(); i++)
	{
		funnelType funnel = *i;
		std::size_t encoded = 0;
		bool isZero = true;
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encoded += (((std::size_t)funnel.val[founderCounter]-(std::size_t)1) << (4 * founderCounter));
			isZero &= (funnel.val[founderCounter] == 0);
		}
		//If it's a dummy value, ignore it
		if(isZero)
		{
			lineFunnelIDs.push_back(-1);
		}
		else
		{
			if(funnelTranslation.find(encoded) == funnelTranslation.end())
			{
				lineFunnelIDs.push_back((int)funnelTranslation.size());
				funnelTranslation.insert(std::make_pair(encoded, (int)funnelTranslation.size()));
				lineFunnelEncodings.push_back(encoded);
			}
			else
			{
				lineFunnelIDs.push_back(funnelTranslation.find(encoded)->second);
			}
		}
	}
	//Then the set of all funnels involved
	funnelTranslation.clear();
	for(std::vector<funnelType>::iterator i = allFunnels.begin(); i != allFunnels.end(); i++)
	{
		funnelType funnel = *i;
		std::size_t encoded = 0;
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encoded += (((std::size_t)funnel.val[founderCounter]-(std::size_t)1) << (4 * founderCounter));
		}
		if(funnelTranslation.find(encoded) == funnelTranslation.end())
		{
			funnelTranslation.insert(std::make_pair(encoded, (int)funnelTranslation.size()));
			allFunnelEncodings.push_back(encoded);
		}
	}

}
