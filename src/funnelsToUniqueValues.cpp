#include "funnelsToUniqueValues.h"
void funnelsToUniqueValues(std::map<funnelEncoding, funnelID>& funnelTranslation, std::vector<funnelID>& funnelIDs, std::vector<funnelEncoding>& funnelEncodings, std::vector<funnelType>& allFunnels, int nFounders)
{
	for(std::vector<funnelType>::iterator i = allFunnels.begin(); i != allFunnels.end(); i++)
	{
		funnelType funnel = *i;
		int encoded = 0;
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encoded += ((funnel.val[founderCounter]-1) << (3 * founderCounter));
		}
		if(funnelTranslation.find(encoded) == funnelTranslation.end())
		{
			funnelIDs.push_back((int)funnelTranslation.size());
			funnelTranslation.insert(std::make_pair(encoded, (int)funnelTranslation.size()));
			funnelEncodings.push_back(encoded);
		}
		else
		{
			funnelIDs.push_back(funnelTranslation.find(encoded)->second);
		}
	}
}