#ifndef FUNNELS_TO_UNIQUE_VALUES_HEADER_GUARD
#define FUNNELS_TO_UNIQUE_VALUES_HEADER_GUARD
#include <map>
#include <vector>
#include "estimateRF.h"
#include "unitTypes.h"
void funnelsToUniqueValues(std::map<funnelEncoding, funnelID>& funnelTranslation, std::vector<funnelID>& lineFunnelIDs, std::vector<funnelEncoding>& lineFunnelEncodings, std::vector<funnelEncoding>& allFunnelEncodings, std::vector<funnelType>& lineFunnels, std::vector<funnelType>& allFunnels, int nFounders);
#endif
