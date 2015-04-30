#ifndef FUNNELS_TO_UNIQUE_VALUES_HEADER_GUARD
#define FUNNELS_TO_UNIQUE_VALUES_HEADER_GUARD
#include <map>
#include <vector>
#include "estimateRf.h"
#include "unitTypes.hpp"
void funnelsToUniqueValues(std::map<funnelEncoding, funnelID>& funnelTranslation, std::vector<funnelID>& funnelIDs, std::vector<funnelEncoding>& funnelEncodings, std::vector<funnelType>& allFunnels, int nFounders);
#endif