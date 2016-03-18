#include "preClusterStep.h"
SEXP preClusterStep(SEXP mpcrossRF_)
{
BEGIN_RCPP
	Rcpp::S4 mpcrossRF = mpcrossRF_;
	Rcpp::S4 rf = mpcrossRF.slot("rf");
	Rcpp::S4 theta = rf.slot("theta");
	Rcpp::RawVector data = theta.slot("data");
	Rcpp::CharacterVector markers = theta.slot("markers");
	Rcpp::NumericVector levels = theta.slot("levels");

	Rcpp::NumericVector::iterator zeroIterator = std::find(levels.begin(), levels.end(), 0);
	if(zeroIterator == levels.end())
	{
		throw std::runtime_error("Slot levels in mpcrossRF@rf@theta must contain the value 0");
	}
	Rbyte zeroLevel = (Rbyte)std::distance(levels.begin(), zeroIterator);
	R_xlen_t nMarkers = markers.size();
	struct group
	{
	public:
		group()
		{
			data.reserve(1);
		}
		group(group&& other)
			: data(std::move(other.data))
		{}
		group& operator=(group&& other)
		{
			data.swap(other.data);
			return *this;
		}
		std::vector<int> data;
	private:
		group(const group& other);
	};
	//The groups that are now fixed
	std::vector<group> finalisedGroups;
	//The groups that we're currently considering
	std::vector<group> continuingGroups(nMarkers);
	//The groups that we will consider in the next step
	std::vector<group> newContinuingGroups(nMarkers);

	//Initially every marker is in its own group
	for(R_xlen_t i = 0; i < nMarkers; i++) continuingGroups[i].data.push_back((int)i);
	while(continuingGroups.size() > 0)
	{
		for(std::size_t i = 0; i < continuingGroups.size(); i++)
		{
			if(continuingGroups[i].data.size() == 0) continue;
			for(std::size_t j = i + 1; j < continuingGroups.size(); j++)
			{
				const std::vector<int>& iData = continuingGroups[i].data, &jData = continuingGroups[j].data;
				if(jData.size() == 0) continue;
				//Check that every marker has recombination 0 with every marker in group j
				for(std::size_t i_ = 0; i_ < iData.size(); i_++)
				{
					R_xlen_t iMarker = iData[i_];
					for(std::size_t j_ = 0; j_ < jData.size(); j_++)
					{
						R_xlen_t jMarker = jData[j_];
						R_xlen_t rowMarker = std::min(jMarker, iMarker);
						R_xlen_t columnMarker = std::max(jMarker, iMarker);
						if(data[(columnMarker * (columnMarker+(R_xlen_t)1))/(R_xlen_t)2 + rowMarker] != zeroLevel)
						{
							//If it doesn't, continue to the next group j
							goto nextJ;
						}
					}
				}
				//If it does, combine the groups and add them to the set to be considered in the next step
				continuingGroups[i].data.insert(continuingGroups[i].data.end(), continuingGroups[j].data.begin(), continuingGroups[j].data.end());
				continuingGroups[j].data.clear();
				newContinuingGroups.emplace_back(std::move(continuingGroups[i]));
				goto nextI;
	nextJ:
				;
			}
			finalisedGroups.emplace_back(std::move(continuingGroups[i]));
	nextI:
			;
		}
		continuingGroups.swap(newContinuingGroups);
		newContinuingGroups.clear();
	}
	Rcpp::List result(finalisedGroups.size());
	for(std::size_t i = 0; i < finalisedGroups.size(); i++)
	{
		Rcpp::IntegerVector currentGroup = Rcpp::wrap(finalisedGroups[i].data);
		//Add 1, because these are going to be R indices
		for(int j = 0; j < currentGroup.size(); j++) currentGroup[j]++;
		result[i] = currentGroup;
	}
	return result;
END_RCPP
}

