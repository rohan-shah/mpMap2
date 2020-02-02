#include "joinMapWithExtra.h"
void joinMapWithExtra(Rcpp::List map, Rcpp::List extraPositions, positionData& allPositions)
{
	Rcpp::Function sort("sort");
	allPositions.names.clear();
	allPositions.positions.clear();
	allPositions.markerIndices.clear();
	allPositions.chromosomes.clear();

	std::vector<std::string> chromosomes = Rcpp::as<std::vector<std::string> >(map.names());
	std::size_t cumulativeMarkerCount = 0;
	for(std::vector<std::string>::iterator i = chromosomes.begin(); i != chromosomes.end(); i++)
	{
		std::size_t currentSize = allPositions.names.size();

		Rcpp::NumericVector currentChromosomeMap = Rcpp::as<Rcpp::NumericVector>(sort(map(*i)));
		Rcpp::CharacterVector currentChromosomeMapMarkerNames = Rcpp::as<Rcpp::CharacterVector>(currentChromosomeMap.names());
		if(extraPositions.containsElementNamed(i->c_str()))
		{
			Rcpp::NumericVector currentChromosomeExtra = Rcpp::as<Rcpp::NumericVector>(sort(extraPositions(*i)));
			Rcpp::CharacterVector currentChromosomeExtraNames = Rcpp::as<Rcpp::CharacterVector>(currentChromosomeExtra.names());
			int markerIndex = 0, extraIndex = 0;
			int totalCount = (int)(currentChromosomeExtra.size() + currentChromosomeMap.size());
			for(int j = 0; j < totalCount; j++)
			{
				//The next position has to be one of the extra ones
				if(markerIndex == (int)currentChromosomeMap.size())
				{
					allPositions.names.push_back(Rcpp::as<std::string>(currentChromosomeExtraNames[extraIndex]));
					allPositions.positions.push_back(currentChromosomeExtra[extraIndex]);
					allPositions.markerIndices.push_back(-1);
					extraIndex++;
				}
				else if(extraIndex == (int)currentChromosomeExtra.size())
				{
					allPositions.names.push_back(Rcpp::as<std::string>(currentChromosomeMapMarkerNames[markerIndex]));
					allPositions.positions.push_back(currentChromosomeMap[markerIndex]);
					allPositions.markerIndices.push_back((int)cumulativeMarkerCount + markerIndex);
					markerIndex++;
				}
				//Either is possible, so take the *next* one in terms of order.
				else if(currentChromosomeExtra[extraIndex] < currentChromosomeMap[markerIndex])
				{
					allPositions.names.push_back(Rcpp::as<std::string>(currentChromosomeExtraNames[extraIndex]));
					allPositions.positions.push_back(currentChromosomeExtra[extraIndex]);
					allPositions.markerIndices.push_back(-1);
					extraIndex++;
				}
				else
				{
					allPositions.names.push_back(Rcpp::as<std::string>(currentChromosomeMapMarkerNames[markerIndex]));
					allPositions.positions.push_back(currentChromosomeMap[markerIndex]);
					allPositions.markerIndices.push_back((int)cumulativeMarkerCount + markerIndex);
					markerIndex++;
				}
			}
			positionData::chromosomeDescriptor currentChromosome;
			currentChromosome.start = (int)currentSize;
			currentChromosome.end = (int)currentSize + (int)(currentChromosomeMap.size() + currentChromosomeExtra.size());
			currentChromosome.name = *i;
			allPositions.chromosomes.push_back(currentChromosome);
		}
		//No extra positions for this chromosome
		else
		{
			std::transform(currentChromosomeMapMarkerNames.begin(), currentChromosomeMapMarkerNames.end(), std::back_inserter(allPositions.names), [](Rcpp::String::StringProxy& x){return Rcpp::as<std::string>(x);});
			allPositions.positions.insert(allPositions.positions.end(), currentChromosomeMap.begin(), currentChromosomeMap.end());
			for(int j = 0; j < (int)currentChromosomeMap.size(); j++) allPositions.markerIndices.push_back((int)cumulativeMarkerCount+j);

			positionData::chromosomeDescriptor currentChromosome;
			currentChromosome.start = (int)currentSize;
			currentChromosome.end = (int)currentSize + (int)currentChromosomeMap.size();
			currentChromosome.name = *i;
			allPositions.chromosomes.push_back(currentChromosome);
		}
		cumulativeMarkerCount += currentChromosomeMap.size();
	}

	//Check marker names
	std::vector<std::string> allMarkerNames = allPositions.names;
	std::sort(allMarkerNames.begin(), allMarkerNames.end());
	if(std::unique(allMarkerNames.begin(), allMarkerNames.end()) != allMarkerNames.end())
	{
		throw std::runtime_error("Extra positions cannot have the same name as a marker");
	}
}
Rcpp::List positionData::makeUnifiedMap()
{
	//Put together a unified map of all the markers and extra locations
	Rcpp::List unifiedMap(chromosomes.size());
	for(int i = 0; i < (int)chromosomes.size(); i++)
	{
		const positionData::chromosomeDescriptor& currentChromosome = chromosomes[i];
		Rcpp::NumericVector unifiedMapCurrentChromosome(currentChromosome.end - currentChromosome.start);
		Rcpp::CharacterVector unifiedMapCurrentChromosomeNames(currentChromosome.end - currentChromosome.start);
		for(int j = 0; j < currentChromosome.end - currentChromosome.start; j++)
		{
			unifiedMapCurrentChromosome[j] = positions[j + currentChromosome.start];
			unifiedMapCurrentChromosomeNames[j] = names[j + currentChromosome.start];
		}
		unifiedMapCurrentChromosome.names() = unifiedMapCurrentChromosomeNames;
		unifiedMap[i] = unifiedMapCurrentChromosome;
	}
	return unifiedMap;
}
