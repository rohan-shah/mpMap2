	void applyIntercrossing(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(errorProb == 0)
		{
			applyIntercrossingNoError(startPosition, endPosition, finalCounter, intercrossingGeneration, selfingGenerations);
		}
		else
		{
			applyIntercrossingWithError(startPosition, endPosition, finalCounter, intercrossingGeneration, selfingGenerations);
		}
	}
	void applyIntercrossingNoError(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(logIntercrossingSingleLociHaplotypeProbabilities == NULL || logIntercrossingHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		const double log2 = log(2.0);
	
		//Initialise the algorithmi
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());

		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = multiple + (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				}
			}
		}
		else
		{
			int startMarkerValue = recodedFinals(finalCounter, startMarkerIndex);
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths1[encodingTheseFounders] = multiple + (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					if(markerEncodingTheseFounders == startMarkerValue)
					{}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, startMarkerIndex) == recodedFounders(founderCounter, startMarkerIndex))
					{
						if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, startMarkerIndex) != recodedFounders(founderCounter, startMarkerIndex))
					{
						if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
					else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
				}
			}
		}
		int identicalIndex = 0;
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			//Is the current marker just a placeholder, a location where there isn't any data?
			if(markerIndex == -1)
			{
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logIntercrossingHaplotypeProbabilities)(positionCounter-startPosition, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
						intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
						pathLengths2[encodingTheseFounders] = *longest;
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingMarker = currentMarkerData.hetData(founderCounter, founderCounter2);
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						if(encodingMarker == markerValue)
						{}
						else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerIndex) == recodedFounders(founderCounter, markerIndex))
						{
							if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, markerIndex) != recodedFounders(founderCounter, markerIndex))
						{
							if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else multipleNextMarker = -std::numeric_limits<double>::infinity();
						if(multipleNextMarker != -std::numeric_limits<double>::infinity())
						{
							std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
							//Founder at the previous marker. 
							for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
							{
								for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
								{
									int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
									if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
									{
										working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logIntercrossingHaplotypeProbabilities)(positionCounter-startPosition, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
									}
								}
							}
							//Get the shortest one, and check that it's not negative infinity.
							std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
							int bestPrevious = (int)std::distance(working.begin(), longest);
							
							memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
							intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
							pathLengths2[encodingTheseFounders] = *longest;
						}
						else
						{
							pathLengths2[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						}
					}
				}
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(positionCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			while(identicalIndex != positionCounter - startPosition + 1)
			{
				int value = intermediate1(key(0,0)-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
	void applyIntercrossingWithError(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration, int selfingGenerations)
	{
		if(logIntercrossingSingleLociHaplotypeProbabilities == NULL || logIntercrossingHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);
		const double log2 = log(2.0);

		//Initialise the algorithm
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());

		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = multiple + (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
				}
			}
		}
		else
		{
			int startMarkerValue = recodedFinals(finalCounter, startMarkerIndex);
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			
			double errorTermStart1 = log((1 - errorProb) + errorProb * 1.0 / (double) startMarkerData.nObservedValues);
			double errorTermStart2 = log(errorProb * 1.0 / (double) startMarkerData.nObservedValues);
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(founderCounter, founderCounter2);
					int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths1[encodingTheseFounders] = multiple + (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					bool isError;
					if(markerEncodingTheseFounders == startMarkerValue)
					{
						pathLengths1[encodingTheseFounders] += errorTermStart1;
						isError = false;
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, startMarkerIndex) == recodedFounders(founderCounter, startMarkerIndex))
					{
						if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						isError = false;
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(founderCounter2, startMarkerIndex) != recodedFounders(founderCounter, startMarkerIndex))
					{
						if(heterozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHeterozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						isError = false;
					}
					else
					{
						pathLengths1[encodingTheseFounders] += errorTermStart2;
						isError = true;
					}
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders];
					error2(encodingTheseFounders, 0) = error1(encodingTheseFounders, 0) = isError;
				}
			}
		}
		int identicalIndex = 0;
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			//Is the current marker just a placeholder, a location where there isn't any data?
			if(markerIndex == -1)
			{
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						//Founder at the previous marker. 
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
								{
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logIntercrossingHaplotypeProbabilities)(positionCounter-startPosition, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
						std::copy(error1.iterator(bestPrevious, identicalIndex), error1.iterator(bestPrevious, positionCounter - startPosition + 1), error2.iterator(encodingTheseFounders, identicalIndex));
						intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
						pathLengths2[encodingTheseFounders] = *longest;
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];

				double errorTermCurrentMarker1 = log((1 - errorProb) + errorProb * 1.0 / (double) currentMarkerData.nObservedValues);
				double errorTermCurrentMarker2 = log(errorProb * 1.0 / (double) currentMarkerData.nObservedValues);
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingMarker = currentMarkerData.hetData(founderCounter, founderCounter2);
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						bool isError;
						if(encodingMarker == markerValue)
						{
							multipleNextMarker += errorTermCurrentMarker1;
							isError = false;
						}
						else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, positionCounter) == recodedFounders(founderCounter, positionCounter))
						{
							if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
							isError = false;
						}
						else if(markerValue == NA_INTEGER && recodedFounders(founderCounter2, positionCounter) != recodedFounders(founderCounter, positionCounter))
						{
							if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
							isError = false;
						}
						else
						{
							multipleNextMarker += errorTermCurrentMarker2;
							isError = true;
						}
						if(multipleNextMarker != -std::numeric_limits<double>::infinity())
						{
							//Founder at the previous marker. 
							std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
							for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
							{
								for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
								{
									int encodingPreviousTheseFounders = key(founderPreviousCounter, founderPreviousCounter2)-1;
									if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity())
									{
										working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logIntercrossingHaplotypeProbabilities)(positionCounter-startPosition, intercrossingGeneration - minAIGenerations, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - (*logIntercrossingSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
									}
								}
							}
							//Get the shortest one, and check that it's not negative infinity.
							std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
							int bestPrevious = (int)std::distance(working.begin(), longest);
							
							memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
							std::copy(error1.iterator(bestPrevious, identicalIndex), error1.iterator(bestPrevious, positionCounter - startPosition + 1), error2.iterator(encodingTheseFounders, identicalIndex));
							intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
							error2(encodingTheseFounders, positionCounter-startPosition+1) = isError;
							pathLengths2[encodingTheseFounders] = *longest;
						}
						else
						{
							pathLengths2[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						}
					}
				}
			}
			//If this condition throws, it's almost guaranteed to be because the map contains two markers at the same location, but the data implies a non-zero distance because recombinations are observed to occur between them.
			std::vector<double>::iterator longest = std::max_element(pathLengths2.begin(), pathLengths2.end());
			if(*longest == -std::numeric_limits<double>::infinity()) throw impossibleDataException(positionCounter, finalCounter);

			intermediate1.swap(intermediate2);
			pathLengths1.swap(pathLengths2);
			error1.swap(error2);
			while(identicalIndex != positionCounter-startPosition + 1)
			{
				int value = intermediate1(key(0,0)-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(founderCounter, founderCounter2)-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
					error2(i, identicalIndex) = error1(key(0,0)-1, identicalIndex);
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}
