	void applyFunnel(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(errorProb == 0)
		{
			applyFunnelNoError(startPosition, endPosition, finalCounter, funnelID, selfingGenerations);
		}
		else 
		{
			applyFunnelWithError(startPosition, endPosition, finalCounter, funnelID, selfingGenerations);
		}
	}
	void applyFunnelNoError(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(logFunnelSingleLociHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		const double log2 = log(2.0);
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);

		//Initialise the algorithm
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		
		//Some values are never touched, so just mark those as negative infinity
		std::fill(pathLengths1.begin(), pathLengths1.end(), -std::numeric_limits<double>::infinity());
		std::fill(pathLengths2.begin(), pathLengths2.end(), -std::numeric_limits<double>::infinity());

		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			//If so things are a bit more simple. 
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = multiple + (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
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
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths1[encodingTheseFounders] = multiple + (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					if(markerEncodingTheseFounders == startMarkerValue)
					{}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) == recodedFounders(funnel[founderCounter], startMarkerIndex))
					{
						if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) != recodedFounders(funnel[founderCounter], startMarkerIndex))
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
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						//Founder at the previous marker. 
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
								double singleLocProb = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity() && singleLocProb != -std::numeric_limits<double>::infinity())
								{
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - singleLocProb;
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
						//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
						intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
						pathLengths2[encodingTheseFounders] = *longest;
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter,markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the next marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingMarker = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						//Account for the fact that each heterozygote is only counted once, so the probabilities are half what they should really be. 
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						if(markerValue == encodingMarker)
						{}
						else if(markerValue == NA_INTEGER && founderCounter == founderCounter2)
						{
							if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else if(markerValue == NA_INTEGER && founderCounter != founderCounter2)
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
									int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
									double singleLocProb = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
									if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity() && singleLocProb != -std::numeric_limits<double>::infinity())
									{
										working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - singleLocProb;
									}
								}
							}
							//Get the longest one, and check that it's not negative infinity.
							std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
							//It's not an error for the longest path to be -Inf, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
							//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
							int bestPrevious = (int)std::distance(working.begin(), longest);
							
							memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
							intermediate2(encodingTheseFounders,positionCounter-startPosition+1) = encodingTheseFounders;
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
			while(identicalIndex != positionCounter-startPosition + 1)
			{
				int value = intermediate1(key(funnel[0], funnel[0])-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
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
	void applyFunnelWithError(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(logFunnelSingleLociHaplotypeProbabilities == NULL || logFunnelHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		const double log2 = log(2.0);

		//Initialise the algorithm
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}

		int startMarkerIndex = allPositions.markerIndices[startPosition];
		double logHomozygoteMissingProb = log(homozygoteMissingProb);
		double logHeterozygoteMissingProb = log(heterozygoteMissingProb);

		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				//We don't need to consider the other triangle of the matrix startMarkerData.hetData(founderCounter, founderCounter2)
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths2[encodingTheseFounders] = pathLengths1[encodingTheseFounders] = multiple + (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
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
				//We don't need to consider the other triangle of the matrix startMarkerData.hetData(founderCounter, founderCounter2)
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					intermediate1(encodingTheseFounders, 0) = encodingTheseFounders;
					double multiple = 0;
					if(founderCounter != founderCounter2) multiple = log2;
					pathLengths1[encodingTheseFounders] = multiple + (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					bool isError;
					if(markerEncodingTheseFounders == startMarkerValue)
					{
						pathLengths1[encodingTheseFounders] += errorTermStart1;
						isError = false;
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) == recodedFounders(funnel[founderCounter], startMarkerIndex))
					{
						if(homozygoteMissingProb != 0) pathLengths1[encodingTheseFounders] += logHomozygoteMissingProb;
						else pathLengths1[encodingTheseFounders] = -std::numeric_limits<double>::infinity();
						isError = false;
					}
					else if(startMarkerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) != recodedFounders(funnel[founderCounter], startMarkerIndex))
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
					error1(encodingTheseFounders, 0) = isError;
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
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						std::fill(working.begin(), working.end(), -std::numeric_limits<double>::infinity());
						//Founder at the previous marker. 
						for(int founderPreviousCounter = 0; founderPreviousCounter < nFounders; founderPreviousCounter++)
						{
							for(int founderPreviousCounter2 = 0; founderPreviousCounter2 <= founderPreviousCounter; founderPreviousCounter2++)
							{
								int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
								double singleLocProb = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
								if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity() && singleLocProb != -std::numeric_limits<double>::infinity())
								{
									working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - singleLocProb;
								}
							}
						}
						//Get the shortest one, and check that it's not negative infinity.
						std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
						//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
						//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
						int bestPrevious = (int)std::distance(working.begin(), longest);
						
						memcpy(&(intermediate2(encodingTheseFounders, identicalIndex)), &(intermediate1(bestPrevious, identicalIndex)), sizeof(int)*(positionCounter - startPosition + 1 - identicalIndex));
						std::copy(error1.iterator(bestPrevious, identicalIndex), error1.iterator(bestPrevious, positionCounter - startPosition + 1), error2.iterator(encodingTheseFounders, identicalIndex));
						intermediate2(encodingTheseFounders, positionCounter-startPosition+1) = encodingTheseFounders;
						error2(encodingTheseFounders, positionCounter-startPosition+1) = false;
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
						int encodingMarker = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						double multipleNextMarker = 0;
						if(founderCounter != founderCounter2) multipleNextMarker = log2;
						bool isError = false;
						if(encodingMarker == markerValue)
						{
							multipleNextMarker += errorTermCurrentMarker1;
						}
						else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) == recodedFounders(funnel[founderCounter], markerIndex))
						{
							if(homozygoteMissingProb != 0) multipleNextMarker += logHomozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
						}
						else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) != recodedFounders(funnel[founderCounter], markerIndex))
						{
							if(heterozygoteMissingProb != 0) multipleNextMarker += logHeterozygoteMissingProb;
							else multipleNextMarker = -std::numeric_limits<double>::infinity();
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
									int encodingPreviousTheseFounders = key(funnel[founderPreviousCounter], funnel[founderPreviousCounter2])-1;
									double singleLocProb = (*logFunnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderPreviousCounter][founderPreviousCounter2];
									if(pathLengths1[encodingPreviousTheseFounders] != -std::numeric_limits<double>::infinity() && singleLocProb != -std::numeric_limits<double>::infinity())
									{
										working[encodingPreviousTheseFounders] = pathLengths1[encodingPreviousTheseFounders] + multipleNextMarker + (*logFunnelHaplotypeProbabilities)(positionCounter-startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderPreviousCounter][founderPreviousCounter2] - singleLocProb;
									}
								}
							}
							//Get the shortest one, and check that it's not negative infinity.
							std::vector<double>::iterator longest = std::max_element(working.begin(), working.end());
							//This error is no longer valid, because some states are impossible to ever be in - E.g. heterozygote {1,2} with funnel {1,2,3,4} and no intercrossing. In this case all the probabilities are zero and all the log probabilities are -inf. So *longest == -std::numeric_limits<double>::infinity() doesn't indicate that there is no valid next state. It indicates that the state for the previous marker is impossible. 
							//if(*longest == -std::numeric_limits<double>::infinity()) throw std::runtime_error("Internal error");
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
				int value = intermediate1(key(funnel[0], funnel[0])-1, identicalIndex);
				for(int founderCounter = 1; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						if(value != intermediate1(encodingTheseFounders, identicalIndex)) goto stopIdenticalSearch;
					}
				}
				//We don't care about the correct indexing here. Put the correct value in every row. 
				for(int i = 0; i < (nFounders*(nFounders+1))/2; i++)
				{
					intermediate2(i, identicalIndex) = value;
					error2(i, identicalIndex) = error1(key(funnel[0], funnel[0])-1, identicalIndex);
				}
				identicalIndex++;
			}
stopIdenticalSearch:
			;
		}
	}

