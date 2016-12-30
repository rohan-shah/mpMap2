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
		if(funnelSingleLociHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		//Compute forward probabilities
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerIndex = allPositions.markerIndices[startPosition];

		double sum = 0;
		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] ;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
		}
		else
		{
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			int markerValue = recodedFinals(finalCounter, startMarkerIndex);
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					if(markerValue == markerEncodingTheseFounders)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) == recodedFounders(funnel[founderCounter], startMarkerIndex) && homozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * homozygoteMissingProb;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) != recodedFounders(funnel[founderCounter], startMarkerIndex) && heterozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * heterozygoteMissingProb;
					}
					else forwardProbabilities(encodingTheseFounders, 0) = 0;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
		}
		for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
		{
			forwardProbabilities(counter, 0) /= sum;
		}
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			//Is the current marker just a placeholder, a location where there isn't any data?
			if(markerIndex == -1)
			{
				//The founders at the new marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) = 0;
						//Founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) += forwardProbabilities(encodingPreviousFounders, positionCounter - startPosition) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
						}
						sum += forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1);
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founders at the new marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int markerEncodingTheseFounders = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) = 0;
						bool markerMatches = markerValue == markerEncodingTheseFounders;
						bool missingHet = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) != recodedFounders(funnel[founderCounter], markerIndex) && heterozygoteMissingProb != 0;
						bool missingHomo = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) == recodedFounders(funnel[founderCounter], markerIndex) && homozygoteMissingProb != 0;
						if(markerMatches || missingHet || missingHomo)
						{
							double factor = 1;
							if(missingHet) factor = heterozygoteMissingProb;
							if(missingHomo) factor = homozygoteMissingProb;
							//Founders at the previous marker
							for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
							{
								for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
								{
									int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
									forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) += forwardProbabilities(encodingPreviousFounders, positionCounter - startPosition) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * factor;
								}
							}
						}
						sum += forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1);
					}
				}
			}
			for(int counter = 0; counter < nFounders*(nFounders+1)/2; counter++)
			{
				forwardProbabilities(counter, positionCounter - startPosition + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
				backwardProbabilities(encodingTheseFounders, endPosition - startPosition - 1) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
			}
		}
		for(int positionCounter = endPosition - 2; positionCounter >= startPosition; positionCounter--)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			if(markerIndex == -1)
			{
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) = 0;
						//The founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
						}
						sum += backwardProbabilities(encodingTheseFounders, positionCounter - startPosition);
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) = 0;
						//The founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int markerEncodingPreviousFounders = currentMarkerData.hetData(funnel[founderCounterPrevious], funnel[founderCounterPrevious2]);
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								if(markerValue == markerEncodingPreviousFounders)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
								}
								else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerIndex) == recodedFounders(founderCounterPrevious, markerIndex) && homozygoteMissingProb != 0)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * homozygoteMissingProb;
								}
								else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerIndex) != recodedFounders(founderCounterPrevious, markerIndex) && heterozygoteMissingProb != 0)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * heterozygoteMissingProb;
								}
							}
						}
						sum += backwardProbabilities(encodingTheseFounders, positionCounter - startPosition);
					}
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				backwardProbabilities(counter, positionCounter - startPosition) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int positionCounter = startPosition; positionCounter < endPosition; positionCounter++)
		{
			double sum = 0;
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter) = backwardProbabilities(counter, positionCounter - startPosition) * forwardProbabilities(counter, positionCounter - startPosition);
				sum += results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter);
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter) /= sum;
			}
		}
	}
	void applyFunnelWithError(int startPosition, int endPosition, int finalCounter, int funnelID, int selfingGenerations)
	{
		if(funnelSingleLociHaplotypeProbabilities == NULL || errorProb != errorProb || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		//Compute forward probabilities
		funnelEncoding enc = (*lineFunnelEncodings)[(*lineFunnelIDs)[finalCounter]];
		int funnel[16];
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			funnel[founderCounter] = ((enc & (15 << (4*founderCounter))) >> (4*founderCounter));
		}
		int startMarkerIndex = allPositions.markerIndices[startPosition];

		double sum = 0;
		//Is the first marker just a placeholder - A place where we want to impute, but for which there is no actual data?
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] ;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
		}
		else
		{
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			int markerValue = recodedFinals(finalCounter, startMarkerIndex);

			double errorTermStart1 = (1 - errorProb) + errorProb * 1.0 / (double) startMarkerData.nObservedValues;
			double errorTermStart2 = errorProb * 1.0 / (double) startMarkerData.nObservedValues;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
				{
					int markerEncodingTheseFounders = startMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
					int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
					if(markerValue == markerEncodingTheseFounders)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * errorTermStart1;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) == recodedFounders(funnel[founderCounter], startMarkerIndex) && homozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * homozygoteMissingProb;
					}
					else if(markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], startMarkerIndex) != recodedFounders(funnel[founderCounter], startMarkerIndex) && heterozygoteMissingProb != 0)
					{
						forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * heterozygoteMissingProb;
					}
					else forwardProbabilities(encodingTheseFounders, 0) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2] * errorTermStart2;
					sum += forwardProbabilities(encodingTheseFounders, 0);
				}
			}
		}
		for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
		{
			forwardProbabilities(counter, 0) /= sum;
		}
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			//Is the current marker just a placeholder, a location where there isn't any data?
			if(markerIndex == -1)
			{
				//The founders at the new marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) = 0;
						//Founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) += forwardProbabilities(encodingPreviousFounders, positionCounter - startPosition) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
						}
						sum += forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1);
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];

				double errorTermCurrentMarker1 = (1 - errorProb) + errorProb * 1.0 / (double) currentMarkerData.nObservedValues;
				double errorTermCurrentMarker2 = errorProb * 1.0 / (double) currentMarkerData.nObservedValues;
				//The founders at the new marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int markerEncodingTheseFounders = currentMarkerData.hetData(funnel[founderCounter], funnel[founderCounter2]);
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) = 0;
						bool markerMatches = markerValue == markerEncodingTheseFounders;
						bool missingHet = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) != recodedFounders(funnel[founderCounter], markerIndex) && heterozygoteMissingProb != 0;
						bool missingHomo = markerValue == NA_INTEGER && recodedFounders(funnel[founderCounter2], markerIndex) == recodedFounders(funnel[founderCounter], markerIndex) && homozygoteMissingProb != 0;
						double factor = 1;
						if(markerMatches) factor = errorTermCurrentMarker1;
						else if(missingHet) factor = heterozygoteMissingProb;
						else if(missingHomo) factor = homozygoteMissingProb;
						else factor = errorTermCurrentMarker2;
						//Founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1) += forwardProbabilities(encodingPreviousFounders, positionCounter - startPosition) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * factor;
							}
						}
						sum += forwardProbabilities(encodingTheseFounders, positionCounter - startPosition + 1);
					}
				}
			}
			for(int counter = 0; counter < nFounders*(nFounders+1)/2; counter++)
			{
				forwardProbabilities(counter, positionCounter - startPosition + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
			{
				int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
				backwardProbabilities(encodingTheseFounders, endPosition - startPosition - 1) = (*funnelSingleLociHaplotypeProbabilities)[selfingGenerations - minSelfingGenerations].values[founderCounter][founderCounter2];
			}
		}
		for(int positionCounter = endPosition - 2; positionCounter >= startPosition; positionCounter--)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			if(markerIndex == -1)
			{
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) = 0;
						//The founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2];
							}
						}
						sum += backwardProbabilities(encodingTheseFounders, positionCounter - startPosition);
					}
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
		
				double errorTermCurrentMarker1 = (1 - errorProb) + errorProb * 1.0 / (double) currentMarkerData.nObservedValues;
				double errorTermCurrentMarker2 = errorProb * 1.0 / (double) currentMarkerData.nObservedValues;
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					for(int founderCounter2 = 0; founderCounter2 <= founderCounter; founderCounter2++)
					{
						int encodingTheseFounders = key(funnel[founderCounter], funnel[founderCounter2])-1;
						backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) = 0;
						//The founders at the previous marker
						for(int founderCounterPrevious = 0; founderCounterPrevious < nFounders; founderCounterPrevious++)
						{
							for(int founderCounterPrevious2 = 0; founderCounterPrevious2 <= founderCounterPrevious; founderCounterPrevious2++)
							{
								int markerEncodingPreviousFounders = currentMarkerData.hetData(funnel[founderCounterPrevious], funnel[founderCounterPrevious2]);
								int encodingPreviousFounders = key(funnel[founderCounterPrevious], funnel[founderCounterPrevious2])-1;
								if(markerValue == markerEncodingPreviousFounders)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * errorTermCurrentMarker1;
								}
								else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerIndex) == recodedFounders(founderCounterPrevious, markerIndex) && homozygoteMissingProb != 0)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * homozygoteMissingProb;
								}
								else if(markerValue == NA_INTEGER && recodedFounders(founderCounterPrevious2, markerIndex) != recodedFounders(founderCounterPrevious, markerIndex) && heterozygoteMissingProb != 0)
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * heterozygoteMissingProb;
								}
								else
								{
									backwardProbabilities(encodingTheseFounders, positionCounter - startPosition) += backwardProbabilities(encodingPreviousFounders, positionCounter - startPosition + 1) * funnelHaplotypeProbabilities(positionCounter - startPosition, selfingGenerations - minSelfingGenerations).values[founderCounter][founderCounter2][founderCounterPrevious][founderCounterPrevious2] * errorTermCurrentMarker2;
								}
							}
						}
						sum += backwardProbabilities(encodingTheseFounders, positionCounter - startPosition);
					}
				}
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				backwardProbabilities(counter, positionCounter - startPosition) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int positionCounter = startPosition; positionCounter < endPosition; positionCounter++)
		{
			double sum = 0;
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter) = backwardProbabilities(counter, positionCounter - startPosition) * forwardProbabilities(counter, positionCounter - startPosition);
				sum += results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter);
			}
			for(int counter = 0; counter < (nFounders*(nFounders+1))/2; counter++)
			{
				results(((nFounders*(nFounders+1))/2)*finalCounter + counter, positionCounter) /= sum;
			}
		}
	}

