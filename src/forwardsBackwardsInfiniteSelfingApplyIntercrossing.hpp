	void applyIntercrossing(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration)
	{
		if(errorProb == 0)
		{
			applyIntercrossingNoError(startPosition, endPosition, finalCounter, intercrossingGeneration);
		}
		else
		{
			applyIntercrossingWithError(startPosition, endPosition, finalCounter, intercrossingGeneration);
		}
	}
	void applyIntercrossingNoError(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration)
	{
		if(errorProb != errorProb || errorProb != 0)
		{
			throw std::runtime_error("Internal error");
		}
		//Compute forward probabilities
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) = 1.0/nFounders;
			}
		}
		else
		{
			int markerValue = recodedFinals(finalCounter, startMarkerIndex);
			int validInitial = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(founderCounter, startMarkerIndex) == markerValue || markerValue == NA_INTEGER)
				{
					forwardProbabilities(founderCounter, 0) = 1;
					validInitial++;
				}
				else forwardProbabilities(founderCounter, 0) = 0;
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) /= (double)validInitial;
			}
		}
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			if(markerIndex == -1)
			{
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					forwardProbabilities(founderCounter, positionCounter - startPosition + 1) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
					sum += forwardProbabilities(founderCounter, positionCounter - startPosition + 1);
				}
			}
			else
			{
				//The founder at the new marker
				int markerValue = recodedFinals(finalCounter, markerIndex);
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					forwardProbabilities(founderCounter, positionCounter - startPosition + 1) = 0;
					if(recodedFounders(founderCounter, markerIndex) == markerValue || markerValue == NA_INTEGER)
					{
						//The founder at the previous marker
						for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
						{
							forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
						}
					}
					sum += forwardProbabilities(founderCounter, positionCounter - startPosition + 1);
				}
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, positionCounter - startPosition + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			backwardProbabilities(founderCounter, endPosition - startPosition - 1) = 1/(double)nFounders;
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
					backwardProbabilities(founderCounter, positionCounter - startPosition) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
					sum += backwardProbabilities(founderCounter, positionCounter - startPosition);
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					backwardProbabilities(founderCounter, positionCounter - startPosition) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						if(recodedFounders(founderCounter2, markerIndex) == markerValue || markerValue == NA_INTEGER)
						{
							backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
						}
					}
					sum += backwardProbabilities(founderCounter, positionCounter - startPosition);
				}
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(founderCounter, positionCounter - startPosition) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int positionCounter = startPosition; positionCounter < endPosition; positionCounter++)
		{
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(nFounders*finalCounter + founderCounter, positionCounter) = backwardProbabilities(founderCounter, positionCounter - startPosition) * forwardProbabilities(founderCounter, positionCounter - startPosition);
				sum += results(nFounders*finalCounter + founderCounter, positionCounter);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(nFounders*finalCounter + founderCounter, positionCounter) /= sum;
			}
		}
	}
	void applyIntercrossingWithError(int startPosition, int endPosition, int finalCounter, int intercrossingGeneration)
	{
		if(errorProb != errorProb || errorProb <= 0 || errorProb >= 1)
		{
			throw std::runtime_error("Internal error");
		}
		//Compute forward probabilities
		int startMarkerIndex = allPositions.markerIndices[startPosition];
		if(startMarkerIndex == -1)
		{
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) = 1.0/nFounders;
			}
		}
		else
		{
			int markerValue = recodedFinals(finalCounter, startMarkerIndex);
			::markerData& startMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[startMarkerIndex]];
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				if(recodedFounders(founderCounter, startMarkerIndex) == markerValue) forwardProbabilities(founderCounter, 0) = (1.0 / (double)nFounders) * ((1 - errorProb) + errorProb / (double)startMarkerData.nObservedValues);
				else if(markerValue == NA_INTEGER) forwardProbabilities(founderCounter, 0) = 1.0 / (double)nFounders;
				else forwardProbabilities(founderCounter, 0) = (1.0 / (double)nFounders) * errorProb / (double)startMarkerData.nObservedValues;
				sum += forwardProbabilities(founderCounter, 0);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, 0) /= sum;
			}
		}
		for(int positionCounter = startPosition; positionCounter < endPosition - 1; positionCounter++)
		{
			int markerIndex = allPositions.markerIndices[positionCounter+1];
			double sum = 0;
			if(markerIndex == -1)
			{
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					forwardProbabilities(founderCounter, positionCounter - startPosition + 1) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
					sum += forwardProbabilities(founderCounter, positionCounter - startPosition + 1);
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the new marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					forwardProbabilities(founderCounter, positionCounter - startPosition + 1) = 0;
					if(recodedFounders(founderCounter, markerIndex) == markerValue)
					{
						//The founder at the previous marker
						for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
						{
							forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter] * ((1 - errorProb) + errorProb / (double)currentMarkerData.nObservedValues);
						}
					}
					else if(markerValue == NA_INTEGER)
					{
						//The founder at the previous marker
						for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
						{
							forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
						}
					}
					else
					{
						//The founder at the previous marker
						for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
						{
							forwardProbabilities(founderCounter, positionCounter - startPosition + 1) += forwardProbabilities(founderCounter2, positionCounter - startPosition) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter] * errorProb / (double)currentMarkerData.nObservedValues;
						}
					}
					sum += forwardProbabilities(founderCounter, positionCounter - startPosition + 1);
				}
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				forwardProbabilities(founderCounter, positionCounter - startPosition + 1) /= sum;
			}
		}
		//Now the backwards probabilities
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			backwardProbabilities(founderCounter, endPosition - startPosition - 1) = 1/(double)nFounders;
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
					backwardProbabilities(founderCounter, positionCounter - startPosition) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
					}
					sum += backwardProbabilities(founderCounter, positionCounter - startPosition);
				}
			}
			else
			{
				int markerValue = recodedFinals(finalCounter, markerIndex);
				::markerData& currentMarkerData = markerData.allMarkerPatterns[markerData.markerPatternIDs[markerIndex]];
				//The founder at the current marker
				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					backwardProbabilities(founderCounter, positionCounter - startPosition) = 0;
					//The founder at the previous marker
					for(int founderCounter2 = 0; founderCounter2 < nFounders; founderCounter2++)
					{
						if(recodedFounders(founderCounter2, markerIndex) == markerValue || markerValue == NA_INTEGER)
						{
							backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter] * ((1 - errorProb) + errorProb / (double)currentMarkerData.nObservedValues);
						}
						else if(markerValue == NA_INTEGER)
						{
							backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter];
						}
						else
						{
							backwardProbabilities(founderCounter, positionCounter - startPosition) += backwardProbabilities(founderCounter2, positionCounter - startPosition + 1) * intercrossingHaplotypeProbabilities(positionCounter - startPosition, intercrossingGeneration - minAIGenerations, 0).values[founderCounter2][founderCounter] * errorProb / (double)currentMarkerData.nObservedValues;
						}
					}
					sum += backwardProbabilities(founderCounter, positionCounter - startPosition);
				}
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				backwardProbabilities(founderCounter, positionCounter - startPosition) /= sum;
			}
		}
		//Now we can compute the marginal probabilities
		for(int positionCounter = startPosition; positionCounter < endPosition; positionCounter++)
		{
			double sum = 0;
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(nFounders*finalCounter + founderCounter, positionCounter) = backwardProbabilities(founderCounter, positionCounter - startPosition) * forwardProbabilities(founderCounter, positionCounter - startPosition);
				sum += results(nFounders*finalCounter + founderCounter, positionCounter);
			}
			for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
			{
				results(nFounders*finalCounter + founderCounter, positionCounter) /= sum;
			}
		}
	}

