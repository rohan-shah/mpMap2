#include "order.h"
SEXP order(SEXP mpcrossLG_sexp, SEXP groupsToOrder_sexp)
{
	Rcpp::S4 mpcrossLG;
	try
	{
		mpcrossLG = mpcrossLG_sexp;
	}
	catch(...)
	{
		throw std::runtime_error("Input mpcrossLG must be an S4 object");
	}

	Rcpp::S4 lg;
	try
	{
		lg = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("lg"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg must be an S4 object");
	}

	Rcpp::S4 rf;
	try
	{
		rf = Rcpp::as<Rcpp::S4>(mpcrossLG.slot("rf"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf must be an S4 object");
	}

	Rcpp::S4 theta;
	try
	{
		theta = Rcpp::as<Rcpp::S4>(rf.slot("theta"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta must be an S4 object");
	}

	Rcpp::RawVector originalRawData;
	try
	{
		originalRawData = theta.slot("data");
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@data must be a raw vector");
	}

	Rcpp::NumericVector levels;
	try
	{
		levels = Rcpp::as<Rcpp::NumericVector>(theta.slot("levels"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@rf@theta@levels must be a numeric vector");
	}

	std::vector<int> groups;
	try
	{
		groups = Rcpp::as<std::vector<int> >(lg.slot("groups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg@groups must be an integer vector");
	}

	std::vector<int> allGroups;
	try
	{
		allGroups = Rcpp::as<std::vector<int> >(lg.slot("allGroups"));
	}
	catch(...)
	{
		throw std::runtime_error("Slot mpcrossLG@lg@allGroups must be an integer vector");
	}

	std::vector<int> groupsToOrder;
	try
	{
		groupsToOrder = Rcpp::as<std::vector<int> >(groupsToOrder_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input groupsToOrder must be an integerVector");
	}

	//Check that groupsToOrder contains unique values, which are contained in allGroups
	std::sort(groupsToOrder.begin(), groupsToOrder.end());
	std::sort(allGroups.begin(), allGroups.end());
	if(std::unique(groupsToOrder.begin(), groupsToOrder.end()) != groupsToOrder.end())
	{
		throw std::runtime_error("Input groupsToOrder must contain unique values");
	}
	for(std::vector<int>::iterator groupToOrder = groupsToOrder.begin(); groupToOrder != groupsToOrder.end(); groupToOrder++)
	{
		std::vector<int>::iterator bound = std::lower_bound(allGroups.begin(), allGroups.end(), *groupToOrder);
		if(bound == allGroups.end() || *bound != *groupToOrder)
		{
			throw std::runtime_error("Input groupsToOrder contained a group that was not present in slot mpcrossLG@lg@allGroups");
		}
	}


	Rcpp::Function asDist("as.dist"), seriate("seriate"), getOrder("get_order"), criterion("criterion");
	int nMarkers = groups.size();

	std::vector<int> permutation;
	permutation.reserve(nMarkers);
	std::vector<int> markersThisGroup;
	for(std::vector<int>::iterator currentGroup = allGroups.begin(); currentGroup != allGroups.end(); currentGroup++)
	{
		markersThisGroup.clear();
		//linkage groups are not required to be contiguous in this case, so we have to scan through to find the number of markers in this group
		for(std::vector<int>::iterator currentMarkerGroup = groups.begin(); currentMarkerGroup != groups.end(); currentMarkerGroup++)
		{
			if(*currentMarkerGroup == *currentGroup)
			{
				markersThisGroup.push_back((int)std::distance(groups.begin(), currentMarkerGroup));
			}
		}
		int nMarkersCurrentGroup = (int)markersThisGroup.size();
		if(nMarkersCurrentGroup == 0) continue;
		Rcpp::NumericMatrix subMatrix(nMarkersCurrentGroup, nMarkersCurrentGroup);

		
		double* mem = &(subMatrix(0,0));
		//Determine whether or not the submatrix is zero, as seriation seems to screw up for all-zero matrices.
		bool nonZero = false;
		for(int i = 0; i < nMarkersCurrentGroup; i++)
		{
			for(int j = 0; j < nMarkersCurrentGroup; j++)
			{
				double* dest = mem + i + j * nMarkersCurrentGroup;
				int row = markersThisGroup[i], column = markersThisGroup[j];
				if(row > column) std::swap(row, column);
				*dest = levels[originalRawData[(column * (column + 1))/2 + row]];
				nonZero |= (*dest != 0);
			}
		}
		if(nonZero)
		{
			Rcpp::RObject distSubMatrix = asDist(Rcpp::Named("m") = subMatrix);

			std::string methods[6] = {"TSP", "OLO", "ARSA", "MDS", "GW", "HC"};
			Rcpp::List nRepsControl = Rcpp::List::create(Rcpp::Named("nreps") = 5);
			std::vector<std::vector<int> > methodResults;
			std::vector<double> objectiveFunctionValues;

			for(int methodIndex = 0; methodIndex < 6; methodIndex++)
			{
				Rcpp::RObject result;
				if(methods[methodIndex] == "ARSA")
				{
					result = seriate(Rcpp::Named("x") = distSubMatrix, Rcpp::Named("method") = methods[methodIndex], Rcpp::Named("control") = nRepsControl);
				}
				else result = seriate(Rcpp::Named("x") = distSubMatrix, Rcpp::Named("method") = methods[methodIndex]);
				std::vector<int> currentMethodResult = Rcpp::as<std::vector<int> >(getOrder(result));
				methodResults.push_back(currentMethodResult);

				Rcpp::NumericVector criterionResult = criterion(Rcpp::Named("x") = distSubMatrix, Rcpp::Named("order") = result, Rcpp::Named("method") = "AR_events");
				objectiveFunctionValues.push_back(Rcpp::as<double>(criterionResult));
			}

			//now consider identity permutation
			std::vector<int> identity; 
			for(int i = 0; i < nMarkersCurrentGroup; i++) identity.push_back(i+1);
			methodResults.push_back(identity);
			objectiveFunctionValues.push_back(Rcpp::as<double>(criterion(Rcpp::Named("x") = distSubMatrix, Rcpp::Named("method") = "AR_events")));

			std::vector<int>& currentGroupPermutation = methodResults[std::distance(objectiveFunctionValues.begin(), std::min_element(objectiveFunctionValues.begin(), objectiveFunctionValues.end()))];
			//Permutation stuff in R is indexed at base 1, but we want it indexed starting at 0.
			for(int i = 0; i < nMarkersCurrentGroup; i++) permutation.push_back(markersThisGroup[currentGroupPermutation[i]-1]+1);
		}
		else
		{
			for(int i = 0; i < nMarkersCurrentGroup; i++) permutation.push_back(markersThisGroup[i]+1);
		}
	}
	return Rcpp::wrap(permutation);
}
