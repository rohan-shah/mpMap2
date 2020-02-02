#include "parsePurdy.h"
#include <set>
bool split(std::string input, Rcpp::Function& gregexpr, std::string& output1, std::string& output2, int& level)
{
	if(input.find('*') != std::string::npos)
	{
		throw std::runtime_error("Cannot handle formats containing character \"*\"");
	}
	Rcpp::IntegerVector result = Rcpp::as<Rcpp::IntegerVector>(Rcpp::as<Rcpp::List>(gregexpr("/[0-9]+/", input))(0));
	Rcpp::IntegerVector matchLength = Rcpp::as<Rcpp::IntegerVector>(result.attr("match.length"));
	int size = -1;
	std::size_t matchPosition = -1;
	if(result.size() == 1 && result[0] == -1)
	{
		if((matchPosition = input.find("///")) != std::string::npos)
		{
			matchLength = 3;
			level = size = 3;
		}
		else if((matchPosition = input.find("//")) != std::string::npos)
		{
			matchLength = 2;
			level = size = 2;
		}
		else if((matchPosition = input.find("/")) != std::string::npos)
		{
			matchLength = 1;
			level = size = 1;
		}
		else return false;
		matchPosition++;
	}
	else if(result.size() == 1)
	{
		size = matchLength[0];
		matchPosition = result[0];
		try
		{
			std::string matched = input.substr(matchPosition, size-2);
			level = std::stoi(matched);
		}
		catch(...)
		{
			throw std::runtime_error("Unable to convert value separated by slashes into a number (E.g /3/, /4/)");
		}
	}
	else if(result.size() > 1)
	{
		level = -1;
		for(int i = 0; i < (int)result.size(); i++)
		{
			int currentLevel = std::stoi(input.substr(result[i]+1, matchLength[i]-1));
			if(currentLevel > level)
			{
				matchPosition = result[i];
				size = matchLength[i];
				level = currentLevel;
			}
		}
	}
	output1 = input.substr(0, matchPosition-1);
	output2 = input.substr(matchPosition-1 + size);
	return true;
}
SEXP parsePurdy(SEXP names_sexp, SEXP format_sexp)
{
BEGIN_RCPP
	Rcpp::CharacterVector names;
	try
	{
		names = Rcpp::as<Rcpp::CharacterVector>(names_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input names must be a character vector");
	}

	Rcpp::CharacterVector format;
	try
	{
		format = Rcpp::as<Rcpp::CharacterVector>(format_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input format must be a character vector");
	}
	std::map<std::string, int> lineIDs;
	std::vector<std::pair<int, int> > parentIDs;
	std::vector<int> levels;
	int nLines = (int)format.size();

	std::list<std::string> remainingParts;
	for(int i = 0; i < nLines; i++)
	{
		remainingParts.push_back(Rcpp::as<std::string>(format[i]));
	}

	Rcpp::Function gregexpr("gregexpr");
	while(remainingParts.size() > 0)
	{
		std::string currentLine = remainingParts.front();
		std::map<std::string, int>::iterator lookupID;
		int currentLineID;
		if((lookupID = lineIDs.find(currentLine)) == lineIDs.end())
		{
			currentLineID = (int)lineIDs.size();
			lineIDs.insert(std::make_pair(currentLine, currentLineID));
			parentIDs.push_back(std::make_pair(-2, -2));
			levels.push_back(-1);
		}
		else currentLineID = lookupID->second;
		std::string part1, part2;
		int level;
		if(split(currentLine, gregexpr, part1, part2, level))
		{
			std::map<std::string, int>::iterator parentID;
			int part1ID, part2ID;
			if((parentID = lineIDs.find(part1)) == lineIDs.end())
			{
				part1ID = (int)lineIDs.size();
				lineIDs.insert(std::make_pair(part1, part1ID));
				parentIDs.push_back(std::make_pair(-2, -2));
				remainingParts.push_back(part1);
				levels.push_back(-1);
			}
			else part1ID = parentID->second;
			if((parentID = lineIDs.find(part2)) == lineIDs.end())
			{
				part2ID = (int)lineIDs.size();
				lineIDs.insert(std::make_pair(part2, part2ID));
				parentIDs.push_back(std::make_pair(-2, -2));
				remainingParts.push_back(part2);
				levels.push_back(-1);
			}
			else part2ID = parentID->second;
			parentIDs[currentLineID] = std::make_pair(part1ID, part2ID);
			levels[currentLineID] = level;
		}
		else
		{
			parentIDs[currentLineID] = std::make_pair(-1, -1);
			levels[currentLineID] = 0;
		}
		remainingParts.pop_front();
	}
	Rcpp::CharacterMatrix results(lineIDs.size() + nLines, 3), resultsReordered(lineIDs.size() + nLines, 3);
	for(std::map<std::string, int>::iterator i = lineIDs.begin(); i != lineIDs.end(); i++)
	{
		results(i->second, 0) = i->first;
	}
	for(int i = 0; i < (int)parentIDs.size(); i++)
	{
		if(parentIDs[i].first != -1 && parentIDs[i].second != -1)
		{
			results(i, 1) = results(parentIDs[i].first, 0);
			results(i, 2) = results(parentIDs[i].second, 0);
		}
		else results(i, 1) = results(i, 2) = "";
	}
	//Also add renamed versions, with the names given as inputs. 
	for(int i = 0; i < nLines; i++)
	{
		std::string encoding = Rcpp::as<std::string>(format[i]);
		int formatLine = lineIDs[encoding];
		if(Rcpp::as<std::string>(results(formatLine, 0)) != encoding) throw std::runtime_error("Internal error");
		results(lineIDs.size() + i, 0) = Rcpp::as<std::string>(names[i]);
		results(lineIDs.size() + i, 1) = results(formatLine, 1);
		results(lineIDs.size() + i, 2) = results(formatLine, 2);
		std::pair<int, int> existingPair = parentIDs[formatLine];
		parentIDs.push_back(existingPair);
		levels.push_back(levels[formatLine]);
	}
	bool continuing = true;
	int currentLevel = 0;
	int currentOutput = 0;
	while(continuing)
	{
		continuing = false;
		for(std::size_t i = 0 ; i < parentIDs.size(); i++)
		{
			if(levels[i] == currentLevel)
			{
				continuing = true;
				resultsReordered(currentOutput, 0) = results(i, 0);
				resultsReordered(currentOutput, 1) = results(i, 1);
				resultsReordered(currentOutput, 2) = results(i, 2);
				currentOutput++;
			}
		}
		currentLevel++;
	}
	return resultsReordered;
END_RCPP
}

