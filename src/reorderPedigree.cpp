#include "reorderPedigree.h"
#ifdef USE_BOOST
#include <boost/graph/adjacency_list.h>
#include <boost/graph/topological_sort.h>
SEXP reorderPedigree(SEXP mother_sexp, SEXP father_sexp, SEXP lineNames_sexp)
{
BEGIN_RCPP
	Rcpp::IntegerVector mother = Rcpp::as<Rcpp::IntegerVector>(mother_sexp);
	Rcpp::IntegerVector father = Rcpp::as<Rcpp::IntegerVector>(father_sexp);
	Rcpp::CharacterVector lineNames = Rcpp::as<Rcpp::CharacterVector>(lineNames_sexp);
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property> graphType;
	int nVertices = mother.size();
	if(father.size() != nVertices)
	{
		throw std::runtime_error("Slots mother and father had different sizes");
	}
	graphType graph(nVertices);
	std::vector<int> founderIndices;
	for(int i = 0; i < nVertices; i++)
	{
		if(mother[i] == NA_INTEGER || father[i] == NA_INTEGER)
		{
			throw std::runtime_error("Values of NA are not allowed in the pedigree");
		}
		if(mother[i] < 0 || father[i] < 0) 
		{
			throw std::runtime_error("Negative values are not allowed in the pedigree");
		}
		if(mother[i] != 0)
		{
			boost::add_edge(i, mother[i]-1, graph);
		}
		if(father[i] != 0)
		{
			boost::add_edge(i, father[i]-1, graph);
		}
		if(mother[i] == 0 && father[i] == 0) founderIndices.push_back(i);
	}
	std::vector<int> order(nVertices);
	boost::topological_sort(graph, order.begin());
	//Ensure that the lines with mother == 0 and father == 0 come first in the ordering
	std::vector<int> newOrder;
	for(int i = 0; i < founderIndices.size(); i++)
	{
		newOrder.push_back(founderIndices[i]);
		for(int j = 0; j < order.size(); j++)
		{
			if(order[j] == founderIndices[i])
			{
				order[j] = -1; break;
			}
		}
	}
	for(int j = 0; j < order.size(); j++)
	{
		if(order[j] != -1) newOrder.push_back(order[j]);
	}
	order.swap(newOrder);

	std::vector<int> orderInverse(nVertices);
	for(int i = 0; i < nVertices; i++)
	{
		orderInverse[order[i]] = i;
	}

	Rcpp::IntegerVector newMother(nVertices), newFather(nVertices);
	Rcpp::CharacterVector newLineNames(nVertices);
	for(int i = 0; i < nVertices; i++)
	{
		if(mother[order[i]] == 0)
		{
			newMother[i] = 0;
		}
		else
		{
			newMother[i] = orderInverse[mother[order[i]]-1]+1;
		}
		if(father[order[i]] == 0)
		{
			newFather[i] = 0;
		}
		else
		{
			newFather[i] = orderInverse[father[order[i]]-1]+1;
		}
		newLineNames[i] = lineNames[order[i]];
	}
	return Rcpp::List::create(Rcpp::Named("mother") = newMother, Rcpp::Named("father") = newFather, Rcpp::Named("lineNames") = newLineNames);
END_RCPP
}
#endif
