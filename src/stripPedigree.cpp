#ifdef USE_BOOST
#include "stripPedigree.h"
#include <boost/graph/adjacency_list.h>
#include <boost/graph/depth_first_search.h>
#include <boost/graph/topological_sort.h>
#include <boost/graph/reverse_graph.h>
enum pedigreeEdgeType
{
	MOTHER, FATHER
};
bool lineNamesSorter(const std::pair<std::string, int>& first, const std::pair<std::string, int>& second)
{
	return first.first < second.first;
}
//Here vertex_descriptor will be the original vertex number in the original pedigree, and vertex_name will be the vertex number in the new pedigree.
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::property<boost::vertex_name_t, std::size_t>, boost::property<boost::edge_name_t, pedigreeEdgeType> > graphType;
typedef boost::color_traits<boost::default_color_type> Colour;
//Visitor used to work out the number of generations for each line in the pedigree. 
class generationsVisitor : public boost::default_dfs_visitor
{
public:
	generationsVisitor(Rcpp::IntegerVector& generations)
		: generations(generations), currentGeneration(-1)
	{}
	template<typename Edge, typename Graph> void forward_or_cross_edge(Edge e, Graph& g)
	{
		int& targetGeneration = generations[boost::get(boost::vertex_name, g, boost::target(e, g))];
		if(targetGeneration != -2 && targetGeneration != currentGeneration+1)
		{
			targetGeneration = -1;
		}
	}
	template<typename Vertex, typename Graph> void discover_vertex(Vertex v, Graph& g)
	{
		currentGeneration++;
		generations[boost::get(boost::vertex_name, g, v)] = currentGeneration;
	}
	template<typename Vertex, typename Graph> void finish_vertex(Vertex v, Graph& g)
	{
		currentGeneration--;
	}
	int currentGeneration;
	Rcpp::IntegerVector& generations;
};
SEXP stripPedigree(SEXP pedigree_sexp, SEXP finalLines_sexp)
{
BEGIN_RCPP
	Rcpp::S4 pedigree = Rcpp::as<Rcpp::S4>(pedigree_sexp);
	Rcpp::IntegerVector mother = pedigree.slot("mother");
	Rcpp::IntegerVector father = pedigree.slot("father");
	Rcpp::CharacterVector lineNames = pedigree.slot("lineNames");

	std::vector<std::pair<std::string, int> > sortedLineNames;
	graphType graph(mother.size());
	int nFounders = 0;
	for(std::size_t pedigreeRow = 0; pedigreeRow < lineNames.size(); pedigreeRow++)
	{
		if(mother[pedigreeRow] != 0)
		{
			boost::add_edge(pedigreeRow, mother[pedigreeRow]-1, MOTHER, graph);
		}
		if(father[pedigreeRow] != 0)
		{
			boost::add_edge(pedigreeRow, father[pedigreeRow]-1, FATHER, graph);
		}
		if(mother[pedigreeRow] == 0 && father[pedigreeRow] == 0) nFounders++;
		if((mother[pedigreeRow] == 0) ^ (father[pedigreeRow] == 0)) throw std::runtime_error("Internal error");

		sortedLineNames.push_back(std::make_pair(Rcpp::as<std::string>(lineNames[pedigreeRow]), pedigreeRow));
	}
	std::sort(sortedLineNames.begin(), sortedLineNames.end(), lineNamesSorter);

	Rcpp::CharacterVector finalLines = Rcpp::as<Rcpp::CharacterVector>(finalLines_sexp);
	//Run a depth first search for every entry in finalLines. Use sortedLineNames to look up the vertex to start at
	std::vector<boost::default_color_type> colourVector(mother.size(), Colour::white());
	for(std::size_t finalCounter = 0; finalCounter < finalLines.size(); finalCounter++)
	{
		std::vector<std::pair<std::string, int> >::iterator i = std::lower_bound(sortedLineNames.begin(), sortedLineNames.end(), std::make_pair(Rcpp::as<std::string>(finalLines[finalCounter]), 0), lineNamesSorter);
		if(i == sortedLineNames.end()) throw std::runtime_error("Internal error");
		boost::depth_first_visit(graph, i->second, boost::make_dfs_visitor(boost::null_visitor()), &(colourVector[0]));
	}

	//We want to do a topological sort on the marked subset. So swap white and black in the colour vector. 
	for(std::vector<boost::default_color_type>::iterator i = colourVector.begin(); i != colourVector.end(); i++)
	{
		if(*i == Colour::black()) *i = Colour::white();
		else *i = Colour::black();
	}
	//A copy of the colour vector, which we can use later to get out the generation for every vertex. We can't combine this graph search with the topological sort one, because that one starts from the end of the pedigree, and we want to start from the founders (the beginning)
	std::vector<boost::default_color_type> copiedColourVector = colourVector;

	std::vector<std::size_t> topologicalSortOutput;
	graphType::vertex_iterator current, end;
	boost::tie(current, end) = boost::vertices(graph);
	typedef boost::topo_sort_visitor<std::back_insert_iterator<std::vector<std::size_t> > > TopoVisitor;
	for(; current != end; current++)
	{
		if(colourVector[*current] == Colour::white())
		{
			boost::depth_first_visit(graph, *current, TopoVisitor(std::back_inserter(topologicalSortOutput)), &(colourVector[0]));
		}
	}

	//Now create the new pedigree. 
	for(std::size_t counter = 0; counter < topologicalSortOutput.size(); counter++)
	{
		boost::put(boost::vertex_name, graph, topologicalSortOutput[counter], counter);
	}
	std::vector<std::string> newLineNames(topologicalSortOutput.size());
	std::vector<int> newMother(topologicalSortOutput.size()), newFather(topologicalSortOutput.size());
	for(std::size_t counter = 0; counter < topologicalSortOutput.size(); counter++)
	{
		newLineNames[counter] = lineNames[topologicalSortOutput[counter]];
		graphType::out_edge_iterator current, end;
		boost::tie(current, end) = boost::out_edges(topologicalSortOutput[counter], graph);
		if(current != end)
		{
			if(boost::get(boost::edge_name, graph, *current) == MOTHER)
			{
				newMother[counter] = boost::get(boost::vertex_name, graph, boost::target(*current, graph)) + 1;
				current++;
				newFather[counter] = boost::get(boost::vertex_name, graph, boost::target(*current, graph)) + 1;
			}
			else
			{
				newFather[counter] = boost::get(boost::vertex_name, graph, boost::target(*current, graph)) + 1;
				current++;
				newMother[counter] = boost::get(boost::vertex_name, graph, boost::target(*current, graph)) + 1;
			}
		}
	}

	//Get out the generations for every line. A value of -1 means inconsistent generations, a value of -2 means unknown. 
	Rcpp::IntegerVector generations(topologicalSortOutput.size(), -2);
	for(std::size_t counter = 0; counter < nFounders; counter++)
	{
		colourVector = copiedColourVector;
		boost::depth_first_visit(boost::make_reverse_graph(graph), counter, generationsVisitor(generations), &(colourVector[0]));
	}
	for(std::size_t counter = 0; counter < nFounders; counter++)
	{
		generations[counter] = 0;
	}
	Rcpp::S4 newPedigree("pedigree");
	newPedigree.slot("mother") = Rcpp::wrap(newMother);
	newPedigree.slot("father") = Rcpp::wrap(newFather);
	newPedigree.slot("lineNames") = Rcpp::wrap(newLineNames);
	newPedigree.slot("selfing") = pedigree.slot("selfing");
	newPedigree.slot("warnImproperFunnels") = pedigree.slot("warnImproperFunnels");
	newPedigree.attr("generations") = generations;
	return newPedigree;
END_RCPP
}
#endif
