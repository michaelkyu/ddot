#include "collapseRedundantNodes.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "Needs 1 argument - ontology to collapse" << endl;
    cout << "Optional second argument - minimum node size to keep (default is 1)" << endl;
    return 0;
  }
  map<string, unsigned> geneNamesToIDs;
  string ontology1 = argv[1];
  DAGraph original_graph(ontology1, geneNamesToIDs, false);
  DAGraph collapsed_graph(ontology1, geneNamesToIDs, false);
  unsigned minToKeep = 1;
  if (argc >= 3) {
    minToKeep = atoi(argv[2]);
  }

  vector<string> geneIDsToNames;
  original_graph.getGeneIDsToNames(geneNamesToIDs, geneIDsToNames);

  vector<bool> nodesToEliminate(original_graph.numNodes(),false);

  for (vector<Node>::iterator nodesIt = original_graph.nodesBegin(); nodesIt != original_graph.nodesEnd(); ++nodesIt) {
    if (nodesIt->numGenes() == 0) {
      nodesToEliminate[nodesIt->getID()] = true;
    } else if ((!nodesIt->isGene()) && (nodesIt->numGenes() < minToKeep)) {
      for (vector<unsigned>::iterator parentsIt = nodesIt->getParentsBegin(); parentsIt != nodesIt->getParentsEnd(); ++parentsIt) {
	//cout << "Adding " << original_graph.getName(*parentsIt) << " to " << original_graph.getName(original_graph.geneIDToGraphID(*(nodesIt->getGenesBegin()), geneIDsToNames)) << endl;
	for (vector<unsigned>::iterator geneIt = nodesIt->getGenesBegin(); geneIt != nodesIt->getGenesEnd(); ++geneIt) {
	  collapsed_graph.addEdge(*parentsIt, original_graph.geneIDToGraphID(*geneIt, geneIDsToNames), "gene");
	}
      }
      nodesToEliminate[nodesIt->getID()] = true;
    } else {
      set< pair<unsigned,string> > nodesToCollapseHere;
      getParentsToCollapseHere(nodesIt->getID(),original_graph,nodesToCollapseHere);

      for (set< pair<unsigned, string> >::iterator collapseIt = nodesToCollapseHere.begin();
	   collapseIt != nodesToCollapseHere.end(); ++collapseIt) {
	for (vector<unsigned>::iterator nodesAboveCollapsingNodesIt = collapsed_graph.getParentsBegin(collapseIt->first);
	     nodesAboveCollapsingNodesIt != collapsed_graph.getParentsEnd(collapseIt->first);
	     ++nodesAboveCollapsingNodesIt) {
	  
	  for (vector<unsigned>::iterator nodesBelowCollapsingNodesIt = collapsed_graph.getChildrenBegin(collapseIt->first);
	       nodesBelowCollapsingNodesIt != collapsed_graph.getChildrenEnd(collapseIt->first);
	       ++nodesBelowCollapsingNodesIt) {
	    
	    //string oldEdgeType1 = collapseIt->second;
	    string oldEdgeType1 = collapsed_graph.getEdgeType(collapseIt->first, *nodesBelowCollapsingNodesIt);
	    string oldEdgeType2 = collapsed_graph.getEdgeType(*nodesAboveCollapsingNodesIt, collapseIt->first);
	    string newEdgeType;
	    if (original_graph.isGene(*nodesBelowCollapsingNodesIt)) {
	      newEdgeType = "gene";
	    } else if ((oldEdgeType1 == "regulates") || (oldEdgeType2 == "regulates")) {
	      newEdgeType = "regulates";
	    } else if ((oldEdgeType1 == "part_of") || (oldEdgeType2 == "part_of")) {
	      newEdgeType = "part_of";
	    } else if ((oldEdgeType1 == "has_part") || (oldEdgeType2 == "has_part")) {
	      newEdgeType = "has_part";
	    } else {
	      newEdgeType = "is_a";
	    }
	    //collapsed_graph.addEdge(*nodesAboveCollapsingNodesIt, nodesIt->getID(), newEdgeType);
	    collapsed_graph.addEdge(*nodesAboveCollapsingNodesIt, *nodesBelowCollapsingNodesIt, newEdgeType);
	  }
	}
	
	nodesToEliminate[collapseIt->first] = true;
      }
    }
    
  }

  for (unsigned i = 0; i < collapsed_graph.numNodes(); ++i) {
    if (!nodesToEliminate[i]) {
      for (vector<unsigned>::iterator childIt = collapsed_graph.getChildrenBegin(i);
	   childIt != collapsed_graph.getChildrenEnd(i); ++childIt) {
	if (!nodesToEliminate[*childIt]) {
	  cout << collapsed_graph.getName(i) << "\t" << collapsed_graph.getName(*childIt) << "\t" << collapsed_graph.getEdgeType(i,*childIt) << endl;
	}
      }
    }
  }
  
  return 1;
}
