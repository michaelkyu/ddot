#include <iostream>
#include "graph.h"

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "Needs 1 argument - ontology" << endl;
    cout << "Optional second argument: 'size', 'descendents', 'ancestors', 'genes', 'both' or 'tree'" << endl;
    cout << "Optional third argument: identifier of leaf node - default is 'gene'" << endl;
    return 0;
  }
  map<string, unsigned> geneNamesToIDs;
  string ontology1 = argv[1];
  string opt = "size";
  string terminalName = "gene";
  if (argc >= 3) {
    opt = argv[2];
    if (argc >= 4) {
      terminalName = argv[3];
    }
  }

  DAGraph g1(ontology1, geneNamesToIDs, false, terminalName);

  for (unsigned i = 0; i < g1.numNodes(); ++i) {
    if (opt == "tree") {
      vector<unsigned>::iterator parentIt = g1.getParentsBegin(i);
      if (parentIt != g1.getParentsEnd(i)) {
	unsigned smallestParent = *parentIt;
	unsigned sizeSmallestParent = g1.numGenesInNode(*parentIt);
	++parentIt;
	for ( ; parentIt != g1.getParentsEnd(i); ++parentIt) {
	  if (g1.numGenesInNode(*parentIt) < sizeSmallestParent) {
	    smallestParent = *parentIt;
	    sizeSmallestParent = g1.numGenesInNode(*parentIt);
	  }
	}
	cout << g1.getName(smallestParent) << "\t" << g1.getName(i) << endl;
      }
    } else if (!g1.isGene(i)) {
      cout << g1.getName(i) << "\t";
      if (opt == "size") {
	cout << g1.numGenesInNode(i);
      } else if (opt == "descendents") {
	for (unsigned j = 0; j < g1.numNodes(); ++j) {
	  if (!g1.isGene(j) && g1.isDescendent(j,i)) {
	    cout << g1.getName(j) << ",";
	  }
	}
      } else if (opt == "genes") {
	cout << g1.numGenesInNode(i) << "\t";
	for (unsigned j = 0; j < g1.numNodes(); ++j) {
	  if (g1.isGene(j) && ((g1.isDescendent(i,j) || g1.isDescendent(j,i)))) {
	    cout << g1.getName(j) << ",";
	  }
	}
      } else if (opt == "ancestors") {
	for (unsigned j = 0; j < g1.numNodes(); ++j) {
	  if (!g1.isGene(j) && g1.isDescendent(i,j)) {
	    cout << g1.getName(j) << ",";
	  }
	}
      } else if (opt == "both") {
	for (unsigned j = 0; j < g1.numNodes(); ++j) {
	  if (!g1.isGene(j) && ((g1.isDescendent(i,j) || g1.isDescendent(j,i)))) {
	    cout << g1.getName(j) << ",";
	  }
	}
      } else if (opt == "genes_and_descendents") {
	for (unsigned j = 0; j < g1.numNodes(); ++j) {
	  //if (g1.isGene(j) && ((g1.isDescendent(i,j) || g1.isDescendent(j,i)))) {
	  //if (!g1.isGene(j) && ((g1.isDescendent(i,j) || g1.isDescendent(j,i)))) {
	  if (g1.isDescendent(j,i)) {
	    cout << g1.getName(j) << ",";
	  }
	}
      } else {
	cout << "Option input was " << opt << endl;
	cout << "Possible options are 'size', 'descendents', 'ancestors', 'genes' or 'both'";
      }
      cout << endl;
    }
  }
  
  return 1;
}
