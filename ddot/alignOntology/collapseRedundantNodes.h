#include <stdlib.h>
#include <iostream>
#include "graph.h"

inline bool areGenesInNodesSame(unsigned node1, unsigned node2, DAGraph & graph) {
  vector<unsigned>::iterator node1GenesIt = graph.getGenesBegin(node1);
  vector<unsigned>::iterator node2GenesIt = graph.getGenesBegin(node2);
  bool areSame = (!graph.isGene(node1) && !graph.isGene(node2));
  while (areSame && (node1GenesIt != graph.getGenesEnd(node1)) && (node2GenesIt != graph.getGenesEnd(node2))) {
    if (*node1GenesIt != *node2GenesIt) {
      //cout << "Different" << endl;
      areSame = false;
    }
    ++node1GenesIt;
    ++node2GenesIt;
  }
  if (areSame && ((node1GenesIt != graph.getGenesEnd(node1)) || (node2GenesIt != graph.getGenesEnd(node2)))) {
    //cout << "Different" << endl;
    areSame = false;
  }
  return areSame;
}

void getParentsToCollapseHere(unsigned nodeID, DAGraph & graph, set<pair<unsigned, string> > & parentsToCollapseHere, string edgeTypesToHere = "") {
  for (vector<unsigned>::iterator parentIt = graph.getParentsBegin(nodeID);
       parentIt != graph.getParentsEnd(nodeID); ++parentIt) {
    bool areSame = areGenesInNodesSame(nodeID, *parentIt, graph);
    if (areSame) {
      if (edgeTypesToHere != "") {
	edgeTypesToHere += ",";
      }
      edgeTypesToHere += graph.getEdgeType(*parentIt,nodeID);
      parentsToCollapseHere.insert(make_pair(*parentIt,edgeTypesToHere));
      //cout << graph.getName(*parentIt) << "\t" << graph.getName(nodeID) << endl;
      //getParentsToCollapseHere(*parentIt, graph, parentsToCollapseHere,edgeTypesToHere);
    }
  }
}

void getChildrenToCollapseHere(unsigned nodeID, DAGraph & graph, set<unsigned> & childrenToCollapseHere) {
  for (vector<unsigned>::iterator childIt = graph.getChildrenBegin(nodeID);
       childIt != graph.getChildrenEnd(nodeID); ++childIt) {
    //cout << "parent: " << graph.getName(nodeID) << "\tchild: " << graph.getName(*childIt) << endl;
    bool areSame = areGenesInNodesSame(nodeID, *childIt, graph);
    if (areSame) {
      childrenToCollapseHere.insert(*childIt);
      getChildrenToCollapseHere(*childIt, graph, childrenToCollapseHere);
    }
  }
}

void getSiblingsToCollapseHere(unsigned nodeID, DAGraph & graph, set<unsigned> & siblingsToCollapseHere) {
  for (vector<unsigned>::iterator parentIt = graph.getParentsBegin(nodeID);
       parentIt != graph.getParentsEnd(nodeID); ++parentIt) {
    for (vector<unsigned>::iterator siblingIt = graph.getChildrenBegin(*parentIt);
	 siblingIt != graph.getChildrenEnd(*parentIt); ++siblingIt) {
      if (*siblingIt != nodeID) {
	bool areSame = areGenesInNodesSame(nodeID, *siblingIt, graph);
	if (areSame) {
	  siblingsToCollapseHere.insert(*siblingIt);
	  getChildrenToCollapseHere(*siblingIt, graph, siblingsToCollapseHere);
	}
      }
    }
  }
}
