#ifndef ONTOLOGY_MATCHING_GRAPH
#define ONTOLOGY_MATCHING_GRAPH

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include "util.h"

using namespace std;

class Node {
 public:

  inline Node() { };

  inline Node(string node_name, unsigned int node_id, map<string, unsigned> & geneNamesToIDs) {
    name = node_name;
    id = node_id;
    //descendents.insert(node_id);
    Utils::insertInOrder(descendents, node_id);
    addDescendentToVec(node_id);

    if (geneNamesToIDs.count(node_name) == 1) {
      //genes.insert(geneNamesToIDs[node_name]);
      Utils::insertInOrder(genes, geneNamesToIDs[node_name]);
      isTermGene = true;
    } else {
      isTermGene = false;
    }

    /*if ((node_name.substr(0,3) != "S00") && (node_name.substr(0,5) != "GRMZM") && (node_name.substr(0,2) != "AC") 
	&& (node_name.substr(0,2) != "EF") && (node_name.substr(0,2) != "AF")) {
      isTermGene = false;
    } else {
      isTermGene = true;
      if (geneNamesToIDs.count(node_name) == 1) {
	genes.insert(geneNamesToIDs[node_name]);
      } else {
	unsigned newGeneID = geneNamesToIDs.size();
	geneNamesToIDs[node_name] = newGeneID;
	genes.insert(newGeneID);
      }
      }*/
  }

  inline void addParent(unsigned int newParent) {
    //parents.insert(newParent);
    Utils::insertInOrder(parents, newParent);
  }

  inline vector<unsigned int>::iterator getParentsBegin() { 
    return parents.begin();
  }

  inline vector<unsigned int>::iterator getParentsEnd() { 
    return parents.end();
  }
  
  inline const vector<unsigned int>& getParents() {
    return parents;
  }

  inline bool isParent(unsigned possibleParent) {
    //return parents.count(possibleParent);
    return Utils::elementExists(parents, possibleParent);
  }

  inline bool isChild(unsigned possibleChild) {
    //return children.count(possibleChild);
    return Utils::elementExists(children, possibleChild);
  }
  
  inline bool isDescendent(unsigned possibleDescendent) {
    //return descendents.count(possibleDescendent);
    if (possibleDescendent < descendentsVec.size()) {
      return descendentsVec[possibleDescendent];
    } else {
      return false;
    }
  }

  inline bool isSibling(unsigned possibleSibling) {
    //return siblings.count(possibleSibling);
    if (possibleSibling < siblingsVec.size()) {
      return siblingsVec[possibleSibling];
    } else {
      return false;
    }
  }

  inline unsigned numParents() {
    return parents.size();
  }

  inline void addChild(unsigned int newChild) {
    //children.insert(newChild);
    Utils::insertInOrder(children, newChild);
  }

  inline vector<unsigned int>::iterator getChildrenBegin() { 
    return children.begin();
  }

  inline vector<unsigned int>::iterator getChildrenEnd() { 
    return children.end();
  }

  inline const vector<unsigned int>& getChildren() {
    return children;
  }

  inline unsigned numChildren() {
    return children.size();
  }

  // This only inserts genes into node's gene list.  To maintain gene lists of all genes in graph,
  // use DAGraph::addGenesToAncestors
  inline void addGenes(vector<unsigned>::iterator genesBegin, vector<unsigned>::iterator genesEnd) {
    //genes.insert(genesBegin,genesEnd);
    for (vector<unsigned>::iterator it = genesBegin; it != genesEnd; ++it) {
      Utils::insertInOrder(genes, *it);
    }
  }

  inline vector<unsigned>::iterator getGenesBegin() { 
    return genes.begin();
  }

  inline vector<unsigned>::iterator getGenesEnd() { 
    return genes.end();
  }
  
  inline unsigned numGenes() {
    return genes.size();
  }

  inline void addDescendents(vector<unsigned>::iterator descendentsBegin, vector<unsigned>::iterator descendentsEnd) {
    for (vector<unsigned>::iterator it = descendentsBegin; it != descendentsEnd; ++it) {
      addDescendentToVec(*it);
      Utils::insertInOrder(descendents, *it);
    }
    //descendents.insert(descendentsBegin,descendentsEnd);
  }

  inline void addDescendentToVec(unsigned newID) {
    if (newID >= descendentsVec.size()) {
      descendentsVec.resize(newID+1, false);
    }
    descendentsVec[newID] = true;
  }

  inline void addSibling(unsigned newID) {
    if (newID >= siblingsVec.size()) {
      siblingsVec.resize(newID+1, false);
    }
    siblingsVec[newID] = true;
  }

  inline vector<unsigned>::iterator getDescendentsBegin() { 
    return descendents.begin();
  }

  inline vector<unsigned>::iterator getDescendentsEnd() { 
    return descendents.end();
  }
  
  inline int numDescendents() {
    return descendents.size();
  }

  inline string getName() {
    return name;
  }

  inline unsigned int getID() {
    return id;
  }
  
  inline bool isGene() {
    return isTermGene;
  }

 private:
  string name;
  unsigned int id;
  vector<unsigned int> parents;
  vector<unsigned int> children;
  vector<unsigned int> genes;
  vector<unsigned int> descendents;
  vector<bool> descendentsVec;
  vector<bool> siblingsVec;
  bool isTermGene;
};

class DAGraph {
 public:
  
  inline void addNode(string nodeName, map<string,unsigned> & geneNamesToIDs) {
    unsigned int nodeID = nodes.size();
    nodes.push_back(Node(nodeName, nodeID, geneNamesToIDs));
    nodeNamesToIDs[nodeName] = nodeID;
  }
  
  inline vector<Node>::iterator nodesBegin() { 
    return nodes.begin();
  }
  
  inline vector<Node>::iterator nodesEnd() { 
    return nodes.end();
  }

  inline unsigned numNodes() {
    return nodes.size();
  }

  void addEdge(string parentName, string childName, map<string,unsigned> & geneNamesToIDs, string edgeType = "default") {
    unsigned int parentID;
    unsigned int childID;

    //cout << parentName << "\t" << childName << endl;
    // If parent or child doesn't already exist, add it
    map<string,unsigned int>::iterator parentIt = nodeNamesToIDs.find(parentName);
    if (parentIt == nodeNamesToIDs.end()) {
      addNode(parentName, geneNamesToIDs);
      parentID = getID(parentName);
    } else {
      parentID = parentIt->second;
    }
    map<string,unsigned int>::iterator childIt = nodeNamesToIDs.find(childName);
    if (childIt == nodeNamesToIDs.end()) {
      if ((edgeType == terminalName) && (geneNamesToIDs.count(childName) == 0)) {
	unsigned newGeneID = geneNamesToIDs.size();
	geneNamesToIDs[childName] = newGeneID;
      }
      addNode(childName, geneNamesToIDs);
      childID = getID(childName);
    } else {
      childID = childIt->second;
    }

    // Update sibling lists
    if (keepSibs) {
      for (vector<unsigned>::iterator previousChildrenIt = getChildrenBegin(parentID); 
	   previousChildrenIt != getChildrenEnd(parentID); ++previousChildrenIt) {
	nodes[*previousChildrenIt].addSibling(childID);
	nodes[childID].addSibling(*previousChildrenIt);
      }
    }

    nodes[parentID].addChild(childID);
    nodes[childID].addParent(parentID);
    edgeTypes[make_pair(parentID,childID)] = edgeType;
    
    // Add child's gene list to parent's gene list and to all of parent's ancestors' gene lists
    addGenesToAncestors(nodes[parentID],nodes[childID].getGenesBegin(),nodes[childID].getGenesEnd());
    
    // Add child's descendent list to parent's descendent list and to all of parent's ancestors' descendent list
    addDescendentsToAncestors(nodes[parentID],nodes[childID].getDescendentsBegin(),nodes[childID].getDescendentsEnd());
    return;
  }

  // This version does not include error checking to ensure both nodes already exist.
  // It should only be used to add edges between existing nodes.  This version also does not
  // maintain the descendents or genes lists
  inline void addEdge(unsigned parentID, unsigned childID, string edgeType = "default") {
    nodes[parentID].addChild(childID);
    nodes[childID].addParent(parentID);
    edgeTypes[make_pair(parentID,childID)] = edgeType;
  }

  // Recursively adds gene list to parent and all of its ansectors
  void addGenesToAncestors(Node & parent,
			   vector<unsigned>::iterator genesBegin, 
			   vector<unsigned>::iterator genesEnd) {
    for (vector<unsigned int>::iterator grandparentsIt = parent.getParentsBegin(); 
	 grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addGenesToAncestors(nodes[*grandparentsIt], genesBegin, genesEnd);
    }
    parent.addGenes(genesBegin, genesEnd);
  }

  // Recursively adds descendents list to parent and all of its ansectors
  void addDescendentsToAncestors(Node & parent,
			   vector<unsigned>::iterator descendentsBegin, 
			   vector<unsigned>::iterator descendentsEnd) {
    for (vector<unsigned int>::iterator grandparentsIt = parent.getParentsBegin(); 
	 grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addDescendentsToAncestors(nodes[*grandparentsIt], descendentsBegin, descendentsEnd);
    }
    parent.addDescendents(descendentsBegin, descendentsEnd);
  }

  /*
  // Recursively adds gene list to parent and all of its ansectors.  Handles cycles in graph
  void addGenesToAncestors(vector<unsigned int> & nodesTraversed, Node & parent,
			   set<string>::iterator genesBegin, 
			   set<string>::iterator genesEnd) {
    
    // If this node has already been traversed, then return
    if (nodesTraversed.count(parent.getID()) == 1) {
      return;
    }

    // Otherwise, add the id to the traversed list and continue
    nodesTraversed.insert(parent.getID());
    for (vector<unsigned int>::iterator grandparentsIt = parent.getParentsBegin(); 
	 grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addGenesToAncestors(nodesTraversed, nodes[*grandparentsIt], genesBegin, genesEnd);
    }
    parent.addGenes(genesBegin, genesEnd);
  }
  */

  inline vector<unsigned int>::iterator getParentsBegin(unsigned int id) { 
    return nodes[id].getParentsBegin();
  }

  inline vector<unsigned int>::iterator getParentsEnd(unsigned int id) { 
    return nodes[id].getParentsEnd();
  }
  
  inline vector<unsigned int>::iterator getChildrenBegin(unsigned int id) { 
    return nodes[id].getChildrenBegin();
  }

  inline vector<unsigned int>::iterator getChildrenEnd(unsigned int id) { 
    return nodes[id].getChildrenEnd();
  }

  inline vector<unsigned>::iterator getGenesBegin(unsigned int id) { 
    return nodes[id].getGenesBegin();
  }

  inline vector<unsigned>::iterator getGenesEnd(unsigned int id) { 
    return nodes[id].getGenesEnd();
  }

  inline bool isParent(unsigned possibleParent, unsigned possibleChild) {
    return nodes[possibleChild].isParent(possibleParent);
  }

  inline bool isChild(unsigned possibleChild, unsigned possibleParent) {
    return nodes[possibleParent].isChild(possibleChild);
  }

  inline bool isDescendent(unsigned possibleDescendent, unsigned possibleAncestor) {
    return nodes[possibleAncestor].isDescendent(possibleDescendent);
  }

  inline bool areSiblings(unsigned possibleSibling1, unsigned possibleSibling2) {
    return nodes[possibleSibling1].isSibling(possibleSibling2);
  }

  inline string getName(unsigned int id) {
    return nodes[id].getName();
  }

  inline unsigned int getID(string name) {
    return nodeNamesToIDs[name];
  }

  inline bool isGene(unsigned int id) {
    return nodes[id].isGene();
  }

  inline unsigned numGenesInNode(unsigned int id) {
    return nodes[id].numGenes();
  }

  inline Node & getNode(unsigned int id) {
    return nodes[id];
  }

  inline string getEdgeType(unsigned parent, unsigned child) {
    return edgeTypes[make_pair(parent,child)];
  }

  inline map< pair<unsigned,unsigned>, string >::iterator edgesBegin() {
    return edgeTypes.begin();
  }

  inline map< pair<unsigned,unsigned>, string >::iterator edgesEnd() {
    return edgeTypes.end();
  }

  inline void getGeneIDsToNames(map<string, unsigned> & geneNamesToIDs, vector<string> & geneIDsToNames) {
    geneIDsToNames = vector<string>(geneNamesToIDs.size(), "");
    for (map<string, unsigned>::iterator it = geneNamesToIDs.begin(); it != geneNamesToIDs.end(); ++it) {
      geneIDsToNames[it->second] = it->first;
    }
    return;
  }

  inline void getGeneIDsToNames(vector<string> & geneIDsToNames, map<string, unsigned> & geneNamesToIDs) {
    getGeneIDsToNames(geneNamesToIDs, geneIDsToNames);
    return;
  }

  inline unsigned geneIDToGraphID(unsigned gene_id, vector<string> & geneIDsToNames) {
    return getID(geneIDsToNames[gene_id]);
  }

  inline void setTerminalName(string thisTerminalName) {
    terminalName = thisTerminalName;
    return;
  }

  DAGraph() {
    keepSibs = true;
    terminalName = "gene";
  }
  
  DAGraph(string fileName, map<string, unsigned> & geneNamesToIDs, bool keepSibsHere = true, string thisTerminalName = "gene") {
    keepSibs = keepSibsHere;
    terminalName = thisTerminalName;
    string line;
    ifstream file(fileName.c_str());
    if (file.is_open()) {
      while (file.good()) {
	getline(file,line);
	vector<string> tokens;
	Utils::Tokenize(line, tokens, "\t");
	if (tokens.size() == 2) {
	  addEdge(tokens[0], tokens[1], geneNamesToIDs);
	} else if (tokens.size() >= 3) {
	  addEdge(tokens[0], tokens[1], geneNamesToIDs, tokens[2]);
	}
      }
    }
  }

 private:
  vector<Node> nodes;
  map<string, unsigned int> nodeNamesToIDs;
  map< pair<unsigned,unsigned>, string > edgeTypes;
  bool keepSibs;
  string terminalName;
};

#endif // ONTOLOGY_MATCHING_GRAPH
