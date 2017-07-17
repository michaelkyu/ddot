#include "alignmentUtils.h"
#include <algorithm>

void alignmentUtils::printSimMat(boost::numeric::ublas::matrix<float> & m, 
				 DAGraph& g1, DAGraph& g2) 
{
  for (vector<Node>::iterator g2_nodes_it = g2.nodesBegin(); g2_nodes_it != g2.nodesEnd(); ++g2_nodes_it) {
    cout << "\t" << g2_nodes_it->getName();
  }
  cout << endl;
  
  int i = 0;
  for (vector<Node>::iterator g1_nodes_it = g1.nodesBegin(); g1_nodes_it != g1.nodesEnd(); ++g1_nodes_it) {
    int j = 0;
    cout << g1_nodes_it->getName() << "\t";
    for (vector<Node>::iterator g2_nodes_it = g2.nodesBegin(); g2_nodes_it != g2.nodesEnd(); ++g2_nodes_it) {
      cout << m(i,j) << "\t";
      ++j;
    }
    cout << endl;
    ++i;
  }
}

void alignmentUtils::calculateExtSim(boost::numeric::ublas::matrix<float> & m, 
				     DAGraph& g1, DAGraph& g2) {
  int i = 0;
  for (vector<Node>::iterator g1_nodes_it = g1.nodesBegin(); g1_nodes_it != g1.nodesEnd(); ++g1_nodes_it) {
    int j = 0;
    for (vector<Node>::iterator g2_nodes_it = g2.nodesBegin(); g2_nodes_it != g2.nodesEnd(); ++g2_nodes_it) {
      m(i,j) = compGeneSets(*g1_nodes_it,*g2_nodes_it);
      ++j;
    }
    ++i;
    //cout << "i: " << i << endl;
  }
}

// Calculates the relational similarity (similarity of two nodes' parents and children)
void alignmentUtils::addRelSimToSimMat(boost::numeric::ublas::matrix<float> & simMat, 
				       boost::numeric::ublas::matrix<float> & prevSim, 
				       DAGraph & g1, DAGraph & g2, bool genesAsTerms) {
  const float percentRelSim = 0.25;
  const float percentExtSim = (1 - percentRelSim);
  int i = 0;
  for (vector<Node>::iterator g1_nodes_it = g1.nodesBegin(); g1_nodes_it != g1.nodesEnd(); ++g1_nodes_it) {
    int j = 0;
    vector <unsigned> g1_node_children, g1_node_parents;
    if (!genesAsTerms) {

      g1_node_children.reserve(g1_nodes_it->numChildren());
      g1_node_parents.reserve(g1_nodes_it->numParents());
      
      for (vector<unsigned>::iterator it = g1_nodes_it->getChildrenBegin(); it != g1_nodes_it->getChildrenEnd(); ++it) {
	if (!g1.isGene(*it)) {
	  //g1_node_children.insert(*it);
	  Utils::insertInOrder(g1_node_children, *it);
	}
      }
      for (vector<unsigned>::iterator it = g1_nodes_it->getParentsBegin(); it != g1_nodes_it->getParentsEnd(); ++it) {
	if (!g1.isGene(*it)) {
	  //g1_node_parents.insert(*it);
	  Utils::insertInOrder(g1_node_parents, *it);
	}
      }
    }
    
    for (vector<Node>::iterator g2_nodes_it = g2.nodesBegin(); g2_nodes_it != g2.nodesEnd(); ++g2_nodes_it) {
      
      vector<unsigned> g2_node_children, g2_node_parents;
      if (!genesAsTerms) {

	g2_node_children.reserve(g2_nodes_it->numChildren());
	g2_node_parents.reserve(g2_nodes_it->numParents());

	for (vector<unsigned>::iterator it = g2_nodes_it->getChildrenBegin(); it != g2_nodes_it->getChildrenEnd(); ++it) {
	  if (!g2.isGene(*it)) {
	    //g2_node_children.insert(*it);
	    Utils::insertInOrder(g2_node_children, *it);
	  }
	}
	for (vector<unsigned>::iterator it = g2_nodes_it->getParentsBegin(); it != g2_nodes_it->getParentsEnd(); ++it) {
	  if (!g2.isGene(*it)) {
	    //g2_node_parents.insert(*it);
	    Utils::insertInOrder(g2_node_parents, *it);
     	  }
	}
      }
      
      // If both sets don't have parents, just use the child similarity
      bool noParents, noChildren;
      if (genesAsTerms) {
	noParents = ((g1_nodes_it->numParents() == 0) && (g2_nodes_it->numParents() == 0));
	noChildren = ((g1_nodes_it->numChildren() == 0) && (g2_nodes_it->numChildren() == 0));
      } else {
	noParents = ((g1_node_parents.size() == 0) && (g2_node_parents.size() == 0));
	noChildren = ((g1_node_children.size() == 0) && (g2_node_children.size() == 0));
      }
      if ((g1_nodes_it->isGene() && g2_nodes_it->isGene()) || (noParents && noChildren)) {
	// Ignore relSim if both are genes or if there are no parents or children for either node
	
	//simMat(i,j) += 0;	
      } else if (noParents) {
	float childSim;
	if (genesAsTerms) {
	  childSim = calculateSetSim(prevSim, g1_nodes_it->getChildren(), g2_nodes_it->getChildren());
	} else {
	  childSim = calculateSetSim(prevSim, g1_node_children, g2_node_children);
	}
	simMat(i,j) = percentExtSim * simMat(i,j) + childSim * percentRelSim;
      } else if (noChildren) {
	float parentSim;
	if (genesAsTerms) {
	  parentSim = calculateSetSim(prevSim, g1_nodes_it->getParents(), g2_nodes_it->getParents());
	} else {
	  parentSim = calculateSetSim(prevSim, g1_node_parents, g2_node_parents);
	}
	simMat(i,j) = percentExtSim * simMat(i,j) + parentSim * percentRelSim;
      } else {
	float parentSim, childSim;
	if (genesAsTerms) {
	  parentSim = calculateSetSim(prevSim, g1_nodes_it->getParents(), g2_nodes_it->getParents());
	  childSim = calculateSetSim(prevSim, g1_nodes_it->getChildren(), g2_nodes_it->getChildren());
	} else {
	  parentSim = calculateSetSim(prevSim, g1_node_parents, g2_node_parents);
	  childSim = calculateSetSim(prevSim, g1_node_children, g2_node_children);
	}
	simMat(i,j) = percentExtSim * simMat(i,j) + ((parentSim + childSim) / 2) * percentRelSim;
      }
      //}
      ++j;
    }
    //cout << "relsim i: " << i << endl;
    ++i;
  }
}

// Calculates the relational similarity (similarity of two nodes' parents and children)
void alignmentUtils::addRelSimToSimMat(boost::numeric::ublas::matrix<float> & simMat, 
				       boost::numeric::ublas::matrix<float> & prevSim, 
				       DAGraph & g1, DAGraph & g2,
				       boost::numeric::ublas::matrix<bool> & foreverRemoved) {
  const float percentRelSim = 0.25;
  const float percentExtSim = (1 - percentRelSim);
  int i = 0;
  for (vector<Node>::iterator g1_nodes_it = g1.nodesBegin(); g1_nodes_it != g1.nodesEnd(); ++g1_nodes_it) {
    int j = 0;
    for (vector<Node>::iterator g2_nodes_it = g2.nodesBegin(); g2_nodes_it != g2.nodesEnd(); ++g2_nodes_it) {
      //if (foreverRemoved(i,j)) {
      //simMat(i,j) = 0;
      //} else {
      // If both sets don't have parents, just use the child similarity
      bool noParents = ((g1_nodes_it->numParents() == 0) && (g2_nodes_it->numParents() == 0));
      bool noChildren = ((g1_nodes_it->numChildren() == 0) && (g2_nodes_it->numChildren() == 0));
      
      if ((g1_nodes_it->isGene() && g2_nodes_it->isGene()) || (noParents && noChildren)) {
	// Ignore relSim if both are genes or if there are no parents or children for either node
	
	//simMat(i,j) += 0;	
      } else if (noParents) {
	float childSim = calculateSetSim(prevSim, g1_nodes_it->getChildren(), g2_nodes_it->getChildren());
	simMat(i,j) = percentExtSim * simMat(i,j) + childSim * percentRelSim;
      } else if (noChildren) {
	float parentSim = calculateSetSim(prevSim, g1_nodes_it->getParents(), g2_nodes_it->getParents());
	simMat(i,j) = percentExtSim * simMat(i,j) + parentSim * percentRelSim;
      } else {
	float parentSim = calculateSetSim(prevSim, g1_nodes_it->getParents(), g2_nodes_it->getParents());
	float childSim = calculateSetSim(prevSim, g1_nodes_it->getChildren(), g2_nodes_it->getChildren());
	simMat(i,j) = percentExtSim * simMat(i,j) + ((parentSim + childSim) / 2) * percentRelSim;
      }
      //}
      ++j;
    }
    //cout << "relsim i: " << i << endl;
    ++i;
  }
}

bool alignmentUtils::compMappingSimScores(pair< float, mappingType > i,
					  pair< float, mappingType > j) {
  return (i.first > j.first);
}

// This function calculates the similarity between two sets, given a matrix of similarities
// between elements of those sets, and the two sets themselves.  It assumes that the matrix
// of similarities between elements is set up so that simMat(i,j) gives the similarity
// between element i of set1 and element j of set2.  The same indices used for the matrix
// are used to specify the elements of the sets.
float alignmentUtils::calculateSetSim(boost::numeric::ublas::matrix<float> & simMat,
				      const vector<unsigned> & set1, const vector<unsigned> & set2) {

  // A vector which holds the similarity scores associated with each potential mapping between
  // elements of the two sets.  mappingSimScores.second.first gives the index of the element
  // of the first set.  mappingSimScores.second.second gives the index of the element of the
  // second set.  mappingSimScores.first gives the similarity score of the mapping
  vector< pair< float, mappingType > > mappingSimScores;
  mappingSimScores.reserve(set1.size()*set2.size());

  for (vector<unsigned int>::const_iterator set1_it = set1.begin();
       set1_it != set1.end(); ++set1_it) {
    for (vector<unsigned int>::const_iterator set2_it = set2.begin();
	 set2_it != set2.end(); ++set2_it) {
      mappingType mapping(*set1_it,*set2_it);
      float mappingScore = simMat(*set1_it,*set2_it);
      if (mappingScore > 0) {
	pair<float, mappingType > mappingAndScore(mappingScore, mapping);
	mappingSimScores.push_back(mappingAndScore);
      }
    }
  }
  
  // Order mappings from most similar to least similar
  sort(mappingSimScores.begin(), mappingSimScores.end(), compMappingSimScores);

  // Greedily add most similar mapping until all elements of (at least) one
  // of the sets have been mapped.  Allow only one mapping per element

  float sumOfSims = 0;

  
  // use vector<unsigned>
  vector<unsigned> set1_entities_mapped;
  vector<unsigned> set2_entities_mapped;
  for (vector< pair< float, mappingType > >::iterator vecIt = mappingSimScores.begin();
       vecIt != mappingSimScores.end(); ++vecIt) {
    // Check to make sure neither element already has a mapping

    if ((!Utils::elementExists(set1_entities_mapped,vecIt->second.first)) && 
	(!Utils::elementExists(set2_entities_mapped,vecIt->second.second))) {
      //if ((set1_entities_mapped.count(vecIt->second.first) == 0) && 
      //(set2_entities_mapped.count(vecIt->second.second) == 0)) {
      sumOfSims += vecIt->first;
      //set1_entities_mapped.insert(vecIt->second.first);
      //set2_entities_mapped.insert(vecIt->second.second);
      Utils::insertInOrder(set1_entities_mapped, vecIt->second.first);
      Utils::insertInOrder(set2_entities_mapped, vecIt->second.second);
      
      // If all of one of the elements of a set have been mapped, quit looking
      if ((set1_entities_mapped.size() == set1.size()) ||
	  (set2_entities_mapped.size() == set2.size())) {
	break;
      }
    }  
  }
  

  /*
  // Use vector<bool>
  unsigned num_entities_mapped = 0;
  //vector<bool> set1_entities_mapped(simMat.size1(),false);
  //vector<bool> set2_entities_mapped(simMat.size2(),false);
  vector<bool> set1_entities_mapped(set1max+1,false);
  vector<bool> set2_entities_mapped(set2max+1,false);
  
  for (vector< pair< float, mappingType > >::iterator vecIt = mappingSimScores.begin();
       vecIt != mappingSimScores.end(); ++vecIt) {
    // Check to make sure neither element already has a mapping
    if (!set1_entities_mapped[vecIt->second.first] && 
	!set2_entities_mapped[vecIt->second.second]) {
      sumOfSims += vecIt->first;
      set1_entities_mapped[vecIt->second.first] = true;
      set2_entities_mapped[vecIt->second.second] = true;
      ++num_entities_mapped;
      
      // If all of one of the elements of a set have been mapped, quit looking
      if ((num_entities_mapped == set1.size()) ||
	  (num_entities_mapped == set2.size())) {
	break;
      }
    }  
  }
  */

  // Return the Jaccard (size of intersection over size of union of the two sets)
  float result = (sumOfSims / (set1.size() + set2.size() - sumOfSims));
  return result;
}

void alignmentUtils::getPreAlignment(boost::numeric::ublas::matrix<float> & simMat,
				     alignmentType & preAlignment, float threshold,
				     removalListType & removalList,
				     set<unsigned> & nodesIn1ToRealign,
				     set<unsigned> & nodesIn2ToRealign,
				     set<mappingType> & mappingsAdded,
				     DAGraph & g1, DAGraph& g2,
				     float minSimilarityForPrealignment) {
  if (nodesIn1ToRealign.size() == 0) {
    for (unsigned i = 0; i < simMat.size1(); ++i) {
      nodesIn1ToRealign.insert(i);
    }
  }

  if (nodesIn2ToRealign.size() == 0) {
    for (unsigned i = 0; i < simMat.size2(); ++i) {
      nodesIn2ToRealign.insert(i);
    }
  }

  // Select all mappings for each element i of ontology 1 that are within
  // threshold of maximum similarity mapping
  for (set<unsigned>::iterator i = nodesIn1ToRealign.begin(); 
       i != nodesIn1ToRealign.end(); ++i) {
    set<mappingType> mappingsToAdd;
    float maxSim = -10;
    for (unsigned j = 0; j < simMat.size2(); ++j) {
      if (g1.isGene(*i) != g2.isGene(j)) {
	continue;
      }
      if ((simMat(*i,j) > minSimilarityForPrealignment) && !removalList.isInRemovalList(mappingType(*i,j))) {
	if (simMat(*i,j) > (maxSim + threshold)) {
	  mappingsToAdd.clear();
	  mappingsToAdd.insert(mappingType(*i,j));
	  maxSim = simMat(*i,j);
	  
	} else if (simMat(*i,j) > maxSim) {
	  
	  // Keep only mappings with a similarity within the threshold
	  // of new maximum similarity
	  set<mappingType> newMappingsToAdd;
	  for (set<mappingType>::iterator it = mappingsToAdd.begin();
	       it != mappingsToAdd.end(); ++it) {
	    if (simMat(*i,j) <= (simMat(it->first,it->second) + threshold)) {
	      newMappingsToAdd.insert(*it);
	    }
	  }
	  mappingsToAdd = newMappingsToAdd;
	  maxSim = simMat(*i,j);

	} else if (simMat(*i,j) >= (maxSim - threshold)) {
	  mappingsToAdd.insert(mappingType(*i,j));
	}
      }
    }
    //cout << "preAlignment i: " << i << endl;
    //double thisMaxSim = -1000;
    //double minSim = 1000;
    //long numMappings = 0;
    for (set<mappingType>::iterator mapIt = mappingsToAdd.begin(); 
	 mapIt != mappingsToAdd.end(); ++mapIt) {
      
      /*double thisSim = simMat(mapIt->first,mapIt->second);
      if (thisSim > thisMaxSim) {
	thisMaxSim = thisSim;
      }
      if (thisSim < minSim) {
	minSim = thisSim;
      }
      ++numMappings;
      */
      if (preAlignment.insert(*mapIt).second) {
	mappingsAdded.insert(*mapIt);
      }
    }
    //cout << "preAlignment i: " << *i << "\tnumMappings: " << numMappings << "\tmaxSim: " << thisMaxSim << "\tminSim: " << minSim << endl;
      
    //preAlignment.insert(mappingsToAdd.begin(), mappingsToAdd.end());
    //mappingsAdded.insert(mappingsToAdd.begin(), mappingsToAdd.end());
  }

  // Select all mappings for each element j of ontology 2 that are within
  // threshold of maximum similarity mapping
  for (set<unsigned>::iterator j = nodesIn2ToRealign.begin(); 
       j != nodesIn2ToRealign.end(); ++j) {
    set<mappingType> mappingsToAdd;
    float maxSim = -10;
    for (unsigned i = 0; i < simMat.size1(); ++i) {
      if (g1.isGene(i) != g2.isGene(*j)) {
	continue;
      }
      if ((simMat(i,*j) > minSimilarityForPrealignment) && !removalList.isInRemovalList(mappingType(i,*j))) {
	if (simMat(i,*j) > (maxSim + threshold)) {
	  mappingsToAdd.clear();
	  mappingsToAdd.insert(mappingType(i,*j));
	  maxSim = simMat(i,*j);
	} else if (simMat(i,*j) > maxSim) {
	  
	  // Keep only mappings with a similarity within the threshold
	  // of new maximum similarity
	  set<mappingType> newMappingsToAdd;
	  for (set<mappingType>::iterator it = mappingsToAdd.begin();
	       it != mappingsToAdd.end(); ++it) {
	    if (simMat(i,*j) <= (simMat(it->first,it->second) + threshold)) {
	      newMappingsToAdd.insert(*it);
	    }
	  }
	  mappingsToAdd = newMappingsToAdd;
	  maxSim = simMat(i,*j);
	  
	} else if (simMat(i,*j) >= (maxSim - threshold)) {
	  mappingsToAdd.insert(mappingType(i,*j));
	}
      }
    }

    /*
    double thisMaxSim = -1000;
    double minSim = 1000;
    long numMappings = 0;
    */
    for (set<mappingType>::iterator mapIt = mappingsToAdd.begin(); 
	 mapIt != mappingsToAdd.end(); ++mapIt) {
      /*
      double thisSim = simMat(mapIt->first,mapIt->second);
      if (thisSim > thisMaxSim) {
	thisMaxSim = thisSim;
      }
      if (thisSim < minSim) {
	minSim = thisSim;
      }
      ++numMappings;
      */
      if (preAlignment.insert(*mapIt).second) {
	mappingsAdded.insert(*mapIt);
      }
    }
    //cout << "preAlignment j: " << *j << "\tnumMappings: " << numMappings << "\tmaxSim: " << thisMaxSim << "\tminSim: " << minSim << endl;

    //cout << "preAlignment j: " << j << endl;
    //preAlignment.insert(mappingsToAdd.begin(), mappingsToAdd.end());
    //mappingsAdded.insert(mappingsToAdd.begin(), mappingsToAdd.end());
  }
}

bool alignmentUtils::semanticVerification(alignmentType & preAlignment,
					  alignmentType & verifiedAlignment,
					  DAGraph & g1, DAGraph & g2,
					  boost::numeric::ublas::matrix<float> & simMat,
					  removalListType & removalList,
					  set<unsigned> & nodesChangedIn1,
					  set<unsigned> & nodesChangedIn2,
					  set<mappingType> & mappingsAdded,
					  int verificationMode) {
  nodesChangedIn1.clear();
  nodesChangedIn2.clear();
  cout << "mappings in prealignment: " << preAlignment.size() << endl;
  cout << "new mappings added: " << mappingsAdded.size() << endl;
  bool alignmentVerified = true;
  alignmentType lastAlignment = verifiedAlignment;
  verifiedAlignment = preAlignment;
  
  
  for (set<mappingType>::iterator firstMap = preAlignment.begin();
       firstMap != preAlignment.end(); ++firstMap) {
    set<mappingType>::iterator secondMap = firstMap;
    ++secondMap;
    for (/*secondMap already initialized*/; secondMap != preAlignment.end(); ++secondMap) {
      verifyPair(*firstMap,*secondMap,verifiedAlignment,
		 g1,g2,simMat,removalList,
		 nodesChangedIn1,nodesChangedIn2,
		 alignmentVerified,verificationMode);
    }
  }
  

  /*
  set<mappingType> mappingsAddedDuringVerification;
  for (set<mappingType>::iterator firstMap = mappingsAdded.begin();
       firstMap != mappingsAdded.end(); ++firstMap) {
    set<mappingType>::iterator secondMap = firstMap;
    ++secondMap;
    for (/*secondMap already initialized/; secondMap != mappingsAdded.end(); ++secondMap) {
      verifyPair(*firstMap,*secondMap,verifiedAlignment,
		 g1,g2,simMat,removalList,
		 nodesChangedIn1,nodesChangedIn2,mappingsAddedDuringVerification,
		 alignmentVerified);
    }
  }

  
  for (set<mappingType>::iterator firstMap = mappingsAdded.begin(); 
       firstMap != mappingsAdded.end(); ++firstMap) {
    for (alignmentType::iterator secondMap = lastAlignment.begin(); 
	 secondMap != lastAlignment.end(); ++secondMap) {
      verifyPair(*firstMap,*secondMap,verifiedAlignment,
		 g1,g2,simMat,removalList,
      		 nodesChangedIn1,nodesChangedIn2,mappingsAddedDuringVerification,
      		 alignmentVerified);
      //++mappingsVerified;
      //if ((mappingsVerified % 100000) == 0) {
      //	cout << "verified " << mappingsVerified << endl;
      //}
    }
  }
  
  set<mappingType> mappingsAddedDuringReverification = mappingsAddedDuringVerification;
  while (mappingsAddedDuringReverification.size() != 0) {
    lastAlignment = verifiedAlignment;
    mappingsAddedDuringVerification = mappingsAddedDuringReverification;
    mappingsAddedDuringReverification.clear();
    for (set<mappingType>::iterator firstMap = mappingsAddedDuringVerification.begin();
	 firstMap != mappingsAddedDuringVerification.end(); ++firstMap) {
      set<mappingType>::iterator secondMap = firstMap;
      ++secondMap;
      for (/*secondMap already initialized/; secondMap != mappingsAddedDuringVerification.end(); ++secondMap) {
	verifyPair(*firstMap,*secondMap,verifiedAlignment,
		   g1,g2,simMat,removalList,
		   nodesChangedIn1,nodesChangedIn2,mappingsAddedDuringReverification,
		   alignmentVerified);
      }
    }
  
    for (set<mappingType>::iterator firstMap = mappingsAddedDuringVerification.begin(); 
	 firstMap != mappingsAddedDuringVerification.end(); ++firstMap) {
      for (alignmentType::iterator secondMap = lastAlignment.begin(); 
	   secondMap != lastAlignment.end(); ++secondMap) {
	verifyPair(*firstMap,*secondMap,verifiedAlignment,
		   g1,g2,simMat,removalList,
		   nodesChangedIn1,nodesChangedIn2,mappingsAddedDuringReverification,
		   alignmentVerified);
      }
    }
  }
  */

  return alignmentVerified;
}

bool sortMappingsDescending(const mappingAndScore i, const mappingAndScore j) {
  return (i.second > j.second);
}

bool alignmentUtils::getAlignment(boost::numeric::ublas::matrix<float> & simMat,
				  alignmentType & alignment, float threshold,
				  DAGraph& g1,DAGraph& g2,
				  float minSimilarityForPrealignment, 
				  int verificationMode,
				  bool allowMultiples) {
  vector<mappingAndScore> sortedMappings;
  sortedMappings.reserve(simMat.size1()*simMat.size2());
  for (int i = 0; i < simMat.size1(); ++i) {
    for (int j = 0; j < simMat.size2(); ++j) {
      if ((simMat(i,j) > minSimilarityForPrealignment) && (g1.isGene(i) == g2.isGene(j))) {
	sortedMappings.push_back(make_pair(make_pair(i,j),simMat(i,j)));
      }
    }
  }
  sort(sortedMappings.begin(), sortedMappings.end(), sortMappingsDescending);

  vector<bool> nodesMappedinG1(g1.numNodes(),false);
  vector<bool> nodesMappedinG2(g2.numNodes(),false);

  float lastScore = 2;
  vector<pair<int,int> > mappingsToAdd;

  for (vector<mappingAndScore>::iterator possibleMappingsIt = sortedMappings.begin(); possibleMappingsIt != sortedMappings.end(); ++possibleMappingsIt) {
    
    
    if (allowMultiples && (possibleMappingsIt->second < lastScore)) {
      for (vector<pair<int,int> >::iterator mappingsToAddIt = mappingsToAdd.begin(); 
	   mappingsToAddIt != mappingsToAdd.end(); ++mappingsToAddIt) {
	nodesMappedinG1[mappingsToAddIt->first] = true;
	nodesMappedinG2[mappingsToAddIt->second] = true;
	alignment.insert(*mappingsToAddIt);
      }
      lastScore = possibleMappingsIt->second;
      mappingsToAdd.clear();
    }
    
    
    if (!nodesMappedinG1[possibleMappingsIt->first.first] && !nodesMappedinG2[possibleMappingsIt->first.second]) {
      
      bool conflictFound = false;
      for (set<mappingType>::iterator alignIt = alignment.begin(); alignIt != alignment.end(); ++alignIt) {
	if ((possibleMappingsIt->second < simMat(alignIt->first,alignIt->second)) && 
	    verifyPair(possibleMappingsIt->first,*alignIt,g1,g2,verificationMode)) {
	  conflictFound = true;
	  break;
	}
      }
      
      if (!conflictFound) {
	if (allowMultiples) {
	  mappingsToAdd.push_back(possibleMappingsIt->first);
	} else {
	  nodesMappedinG1[possibleMappingsIt->first.first] = true;
	  nodesMappedinG2[possibleMappingsIt->first.second] = true;
	  alignment.insert(possibleMappingsIt->first);
	}
      }
    }
  }
  
  if (allowMultiples) {
    for (vector<pair<int,int> >::iterator mappingsToAddIt = mappingsToAdd.begin(); mappingsToAddIt != mappingsToAdd.end(); ++mappingsToAddIt) {
      nodesMappedinG1[mappingsToAddIt->first] = true;
      nodesMappedinG2[mappingsToAddIt->second] = true;
      alignment.insert(*mappingsToAddIt);
    }
  }
  
  cout << "Mappings in alignment: " << alignment.size() << endl;
  return true;
}
