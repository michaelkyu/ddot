#ifndef ONTOLOGY_MATCHING_REMOVALLIST
#define ONTOLOGY_MATCHING_REMOVALLIST

#include <map>
#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

typedef pair<unsigned, unsigned> mappingType;

class removalListType {
 public:

  inline void printReasonsForRemoval(mappingType mapping) {
    if (isInRemovalList(mapping)) {
      cout << mapping.first << "," << mapping.second << " has the following removal reasons:" << endl;
      pair< multimap<mappingType,mappingType>::iterator, multimap<mappingType,mappingType>::iterator > removalsToReasonsRange;
      removalsToReasonsRange = removalsToReasons.equal_range(mapping);
      
      for (multimap<mappingType,mappingType>::iterator removalIt = removalsToReasonsRange.first;
	   removalIt != removalsToReasonsRange.second; ++removalIt) {
	mappingType reason = removalIt->second;
	if (!isInRemovalList(reason)) {
	  cout << reason.first << "," << reason.second << endl;
	}
      }
    } else {
      cout << mapping.first << "," << mapping.second << " is not in removal list" << endl;
    }
  }

  inline removalListType(int numNodesInFirstOntology, int numNodesInSecondOntology) {
    mappingInRemovalList.resize(numNodesInFirstOntology, numNodesInSecondOntology, false);
    for (unsigned i = 0; i < mappingInRemovalList.size1(); ++i) {
      for (unsigned j = 0; j < mappingInRemovalList.size2(); ++j) {
	mappingInRemovalList(i,j) = false;
      }
    }
  }

  inline bool addMapping(mappingType mappingToRemove, mappingType reasonForRemoval) {
    removalsToReasons.insert(pair<mappingType,mappingType>(mappingToRemove,reasonForRemoval));
    reasonsToRemovals.insert(pair<mappingType,mappingType>(reasonForRemoval,mappingToRemove));
    if (!isInRemovalList(reasonForRemoval) && !isInRemovalList(mappingToRemove)) {
      mappingInRemovalList(mappingToRemove.first,mappingToRemove.second) = true;
      return true;
    }
    return false;
  }

  inline void removeReason(mappingType reason, vector<mappingType> & mappingsToReturn,
			   vector<mappingType> & mappingsToRemove) {

    vector<mappingType> newMappingsToReturn;
    vector<mappingType> newMappingsToRemove;

    pair< multimap<mappingType,mappingType>::iterator, multimap<mappingType,mappingType>::iterator > reasonsToRemovalsRange;
    reasonsToRemovalsRange = reasonsToRemovals.equal_range(reason);
    // Return all mappings from removalsToReasons which have this as their only live reason listed
    for (multimap<mappingType,mappingType>::iterator reasonIt = reasonsToRemovalsRange.first;
	 reasonIt != reasonsToRemovalsRange.second; ++reasonIt) {
      
      pair< multimap<mappingType,mappingType>::iterator, multimap<mappingType,mappingType>::iterator > removalsToReasonsRange;
      removalsToReasonsRange = removalsToReasons.equal_range(reasonIt->second);

      unsigned numLiveReasons = 0;
      for (multimap<mappingType,mappingType>::iterator removalIt = removalsToReasonsRange.first;
	   removalIt != removalsToReasonsRange.second; ++removalIt) {
	if (!isInRemovalList(removalIt->second)) {
	  ++numLiveReasons;
	}
      }
      // If there is only one live reason for the removal and this is it,
      // add the mapping to the mappingsToReturn
      if (numLiveReasons == 0) {
	newMappingsToReturn.push_back(reasonIt->second);
	//cout << "Switching off " << reasonIt->second.first << "," << reasonIt->second.second << endl;
	mappingInRemovalList(reasonIt->second.first,reasonIt->second.second) = false;
      }
    }
    
    // Add the mappings which conflict with the returning mappings to the mappings to remove
    for (vector<mappingType>::iterator mappingReturnedIt = newMappingsToReturn.begin();
	 mappingReturnedIt != newMappingsToReturn.end(); ++mappingReturnedIt) {

      pair< multimap<mappingType,mappingType>::iterator, multimap<mappingType,mappingType>::iterator > reasonsToRemovalsRange;
      reasonsToRemovalsRange = reasonsToRemovals.equal_range(*mappingReturnedIt);

      for (multimap<mappingType,mappingType>::iterator it = reasonsToRemovalsRange.first;
	   it != reasonsToRemovalsRange.second; ++it) {

	// The mapping being returned is once again a reason for the removal of the mapping it->second
	if (!isInRemovalList(it->second)) {
	  //cout << "Switching on " << it->second.first << "," << it->second.second << endl;
	  mappingInRemovalList(it->second.first,it->second.second) = true;
	  removeReason(it->second, mappingsToReturn, mappingsToRemove);
	  newMappingsToRemove.push_back(it->second);
	}
      }
    }
    
    mappingsToRemove.insert(mappingsToRemove.end(),newMappingsToRemove.begin(),newMappingsToRemove.end());
    mappingsToReturn.insert(mappingsToReturn.end(),newMappingsToReturn.begin(),newMappingsToReturn.end());

  }

  inline bool isInRemovalList(mappingType mapping) {
    return mappingInRemovalList(mapping.first,mapping.second);
  }

 private:
  boost::numeric::ublas::matrix<bool> mappingInRemovalList;

  // removalsToReasons.first = removedMapping
  // removalsToReasons.second = reasonForRemoval
  multimap<mappingType, mappingType> removalsToReasons;

  // reasonsToRemovals.first = reasonForRemoval
  // reasonsToRemovals.second = removedMapping
  multimap<mappingType,mappingType> reasonsToRemovals;

};

#endif
