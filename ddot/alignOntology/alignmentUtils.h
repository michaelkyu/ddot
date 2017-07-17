#ifndef ONTOLOGY_MATCHING_ALIGNMENTUTILS
#define ONTOLOGY_MATCHING_ALIGNMENTUTILS

#include <boost/numeric/ublas/matrix.hpp>
#include "graph.h"
#include "removalList.h"

typedef set<mappingType> alignmentType;
typedef pair<mappingType, float> mappingAndScore;

class alignmentUtils {
 public:
  static inline int numGenesInIntersection(Node & n1, Node & n2) {
    if ((n1.numGenes() == 0) && (n2.numGenes() == 0)) {
      return 0;
    }
    
    vector<unsigned>::iterator first1 = n1.getGenesBegin();
    vector<unsigned>::iterator last1 = n1.getGenesEnd();
    vector<unsigned>::iterator first2 = n2.getGenesBegin();
    vector<unsigned>::iterator last2 = n2.getGenesEnd();
    
    int numInIntersection = 0;
    while (first1!=last1 && first2!=last2)
      {
	if (*first1<*first2) ++first1;
	else if (*first2<*first1) ++first2;
	else { ++numInIntersection; ++first1; ++first2;}
      }
    return numInIntersection;
  }

  static inline float compGeneSets(Node & n1, Node & n2) {
    int numInIntersection = numGenesInIntersection(n1,n2);
    return (static_cast<float>(numInIntersection) / (n1.numGenes() + n2.numGenes() - numInIntersection));
  }

  static void printSimMat(boost::numeric::ublas::matrix<float> & m, DAGraph& g1, DAGraph& g2);
  
  static inline void printMat(boost::numeric::ublas::matrix<float> & m) {
    for (unsigned k = 0; k < m.size1(); ++k) {
      for (unsigned l = 0; l < m.size2(); ++l) {
	cout << m(k,l) << "\t";
      }
      cout << endl;
    }
  }

  /*
  static inline bool compAligns(alignmentType & align1, alignmentType & align2) {
    if (align1.size() != align2.size()) {
      return false;
    }
    for (alignmentType::iterator it = align1.begin(); it != align1.end(); ++it) {
      if (align2.count(*it) == 0) {
	return false;
      }
    }
    return true;
  }
  */

  static void calculateExtSim(boost::numeric::ublas::matrix<float> & m, DAGraph& g1, DAGraph& g2);

  static void addRelSimToSimMat(boost::numeric::ublas::matrix<float> & simMat, 
				boost::numeric::ublas::matrix<float> & prevSim, 
				DAGraph & g1, DAGraph & g2,
				boost::numeric::ublas::matrix<bool> & foreverRemoved);

  static void addRelSimToSimMat(boost::numeric::ublas::matrix<float> & simMat, 
				boost::numeric::ublas::matrix<float> & prevSim, 
				DAGraph & g1, DAGraph & g2, bool genesAsTerms = true);
  
  static float calculateSetSim(boost::numeric::ublas::matrix<float> & prevSim,
			       const vector<unsigned> & set1, const vector<unsigned> & set2);

  static void getPreAlignment(boost::numeric::ublas::matrix<float> & simMat,
			      alignmentType & preAlignment, float threshold,
			      removalListType & removalList,
			      set<unsigned> & nodesIn1ToRealign,
			      set<unsigned> & nodesIn2ToRealign,
			      set<mappingType> & mappingsAdded,
			      DAGraph& g1, DAGraph& g2, 
			      float minSimilarityForPrealignment);

  static bool semanticVerification(alignmentType & preAlignment,
				   alignmentType & verifiedAlignment,
				   DAGraph & g1, DAGraph & g2,
				   boost::numeric::ublas::matrix<float> & simMat,
				   removalListType & removalList,
				   set<unsigned> & nodesChangedIn1,
				   set<unsigned> & nodesChangedIn2,
				   set<mappingType> & mappingsAdded,
				   int verificationMode);

  static bool getAlignment(boost::numeric::ublas::matrix<float> & simMat,
			   alignmentType & alignment, float threshold,
			   DAGraph& g1,DAGraph& g2,
			   float minSimilarityForPrealignment, 
			   int verificationMode, bool allowMultiples = false);

  static bool compMappingSimScoresDouble(pair< double, mappingType > i,
					 pair< double, mappingType > j) {
    return i.first < j.first;
  }
  
 private:
  static bool compMappingSimScores(pair< float, mappingType > i,
				   pair< float, mappingType > j);

  static inline void removeConflict(mappingType mappingToRemove, mappingType reason,
				    alignmentType & verifiedAlignment,
				    removalListType & removalList,
				    set<unsigned> & nodesChangedIn1,
				    set<unsigned> & nodesChangedIn2) {
    if (removalList.addMapping(mappingToRemove,reason)) {
      verifiedAlignment.erase(mappingToRemove);
      vector<mappingType> mappingsToReturn;
      vector<mappingType> mappingsToRemove;
      removalList.removeReason(mappingToRemove, mappingsToReturn, mappingsToRemove);
      
      nodesChangedIn1.insert(mappingToRemove.first);
      nodesChangedIn2.insert(mappingToRemove.second);

      for (vector<mappingType>::iterator mappingsReturnedIt = mappingsToReturn.begin();
	   mappingsReturnedIt != mappingsToReturn.end(); ++mappingsReturnedIt) {
	nodesChangedIn1.insert(mappingsReturnedIt->first);
	nodesChangedIn2.insert(mappingsReturnedIt->second);
      }

      verifiedAlignment.insert(mappingsToReturn.begin(), mappingsToReturn.end());

      for (vector<mappingType>::iterator mappingsToRemoveIt = mappingsToRemove.begin();
	   mappingsToRemoveIt != mappingsToRemove.end(); ++mappingsToRemoveIt) {
	nodesChangedIn1.insert(mappingsToRemoveIt->first);
	nodesChangedIn2.insert(mappingsToRemoveIt->second);
	verifiedAlignment.erase(*mappingsToRemoveIt);
      }
    }
  }

  static inline void verifyPair(mappingType firstMap, mappingType secondMap,
				alignmentType & verifiedAlignment,
				DAGraph & g1, DAGraph & g2,
				boost::numeric::ublas::matrix<float> & simMat,
				removalListType & removalList,
				set<unsigned> & nodesChangedIn1,
				set<unsigned> & nodesChangedIn2,
				bool & alignmentVerified,
				int verificationMode) {
    // Maps are (e1, e1_p) and (e2, e2_p)
    const unsigned e1 = firstMap.first;
    const unsigned e2 = secondMap.first;
    const unsigned e1_p = firstMap.second;
    const unsigned e2_p = secondMap.second;
    
    if ((e1 == e2) && (e1_p == e2_p)) {
      return;
    }
    bool conflictFound = false;
    
    // Check to see if the same element is mapped twice.  If yes, remove mapping with
    // lower similarity
    if ((e1 == e2) || (e1_p == e2_p)) {
      conflictFound = true;
    } else {
      /*
      const bool e2_in_e1 = g1.isDescendent(e2,e1);
      const bool e1_in_e2 = g1.isDescendent(e1,e2);
      const bool e2_p_in_e1_p = g2.isDescendent(e2_p,e1_p);
      const bool e1_p_in_e2_p = g2.isDescendent(e1_p,e2_p);
      */
      
      //  cout << "e2_in_e1: " << e2_in_e1 << endl;
      //  cout << "e1_in_e2: " << e1_in_e2 << endl;
      //  cout << "e2_p_in_e1_p: " << e2_p_in_e1_p << endl;
      //  cout << "e1_p_in_e2_p: " << e1_p_in_e2_p << endl;
      
      
      
      if (verificationMode == 1) {
	// Check for parent-child criss cross (isDescendent(possibleDescendent, possibleAncestor))
	if (g1.isDescendent(e2,e1) && g2.isDescendent(e1_p,e2_p)) {
	  conflictFound = true;
	  //cout << "criss cross" << endl;
	  
	} else if (g1.isDescendent(e1,e2) && g2.isDescendent(e2_p,e1_p)) {
	  conflictFound = true;
	  //cout << "criss cross" << endl;
	}
      } else if (verificationMode == 2) {
	// Check for subsumption incompleteness (except for genes, which is too strict)
	if (!g1.isGene(e1) && !g2.isGene(e1_p) && !g1.isGene(e2) && !g2.isGene(e2_p)) {
	  if ((g1.isDescendent(e2,e1) != g2.isDescendent(e2_p,e1_p)) || 
	      (g1.isDescendent(e1,e2) != g2.isDescendent(e1_p,e2_p))) {
	    conflictFound = true;
	    //cout << "Subsumption incompleteness" << endl;
	  }
	}
      } else if (verificationMode == 3) {
	// Check for parent-child criss cross (isDescendent(possibleDescendent, possibleAncestor))
	if (g1.isDescendent(e2,e1) && g2.isDescendent(e1_p,e2_p)) {
	  conflictFound = true;
	  //cout << "criss cross" << endl;
	  
	} else if (g1.isDescendent(e1,e2) && g2.isDescendent(e2_p,e1_p)) {
	  conflictFound = true;
	  //cout << "criss cross" << endl;
	}

	// Check for sibling-sibling to parent-child mapping
	if (g1.areSiblings(e1,e2) && (g2.isDescendent(e1_p,e2_p) || g2.isDescendent(e2_p,e1_p))) {
	  conflictFound = true;
	  //cout << "sibling-sibling to parent-child" << endl;
	} else if (g2.areSiblings(e1_p,e2_p) && (g1.isDescendent(e1,e2) || g1.isDescendent(e2,e1))) {
	  conflictFound = true;
	  //cout << "sibling-sibling to parent-child" << endl;
	}
      }
    }
      
    
    if (conflictFound) {
      const float firstSimVal = simMat(e1,e1_p);
      const float secondSimVal = simMat(e2,e2_p);
      
      if (firstSimVal > secondSimVal) {
	removeConflict(secondMap, firstMap, verifiedAlignment, removalList,
		       nodesChangedIn1, nodesChangedIn2);
	alignmentVerified = false;
      } else if (firstSimVal < secondSimVal) {
	removeConflict(firstMap, secondMap, verifiedAlignment, removalList,
		       nodesChangedIn1, nodesChangedIn2);
	alignmentVerified = false;
      }
    }
  }
  
  static inline bool verifyPair(mappingType firstMap, mappingType secondMap,
				DAGraph & g1, DAGraph & g2,
				int verificationMode) {
    // Maps are (e1, e1_p) and (e2, e2_p)
    const unsigned e1 = firstMap.first;
    const unsigned e2 = secondMap.first;
    const unsigned e1_p = firstMap.second;
    const unsigned e2_p = secondMap.second;
    
    // Check to see if the same element is mapped twice.
    if ((e1 == e2) || (e1_p == e2_p)) {
      return true;
    }

    if (verificationMode == 1) {    
      // Check for parent-child criss cross (isDescendent(possibleDescendent, possibleAncestor))
      if (g1.isDescendent(e2,e1) && g2.isDescendent(e1_p,e2_p)) {
	return true;
	//cout << "criss cross" << endl;
      } else if (g1.isDescendent(e1,e2) && g2.isDescendent(e2_p,e1_p)) {
	return true;
	//cout << "criss cross" << endl;
      } 
    } else if (verificationMode == 2) {
      // Check for subsumption incompleteness (except for genes, which is too strict)
      if (!g1.isGene(e1) && !g2.isGene(e1_p) && !g1.isGene(e2) && !g2.isGene(e2_p)) {
	if ((g1.isDescendent(e2,e1) != g2.isDescendent(e2_p,e1_p)) || 
	    (g1.isDescendent(e1,e2) != g2.isDescendent(e1_p,e2_p))) {
	  return true;
	  //cout << "Subsumption incompleteness" << endl;
	}
      }
    } else if (verificationMode == 3) {
      // Check for parent-child criss cross (isDescendent(possibleDescendent, possibleAncestor))
      if (g1.isDescendent(e2,e1) && g2.isDescendent(e1_p,e2_p)) {
	return true;
	//cout << "criss cross" << endl;
 
      } else if (g1.isDescendent(e1,e2) && g2.isDescendent(e2_p,e1_p)) {
	return true;
	//cout << "criss cross" << endl;
      }
      
      // Check for sibling-sibling to parent-child mapping
      if (!g1.isGene(e1) && !g1.isGene(e2) && g1.areSiblings(e1,e2) && !g1.isDescendent(e1,e2) && !g1.isDescendent(e2,e1) && 
	  (g2.isDescendent(e1_p,e2_p) || g2.isDescendent(e2_p,e1_p))) {
	
	//cout << "sibling-sibling to parent-child found for (" << g1.getName(e1) << "," ;
	//cout << g2.getName(e1_p) << ") and (" << g1.getName(e2) << "," << g2.getName(e2_p) << ")" << endl;
	return true;
      } else if (!g2.isGene(e1_p) && !g2.isGene(e2_p) && g2.areSiblings(e1_p,e2_p) && !g2.isDescendent(e1_p,e2_p) && !g2.isDescendent(e2_p,e1_p) &&
		 (g1.isDescendent(e1,e2) || g1.isDescendent(e2,e1))) {

	//cout << "sibling-sibling to parent-child found for (" << g1.getName(e1) << ",";
	//cout << g2.getName(e1_p) << ") and (" << g1.getName(e2) << "," << g2.getName(e2_p) << ")" << endl;
	return true;
      }
    }
    
    return false;
  }

};

#endif // ONTOLOGY_MATCHING_ALIGNMENTUTILS
