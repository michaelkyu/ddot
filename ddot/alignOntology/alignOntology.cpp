#include <iostream>
#include <time.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <stdlib.h>
#include "graph.h"
#include "alignmentUtils.h"

int main(int argc, char* argv[]) {
  if (argc < 5) {
    cout << "Needs 4 arguments - ontologyComputed, referenceOntology, minimum similarity value accepted, version of semantic verification (options are 'criss_cross', 'strict_hierarchy', 'sib_sib', and 'none')" << endl;
    cout << "Optional 5th argument - allow multiple mappings of same node when scores are identical.  1 will allow multiple mappings, 0 will not (default is 0 if no 5th argument is given" << endl;
    cout << "Optional 6th argument - treat genes as terms.  Default is 1 (yes).  0 is no." << endl;
    cout << "Optional 7th argument - category name of terminal nodes (i.e. gene, patient, etc.). Default is gene" << endl;
    cout << "similarity value can range from 0 to 1 (1 means nodes in ontology are identical)" << endl;
    return 0;
  }
  map<string, unsigned> geneNamesToIDs;
  string ontology1 = argv[1];
  string ontology2 = argv[2];
  string terminalName = "gene";
  if (argc >= 8) {
    terminalName = argv[7];
  }
  cout << "Getting ontology 1" << endl;
  DAGraph g1(ontology1, geneNamesToIDs, true, terminalName);
  cout << "Getting ontology 2" << endl;
  DAGraph g2(ontology2, geneNamesToIDs, true, terminalName);

  bool allowMultiples = false;
  if ((argc >= 6) && (atoi(argv[5]) == 1)) {
    allowMultiples = true;
  }

  bool genesAsTerms = true;
  if ((argc >=7) && (atoi(argv[6]) == 0)) {
    genesAsTerms = false;
  }

  vector<string> geneIDsToNames(geneNamesToIDs.size());
  for (map<string, unsigned>::iterator it = geneNamesToIDs.begin(); it != geneNamesToIDs.end(); ++it) {
    geneIDsToNames[it->second] = it->first;
  }

  float threshold = 0.01;
  float minSimilarityForPrealignment = atof(argv[3]);
  cout << "Minimum similarity value accepted: " << minSimilarityForPrealignment << endl;

  int verificationMode = 0;
  string verMode = argv[4];
  if (verMode == "criss_cross") {
    cout << "No parent child criss cross allowed" << endl;
    verificationMode = 1;
  } else if (verMode == "strict_hierarchy") {
    cout << "Strict hierarchy maintained in mappings" << endl;
    verificationMode = 2;
  } else if (verMode == "none") {
    cout << "No hierarchy information enforced in mappings" << endl;
  } else if (verMode == "sib_sib") {
    cout << "No parent child criss cross allowed and no sibling-sibling to parent-child mapping allowed" << endl;
    verificationMode = 3;
  } else {
    cout << "Options for semantic verification (argument 4) are 'criss_cross', 'strict_hierarchy', 'sib_sib', and 'none'" << endl;
    return 1;
  }
  

  cout << "Terms in first ontology: " << g1.numNodes() << "\t" << "Nodes in second ontology: " << g2.numNodes() << endl;
  boost::numeric::ublas::matrix<float> extSim(g1.numNodes(),g2.numNodes());
  cout << "Calculating extsim" << endl;
  time_t start, end;
  time (&start);
  alignmentUtils::calculateExtSim(extSim, g1, g2);
  time (&end);
  double dif = difftime(end,start);
  cout << "calculateExtSim took " << dif << " seconds" << endl;

  boost::numeric::ublas::matrix<float> lastSim = extSim;
  
  //removalListType removalList(g1.numNodes(), g2.numNodes());
  vector< alignmentType > verifiedAlignments;
  alignmentType finalAlignment;
  //boost::numeric::ublas::matrix<bool> foreverRemoved(g1.numNodes(), g2.numNodes());

  // Mike Yu edit: limit the number of iterations
  int max_iter = 10;
  int curr_iter = 1;
  cout << "Max iterations: " << max_iter << endl;

  bool repeatAlignmentFound = false;
  while (!repeatAlignmentFound) {

    // Mike Yu: print out iteration number
    //cout << "Starting iteration" << endl;   
    cout << "Starting iteration " << curr_iter << endl;   

    boost::numeric::ublas::matrix<float> simMat(extSim);
    cout << "Adding relSim to simMat" << endl;

    time (&start);
    alignmentUtils::addRelSimToSimMat(simMat, lastSim, g1, g2, genesAsTerms);
    time (&end);
    double dif = difftime(end,start);
    cout << "addRelSimToSimMat took " << dif << " seconds" << endl;

    //for (int i = 0; i < relSim.size1(); ++i) {
    //  cout << g1.getName(i) << "\t" << extSim(i,i) << "\t" << relSim(i,i) << endl;
    //}
    //return 0;
    
    alignmentType alignment;
    time(&start);
    cout << "Aligning" << endl;
    alignmentUtils::getAlignment(simMat,alignment,threshold,
				 g1,g2,minSimilarityForPrealignment, verificationMode, allowMultiples);
    time (&end);
    dif = difftime(end,start);
    cout << "Alignment took " << dif << " seconds" << endl;
    cout << "Got an alignment" << endl;
    
    for (vector<alignmentType>::iterator it = verifiedAlignments.begin();
	 it != verifiedAlignments.end(); ++it) {
      if (alignment == *it) {
	repeatAlignmentFound = true;
	finalAlignment = alignment;
	cout << "Matched" << endl;
      }
    }
    lastSim = simMat;

    /*
    for (int i = 0; i < lastSim.size1(); ++i) {
      for (int j = 0; j < lastSim.size2(); ++j) {
	if (removalList.isInRemovalList(make_pair(i,j))) {
	  foreverRemoved(i,j) = true;
	}
      }
    }
    */
    verifiedAlignments.push_back(alignment);
    
    // End early
    //repeatAlignmentFound = true;
    //finalAlignment = verifiedAlignment;
    //cout << "Matched" << endl;

    // Mike Yu edit: limit the number of iterations
    curr_iter += 1;
    if (curr_iter > max_iter) {
      repeatAlignmentFound = true;
      finalAlignment = alignment;
      cout << "Matched" << endl;
    }
  }

  for (alignmentType::iterator it = finalAlignment.begin(); it != finalAlignment.end(); ++it) {
    
    if (!g1.isGene(it->first)) {
      cout << g1.getName(it->first) << "\t" << g2.getName(it->second) << "\t" << lastSim(it->first,it->second) << endl;
    }

    /* THE BELOW CODE WILL PRINT DESCENDENTS OF EACH TERM
    cout << g1.getName(it->first) << "\t" << g2.getName(it->second) << "\t" << lastSim(it->first,it->second) << "\t";    
    for (set<unsigned>::iterator genesIt = g1.getGenesBegin(it->first); genesIt != g1.getGenesEnd(it->first); ++genesIt) {
      cout << geneIDsToNames[*genesIt] << ",";
    }
    cout << "\t";
    for (set<unsigned>::iterator genesIt = g2.getGenesBegin(it->second); genesIt != g2.getGenesEnd(it->second); ++genesIt) {
      cout << geneIDsToNames[*genesIt] << ",";
    }
    cout << endl;
    */
  }
 
  return 1;
}
