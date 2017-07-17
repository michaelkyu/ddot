#include "calculateFDRsForAlignment.h"

int main(int argc, char* argv[]) {

  if (argc < 5) {
    cout << "Needs 4 arguments - computed ontology, alignment, root for rand alignment files, num rand alignment files" << endl;
    //cout << "Optional 5th argument - False Discovery Rate.  Won't print mappings with higher FDR" << endl;
    cout << "Optional 5th argument - Effect size.  Default is 0.  Otherwise will decrease actual score by effect size to ensure effect is not only statistically significant but also large" << endl;
    cout << "Optional 6th argument - Using pvals instead of scores (so lower is better).  1 = using pvals.  Anything else = using scores" << endl;
    return 0;
  }
  map<string, unsigned> geneNamesToIDs;
  string ontology1 = argv[1];
  DAGraph g1(ontology1, geneNamesToIDs, false);
  string alignmentFile = argv[2];
  string randAlignmentFiles = argv[3];
  int numRandAlignmentFiles = atoi(argv[4]);

  /*
  double FDR = 2;
  if (argc >= 6) {
    FDR = atof(argv[5]);
  }
  */

  bool useLowerInstead = false;
  double effectSize = 0;
  if (argc >= 6) {
    effectSize = atof(argv[5]);
    if (argc >= 7) {
      if (atoi(argv[6]) == 1) {
	useLowerInstead = true;
      }
    }
  }

  int numGenes = geneNamesToIDs.size();
  unsigned numBins = getBin(numGenes,numGenes)+1;

  // Count the number of nodes in the ontology of each size
  vector<int> numNodesOfSize(numBins,0);
  for (vector<Node>::iterator nodeIt = g1.nodesBegin(); nodeIt != g1.nodesEnd(); ++nodeIt) {
    if (!nodeIt->isGene()) {
      ++numNodesOfSize[getBin(nodeIt->numGenes(),numGenes)];
    }
  }

  vector<double> empty;
  vector<vector<double> > nodeSizeToScoresRand(numBins, empty);

  for (int i = 0; i < numRandAlignmentFiles; ++i) {
    string line;
    string fileToOpen = randAlignmentFiles + "_" + boost::lexical_cast<string>( i );
    ifstream file(fileToOpen.c_str());
    int lastID = -1;
    if (file.is_open()) {
      while (file.good()) {
	getline(file,line);
	vector<string> tokens;
	Utils::Tokenize(line, tokens, "\t");
	if (tokens.size() == 3) {
	  if (!g1.isGene(g1.getID(tokens[0]))) {
	    int newNodeID = g1.getID(tokens[0]);
	    if (newNodeID != lastID) {
	      lastID = newNodeID;
	      nodeSizeToScoresRand[getBin(g1.numGenesInNode(newNodeID),numGenes)].push_back(atof(tokens[2].c_str()));
	    }
	  }
	}
      }
    }
    file.close();
  }

  for (unsigned i = 0; i < numBins; ++i) {
    if (numNodesOfSize[i]) {
      if (useLowerInstead) {
	sort(nodeSizeToScoresRand[i].begin(), nodeSizeToScoresRand[i].end());
      } else {
	sort(nodeSizeToScoresRand[i].begin(), nodeSizeToScoresRand[i].end(), compDoubleSortDescending);
      }
      
      /*
      cout << "BEGIN " << i << endl;
      for (vector<double>::iterator it = nodeSizeToScoresRand[i].begin(); it != nodeSizeToScoresRand[i].end(); ++it) {
	cout << *it << endl;
      }
      */
      //cout << i << "\t" << nodeSizeToScoresRand[i].size() / (static_cast<double>(numRandAlignmentFiles)*numNodesOfSize[i]) << "\t" << numNodesOfSize[i] << "\t" << nodeSizeToScoresRand[i].size() / static_cast<double>(numRandAlignmentFiles) << endl;
    }
  }

  vector<vector<double> > nodeSizeToScoresActual(numBins, empty);

  string line;
  ifstream file(alignmentFile.c_str());
  if (file.is_open()) {
    int lastID = -1;
    while (file.good()) {
      getline(file,line);
      vector<string> tokens;
      Utils::Tokenize(line, tokens, "\t");
      if (tokens.size() == 3) {
	if (!g1.isGene(g1.getID(tokens[0]))) {
	  int newNodeID = g1.getID(tokens[0]);
	  if (newNodeID != lastID) {
	    lastID = newNodeID;
	    //nodeSizeToScoresActual[getBin(g1.numGenesInNode(newNodeID),numGenes)].push_back(atof(tokens[2].c_str()));
	    nodeSizeToScoresActual[getBin(g1.numGenesInNode(newNodeID),numGenes)].push_back(atof(tokens[2].c_str()) - effectSize);
	    //cout << "Adding " << getBin(g1.numGenesInNode(newNodeID),numGenes) << "\t" << atof(tokens[2].c_str()) << endl;
	  }
	}
      }
    }
  }
  file.close();

  for (unsigned i = 0; i < numBins; ++i) {
    if (numNodesOfSize[i]) {
      if (useLowerInstead) {
	sort(nodeSizeToScoresActual[i].begin(), nodeSizeToScoresActual[i].end());
      } else {
	sort(nodeSizeToScoresActual[i].begin(), nodeSizeToScoresActual[i].end(), compDoubleSortDescending);
      }
    }
  }

  map<double,double> emptyMap;
  vector<map<double, double> > nodeSizeAndScoreToFDR(numBins, emptyMap);
  
  for (unsigned i = 0; i < numBins; ++i) {
    if (numNodesOfSize[i]) {
      vector<double>::iterator randIt = nodeSizeToScoresRand[i].begin();
      unsigned numRandBetter = 0;
      unsigned numActualBetter = 0;
      for (vector<double>::iterator actualIt = nodeSizeToScoresActual[i].begin();
	   actualIt != nodeSizeToScoresActual[i].end(); ++actualIt) {
	++numActualBetter;
	while ((randIt != nodeSizeToScoresRand[i].end()) && 
	       (((!useLowerInstead) && (*randIt > *actualIt)) || ((useLowerInstead) && (*randIt < *actualIt)))) {
	  ++numRandBetter;
	  ++randIt;
	}
	double FDR;
	if ((numActualBetter == 0) && (numRandBetter != 0)) {
	  FDR = 1;
	} else if ((numActualBetter == 0) && (numRandBetter == 0)) {
	  FDR = -1;
	} else {
	  FDR = (numRandBetter / static_cast<double>(numRandAlignmentFiles)) / numActualBetter;
	}
	nodeSizeAndScoreToFDR[i][*actualIt] = FDR;
	//cout << i << "\t" << *actualIt << "\t" << FDR << endl;
      }
    }
  }
  
  cout << "#BinID\tBin Range\tNodes in Bin";
  double maxScore = 0.91;
  for (double score = 0.05; score < maxScore; score += 0.05) {
    cout << "\t" << score;
  }
  cout << endl;
  for (int binID = 0; binID < numBins; ++binID) {
    if (numNodesOfSize[binID] != 0) {
      pair <int,int> binRange = binIDtoRange(binID, numGenes);
      cout << "#" << binID << "\t" << binRange.first << " - " << binRange.second << "\t" << numNodesOfSize[binID];
      for (double score = 0.05; score < maxScore; score += 0.05) {
	unsigned numActualBetter = 0;
	for (vector<double>::iterator actualIt = nodeSizeToScoresActual[binID].begin();
	     actualIt != nodeSizeToScoresActual[binID].end(); ++actualIt) {
	  if (((!useLowerInstead) && (score < *actualIt)) || ((useLowerInstead) && (score > *actualIt))) {
	    ++numActualBetter;
	  } else {
	    break;
	  }
	}
	
	
	unsigned numRandBetter = 0;
	for (vector<double>::iterator randIt = nodeSizeToScoresRand[binID].begin();
	     randIt != nodeSizeToScoresRand[binID].end(); ++randIt) {
	  if (((!useLowerInstead) && (score < *randIt)) || ((useLowerInstead) && (score > *randIt))) {
	    ++numRandBetter;
	  } else {
	    break;
	  }
	}
	
	double FDR;
	if ((numActualBetter == 0) && (numRandBetter != 0)) {
	  FDR = 1;
	} else if ((numActualBetter == 0) && (numRandBetter == 0)) {
	  FDR = -1;
	} else {
	  FDR = (numRandBetter / static_cast<double>(numRandAlignmentFiles)) / numActualBetter;
	}
	if (FDR == -1) {
	  cout << "\tUndef";
	} else {
	  cout << "\t" << FDR;
	}
      }
      cout << endl;
    }
  }


  cout << "#FDRs For Individual Alignments" << endl;
  ifstream file2(alignmentFile.c_str());
  if (file2.is_open()) {
    while (file2.good()) {
      getline(file2,line);
      vector<string> tokens;
      Utils::Tokenize(line, tokens, "\t");
      if (tokens.size() == 3) {
	cout << line << "\t";
	int nodeID = g1.getID(tokens[0]);
	if (!g1.isGene(nodeID)) {
	  int numGenesInNode = g1.numGenesInNode(nodeID);
	  double score = atof(tokens[2].c_str()) - effectSize;
	  cout << nodeSizeAndScoreToFDR[getBin(numGenesInNode, numGenes)][score] << "\t" << numGenesInNode << endl;
	} else {
	  cout << 0 << endl;
	}
      }
    }
  }
  file2.close();
  
  return 1;
}
