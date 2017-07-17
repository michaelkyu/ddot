#include <iostream>
#include <time.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>
#include "graph.h"

bool compDoubleSortDescending(double i, double j) { return i > j; }


/*unsigned getBin(unsigned numGenesInThis, unsigned numGenesInGraph) {
  //return numGenesInThis;
  
  if (numGenesInThis <= 4) {
    return numGenesInThis;
  } else if ((numGenesInThis >= 5) && (numGenesInThis <= 8)) {
    return 5;
  } else {
    return 6;
  }
  
  }*/

unsigned getBin(unsigned numGenesInThis, unsigned numGenesInGraph) {
  //return numGenesInThis;
  
  if (numGenesInThis <= 4) {
    return numGenesInThis;
  } else if ((numGenesInThis >= 5) && (numGenesInThis <= 8)) {
    return 5;
  } else if ((numGenesInThis >= 9) && (numGenesInThis <= 16)) {
    return 6;
  } else if ((numGenesInThis >= 17) && (numGenesInThis <= 32)) {
    return 7;
  } else if ((numGenesInThis >= 33) && (numGenesInThis <= 64)) {
    return 8;
  } else if ((numGenesInThis >= 65) && (numGenesInThis <= 999)) {
    return 9;
  } else {
    return (9 + (numGenesInThis / 1000));
  }
}

/*unsigned getBin(unsigned numGenesInThis, unsigned numGenesInGraph) {
  if (numGenesInThis <= 4) {
    return numGenesInThis;
  }
  unsigned bin = 4;
  unsigned max_bin_size = 4;
  unsigned next_addition = max_bin_size;
  unsigned halfway = numGenesInGraph / 2;
  while (numGenesInThis > max_bin_size) {
    //cout << max_bin_size << endl;
    if ((max_bin_size < halfway) && ((max_bin_size + next_addition) < halfway)) {
      max_bin_size += next_addition;
      next_addition = max_bin_size;
    } else if ((max_bin_size < halfway) && ((max_bin_size + next_addition) >= halfway)){
      max_bin_size = halfway + (halfway - max_bin_size);
      next_addition /= 2;
    } else if ((max_bin_size > halfway) && ((numGenesInGraph - 3) > max_bin_size) &&
	       ((numGenesInGraph - 3) <= (max_bin_size) + next_addition)) {
      max_bin_size = numGenesInGraph - 3;
      next_addition = 1;
    } else if ((max_bin_size > halfway) && ((numGenesInGraph - 3) > (max_bin_size + next_addition))) {
      max_bin_size += next_addition;
      next_addition /= 2;
    } else if ((max_bin_size > halfway) && ((numGenesInGraph - 3) <= max_bin_size)) {
      max_bin_size += next_addition;
      next_addition = 1;
    }
    ++bin;
  }
  return bin;
  }*/


pair<int,int> binIDtoRange(int ID, int numGenesInGraph) {
  pair<int, int> retval;
  bool minFound = false;
  bool maxFound = false;
  for (int i = 0; i <= numGenesInGraph; ++i) {
    if (!minFound && getBin(i, numGenesInGraph) == ID) {
      retval.first = i;
      minFound = true;
    } else if (!maxFound && getBin(i, numGenesInGraph) > ID) {
      retval.second = i-1;
      maxFound = true;
      return retval;
    }
  }
  retval.second = numGenesInGraph;
  maxFound = true;
  return retval;
}
    
