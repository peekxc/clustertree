#ifndef BOUND_H
#define BOUND_H

#include <ANN/ANNx.h> // ANN type defs, ANNorthRect, etc.

// Every node has a set number of bounds related to it that allow pruning of branch queries,
// each of which should only need to be computed once. The Bound struct stores these bounds,
// allowing any recursion-based structure to be memoized.
struct Bound {
  ANNdist rho, lambda;
  ANNpoint centroid; // Centroid computed during tree construction
  Bound() : rho(-1.0), lambda(-1.0), centroid(NULL) { }
  Bound(const ANNorthRect& bb, int d) : rho(-1.0), lambda(-1.0) {
    centroid = new ANNcoord[d];
    for (int i = 0; i < d; ++i){ centroid[i] = ANNcoord((bb.lo[i] + bb.hi[i]) / 2.0); }
  }
  ~Bound(){ delete[] centroid; }
};
#endif