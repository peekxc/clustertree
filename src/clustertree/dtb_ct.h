#ifndef DTB_CT_H
#define DTB_CT_H

#include <Rcpp.h>
using namespace Rcpp;

#include <DT/dt.h> // Dual Tree definitions
#include <DT/MST/dtb.h> // Dual Tree Borvuka definitions
#include <metrics.h> // metric definitions

class DTB_CT : DualTreeBoruvka {
  NumericVector r_k; // contains vector of radii
  double alpha; // extends the connection radius
  DTB_CT(const bool prune, const int dim, const int n, Metric& m, NumericVector& _r_k, const double _alpha);

  // Override the istance calculation to bound the connection radii
  ANNdist computeDistance(const int q_idx, const int r_idx,
                                  ANNdist eps1 = ANN_DIST_INF,
                                  ANNdist eps2 = ANN_DIST_INF) override;
};


#endif