#include <Rcpp.h>
using namespace Rcpp;

#include <DT/dt.h> // Dual Tree definitions
#include <DT/MST/dtb.h> // Dual Tree Borvuka definitions


class DTB_CT{
  DTB_CT(const bool prune, const int dim, Metric& m); // use default constructor

  // Override the istance calculation to bound the connection radii
  ANNdist computeDistance(const int q_idx, const int r_idx,
                                  ANNdist eps1 = ANN_DIST_INF,
                                  ANNdist eps2 = ANN_DIST_INF) override;
};