#ifndef METRICS_H
#define METRICS_H

#include <ANN/ANN.h> // ANN typedefs

// Functors to compute some notion of distance
namespace Metrics{
  // L_2 norm (squared)
  struct L_2 {
    const int d;
    const bool use_inc_distance;
    L_2(const int _d) : d(_d), use_inc_distance(false) { }
    inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
      ANNdist dist = 0;
      for (int i = 0; i < d; ++i) { dist += (x_j[i] - x_i[i]) * (x_j[i] - x_i[i]); }
      return dist;
    }
    // Allow dimension-specific distance calculation for incremental distance updates
    inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i){
      return (x_j[i] - x_i[i]) * (x_j[i] - x_i[i]);
    }
  };
}


#endif