#ifndef METRICS_H
#define METRICS_H

#include <ANN/ANN.h> // ANN typedefs

// Generic metric abstract class to allow polymorphism
struct Metric {
  virtual inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j) = 0;
  virtual inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j, const int i) = 0;
  virtual ~Metric() {}
};

// ----------------------------------------------------------
// --- Functors to compute some notion of metric distance ---
// ----------------------------------------------------------


// L_2 norm (squared)
struct L_2 : Metric {
  const unsigned int d; // dimension of the space
  L_2(const int _d) : d(_d) { }
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


#endif