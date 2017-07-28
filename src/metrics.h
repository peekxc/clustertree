#ifndef METRICS_H
#define METRICS_H

#include <ANN/ANN.h> // ANN typedefs
#include <utilities.h> // REP

// L_p norm
struct L_2 {
  const int d;
  L_2(const int _d) : d(_d) { }
  inline ANNdist operator()(ANNpoint x_i, ANNpoint x_j){
    switch(d){
      case 1:
        return (x_j[0] - x_i[0]) * (x_j[0] - x_i[0]);
      case 2:
        return (x_j[0] - x_i[0]) * (x_j[0] - x_i[0]) + (x_j[1] - x_i[1]) * (x_j[1] - x_i[1]);
      default:
        break;
    }
    ANNdist dist = 0;
    for (int i = 0; i < d; ++i)
    { dist += (x_j[i] - x_i[i]) * (x_j[i] - x_i[i]); }
    return dist;
  }
};

#endif