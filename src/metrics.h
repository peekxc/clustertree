#ifndef METRICS_H
#define METRICS_H

#include <ANN/ANN.h> // ANN typedefs
#include <algorithm>  // std::transform
#include "RcppHeader.h"

/*  metrics.h
*  To facilitate the use of arbitrary metrics in the dual tree-associated algorithms, functors that
*  are able to compute a metric distance are recorded here. These effectively replace the use of 'annDist' and
*  it's associated macros for the ANN library. Each metric class is derived from an abstract metric class which
*  creates the interface that every metric functor is required to have, namely:
*  1) An overloaded function call 'operator()' which recieves as input two ANNpoint's and returns a the distance
*  between them (of type ANNdist)
*  2) An overloaded function call 'operator()' which recieves two ANNpoint's, an index 'i', and a distance 'dist', and
*  returns the incremented distance between 'dist' and the ith component/dimension between the two points. This is required
*  to facilitate incremental distance calculations and dramatically reduce the overall number of FLOPs in cases where
*  the data is sparse or high-dimensional. Assumes each component distance is non-negative.
*  3) A 'finalize' method which transforms a vector of metric distances into their 'final' quantity. This is useful for
*  metrics where it there exists a monotonic relationship between an intermediate distance value and the final distance
*  value, and the intermediate distance is used to reduce computation. Examples of this include reserving the
*  square root calculation until all distances have been computed for euclidean distance, or reserving a unit-conversion
*  operation for previously normalized coordinate distances (e.g. converting from normalized to WGS84 coordinates). If
*  not needed, can simply inherit the parent metric class and do nothing.
*
*  To allow dynamically choosing the metric from a calling function, a static function is given that returns a Metric&
*  type. Note that this implies two requirements:
*  1) There must be a default constructor available for every metric class. It is up to the calling function to
*  initialize any auxiliary member variables needed by the distance function (operator()). Optionally, a more organized
*  way would be to add an 'init' member function.
*  2) All metrics must be handled BY REFERENCE to avoid slicing. It is not required that a derived metric have the
*  same size as the parent base class.
*/

// Generic metric interface. This abstract class may be inherited by, and should be compatible with, any valid metric
// satisfying the normal properties (symmetry, non-negativity, triangle inequality, etc.). Note that by definition,
// the metric need only be defined for every pair of elements in a given set.
struct Metric {
protected:
  ANNpoint x_q, x_r;
  const ANNpointArray Q, R; // query and reference array pointers
public:
  Metric(ANNpointArray const Q_, ANNpointArray const R_) : Q(Q_), R(R_) { }
  virtual inline ANNdist distance(const int q_idx, const int r_idx) = 0;
  virtual inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist) = 0;
  virtual inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold) = 0;
  virtual inline void finalize(std::vector<ANNdist>& dist_x) = 0;
  virtual inline ANNdist finalize(ANNdist dist_x) = 0;
};

// Extended interface for Metrics which may be embedded in the Cartesian plane, i.e. the distance between
// any two points may be readily computed without knowing these points ahead of time.
struct MinkowskiMetric : public Metric {
  using Metric::distance;
  MinkowskiMetric(ANNpointArray const Q_, ANNpointArray const R_) : Metric(Q_, R_) { }
  virtual inline ANNdist distance(ANNpoint qx, ANNpoint rx) = 0;
};

// Available metrics
enum METRIC_TYPES { L2 = 0, L1 = 1, LINF = 2, LP = 3, RSL_M = 4 };

// ----------------------------------------------------------
// --- Functors to compute some notion of metric distance ---
// ----------------------------------------------------------

// L_2 norm (squared)
struct L_2 : MinkowskiMetric {
  const int d; // dimension of the space
  L_2(ANNpointArray Q, ANNpointArray R, const int dim) : MinkowskiMetric(Q, R), d(dim){ }
  inline ANNdist distance(const int q_idx, const int r_idx){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      dist += ANN_POW(*x_q - *x_r); // compute length and adv coordinate
    }
    return dist;
  }
  // Allow dimension-specific distances to be added to given distance by reference
  inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist){
    x_q = Q[q_idx], x_r = R[r_idx];
    dist += ANN_POW(x_q[d_i] - x_r[d_i]); // compute length and adv coordinate
  }

  // Allow calculation of distances subject to an upper-bound threshold. If at any point (dimension) in the distance
  // calculation the distance exceeds this threshold, INF is returned.
  inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      if ((dist += ANN_POW(*x_q - *x_r)) > threshold) { return(ANN_DIST_INF); }
    }
    return(dist);
  }

  // Finalize by taking the square root, changes dist_x by reference
  inline void finalize(std::vector<ANNdist>& dist_x){
    std::transform(dist_x.begin(), dist_x.end(), dist_x.begin(), (double(*)(double)) sqrt);
  }
  inline ANNdist finalize(ANNdist dist_x) { return(std::sqrt(dist_x)); };

  inline ANNdist distance(ANNpoint qx, ANNpoint rx){
    ANNdist dist(0.);
    for (int i = d; i > 0; --i, ++qx, ++rx) {
      dist += ANN_POW(*qx - *rx); // compute length and adv coordinate
    }
    return(dist);
  }
};

// L_1 norm
struct L_1 : MinkowskiMetric {
  const int d; // dimension of the space
  L_1(ANNpointArray Q, ANNpointArray R, const int dim) : MinkowskiMetric(Q, R), d(dim){ }
  inline ANNdist distance(const int q_idx, const int r_idx){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      dist += std::abs(*x_q - *x_r);
    }
    return dist;
  }
  // Allow dimension-specific distances to be added to given distance by reference
  inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist){
    x_q = Q[q_idx], x_r = R[r_idx];
    dist += std::abs(x_q[d_i] - x_r[d_i]); // compute length and adv coordinate
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      if ((dist += std::abs(*x_q - *x_r)) > threshold) { return(ANN_DIST_INF); }
    }
    return(dist);
  }

  // Do nothing
  inline void finalize(std::vector<ANNdist>& dist_x){}
  inline ANNdist finalize(ANNdist dist_x) { return(dist_x); };

  inline ANNdist distance(ANNpoint qx, ANNpoint rx){
    ANNdist dist(0.);
    for (int i = d; i > 0; --i, ++qx, ++rx) {
      dist += std::abs(*qx - *rx); // compute length and adv coordinate
    }
    return(dist);
  }
};

// L_inf norm
struct L_inf : MinkowskiMetric {
  const int d; // dimension of the space
  L_inf(ANNpointArray Q, ANNpointArray R, const int dim) : MinkowskiMetric(Q, R), d(dim) { }
  inline ANNdist distance(const int q_idx, const int r_idx){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      dist = std::max(dist, std::abs(*x_q - *x_r));
    }
    return dist;
  }
  // Allow dimension-specific distances to be added to given distance by reference
  inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist){
    x_q = Q[q_idx], x_r = R[r_idx];
    dist = std::max(dist, std::abs(x_q[d_i] - x_r[d_i])); // compute length and adv coordinate
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      if ((dist = std::max(dist, std::abs(*x_q - *x_r)))) { return(ANN_DIST_INF); }
    }
    return(dist);
  }

  // Do nothing
  inline void finalize(std::vector<ANNdist>& dist_x){}
  inline ANNdist finalize(ANNdist dist_x) { return(dist_x); };

  inline ANNdist distance(ANNpoint qx, ANNpoint rx){
    ANNdist dist(0.);
    for (int i = d; i > 0; --i, ++qx, ++rx) {
      dist += std::max(dist, std::abs(*qx - *rx)); // compute length and adv coordinate
    }
    return(dist);
  }
};


// Lp norm
struct L_p : MinkowskiMetric {
  const int d; // dimension of the space
  const int p; // p-norm
  L_p(ANNpointArray Q, ANNpointArray R, const int dim, const int p_norm) : MinkowskiMetric(Q, R), d(dim), p(p_norm){ }
  inline ANNdist distance(const int q_idx, const int r_idx){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      dist += std::pow(std::abs(*x_q - *x_r), p);
    }
    return dist;
  }
  // Allow dimension-specific distances to be added to given distance by reference
  inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist){
    x_q = Q[q_idx], x_r = R[r_idx];
    dist += std::pow(std::abs(x_q[d_i] - x_r[d_i]), p); // compute length and adv coordinate
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      if ((dist += std::pow(std::abs(*x_q - *x_r), p)) > threshold){
        return(ANN_DIST_INF); // return infinity if distance is above given threshold
      }
    }
    return(dist);
  }

  // Raise to the 1/p power
  inline void finalize(std::vector<ANNdist>& dist_x){
    std::transform(dist_x.begin(), dist_x.end(), dist_x.begin(), [this](ANNdist dif) { return std::pow(dif, 1.0/p); });
  }
  inline ANNdist finalize(ANNdist dist_x) { return(std::pow(dist_x, 1.0/p)); };

  inline ANNdist distance(ANNpoint qx, ANNpoint rx){
    ANNdist dist(0.);
    for (int i = d; i > 0; --i, ++qx, ++rx) {
      dist += std::pow(std::abs(*qx - *rx), p); // compute length and adv coordinate
    }
    return(dist);
  }
};


// Robust Single Linkage distance
struct RSL : MinkowskiMetric {
  const int d; // dimension of the space
  const ANNdist alpha; // alpha
  std::vector<ANNdist> r_k; // radii of balls centered at each x_i containing k points (inclusive of x_i itself)
  RSL(ANNpointArray Q, ANNpointArray R, const int dim, const ANNdist alpha_, const std::vector<ANNdist>& r_k_) : MinkowskiMetric(Q, R), d(dim), alpha(alpha_), r_k(r_k_){ }
  inline ANNdist distance(const int q_idx, const int r_idx){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      dist += ANN_POW(*x_q - *x_r);
    }
    return std::max(dist/alpha, std::max(r_k[q_idx], r_k[r_idx]));
  }
  // Allow dimension-specific distances to be added to given distance by reference
  inline void distance(const int q_idx, const int r_idx, const int d_i, ANNdist& dist){
    // x_q = Q[q_idx], x_r = R[r_idx];
    // dist += ANN_POW(((*x_q - *x_r)/alpha)); // compute length and adv coordinate
    // dist = std::max(dist, std::max(r_k[q_idx], r_k[r_idx]))
  }
  // Allow dimension-specific distance calculation for incremental distance updates
  inline ANNdist distance(const int q_idx, const int r_idx, ANNdist threshold){
    ANNdist dist(0.);
    x_q = Q[q_idx], x_r = R[r_idx];
    for (int i = d; i > 0; --i, ++x_q, ++x_r) {
      if ((dist += std::max(ANN_POW(*x_q - *x_r)/alpha, std::max(r_k[q_idx], r_k[r_idx]))) > threshold) { return(ANN_DIST_INF); }
    }
    return(dist);
  }

  // Finalize by taking the square root, changes dist_x by reference
  inline void finalize(std::vector<ANNdist>& dist_x){
    std::transform(dist_x.begin(), dist_x.end(), dist_x.begin(), (double(*)(double)) sqrt);
  }
  inline ANNdist finalize(ANNdist dist_x) { return(std::sqrt(dist_x)); };

  // Not technically correct, but this distance is used by the box computations, and
  // will always be smaller or equal to the mutual reachability distance, so should be safe to use
  inline ANNdist distance(ANNpoint qx, ANNpoint rx){
    ANNdist dist(0.);
    for (int i = d; i > 0; --i, ++qx, ++rx) {
      dist += ANN_POW(*qx - *rx); // compute length and adv coordinate
    }
    return(dist);
  }
};


// Assume the error handling is done in R
SEXP chooseMetric_int(std::string metric_name);
Metric& getMetric(SEXP metric_ptr);

#endif