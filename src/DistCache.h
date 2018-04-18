#ifndef DISTCACHE_H
#define DISTCACHE_H
// DistCache.h
// This containers sole purpose is to provide quick access, storage, and update utilities for
// storing distances (floating points) which are indexed by pairs of indexes (integers), but these distances are
// assumed to be added in any order, changed frequently, and may in total be far less than then the highest index (n),
// which is also assumed to be known upon instantiation.
// Conceptually, this is similar to a SparseMatrix design, but modified for faster updating and specialized specifically for
// distance 'matrices' over metric distances. Assuming a metric space allows for a linear index.
// Author: Matt Piekenbrock

#include "ANN.h" // ANNdist
#include "utilities.h" // Indexing macros, INDEX_TF
#include "metrics.h" // To compute metric distance
#include <type_traits> // for checking template type traits
#include <unordered_map>

// Quicker hash function to improve unordered_map access
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
template<typename S, typename T> struct std::hash<std::pair<S, T>>
{
  inline size_t operator()(const pair<S, T> & v) const
  {
    size_t seed = 0;
    ::hash_combine(seed, v.first);
    ::hash_combine(seed, v.second);
    return seed;
  }
};

template <typename METRIC_T>
struct DistCache{
  static_assert(std::is_base_of<Metric, METRIC_T>::value, "Must be a valid metric.");
  METRIC_T metric; // metric distance to use, includes dimensionality of space
  const int n; // size of the data set
  const ANNpointArray pts; // pointer to where the n points are stored
  std::set<int> idx; // max size: n x n
  std::vector<ANNdist> dist; // max size: n x n

  // arma::umat sp_idx;
  // arma::vec dist_values;
  // arma::sp_mat mat;



  // Constructor needs size of data set, pointer to data points, and metric to use
  DistCache(const int size, const ANNpointArray pt_ref, METRIC_T& m);

  // Check for the existence of a key of a pair
  inline bool exists(const int from, const int to);

  // Retrieve a cached distance, or computes and stores (caches) a new one if not available
  ANNdist operator()(const int from, const int to);

  // Testing other possibilities
  // unordered_map
  std::unordered_map<std::pair<int, int>, ANNdist, std::hash<std::pair<int, int>>> dist_map;
  ANNdist test2(const int from, const int to);

  // ordered map of lower-triangular indices --> relative sequential indices for dist vector
  std::map<int, int> map_idx;
  ANNdist test3(const int from, const int to);

  // ordered map of lower-triangular indices --> relative sequential indices for dist vector
  std::unordered_map<int, int> umap_idx;
  ANNdist test4(const int from, const int to);

  // Sparse Matrix w/ Armadillo
  arma::sp_mat E;
  ANNdist test5(const int from, const int to);


};

#endif
