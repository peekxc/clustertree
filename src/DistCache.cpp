// // DistCache.cpp
// // This containers sole purpose is to provide quick access, storage, and update utilities for
// // storing distances (floating points) which are indexed by pairs of indexes (integers), but these distances are
// // assumed to be added in any order, changed frequently, and may in total be far less than then the highest index (n),
// // which is also assumed to be known upon instantiation.
// // Conceptually, this is similar to a SparseMatrix design, but modified for faster updating and specialized specifically for
// // distance 'matrices' over metric distances. Assuming a metric space allows for a linear index.
// // Author: Matt Piekenbrock
//
// #include "DistCache.h"
// #include "RcppHeader.h"
//
// template <class METRIC_T>
// DistCache<METRIC_T>::DistCache(const int size, const ANNpointArray pt_ref, METRIC_T& m) : metric(m), n(size), pts(pt_ref) {
//   idx = std::set<int>();
//   dist = std::vector<ANNdist>();
//   dist.reserve(n);// Can reasonably assume at least one distance will be computed for each point.
// }
//
// // Check for the existence of a key of a pair
// template <typename METRIC_T>
// inline bool DistCache<METRIC_T>::exists(const int from, const int to){
//   return(idx.count(INDEX_TF(n, from < to ? from : to, from < to ? to : from)) > 0);
// }
//
// // Takes a pair of indices and returns the distance between the points they correspond too.
// template <typename METRIC_T>
// inline ANNdist DistCache<METRIC_T>::operator()(const int from, const int to){
//   if (from == to){ return(0.0); }
//   const int i = from < to ? INDEX_TF(n, from, to) : INDEX_TF(n, to, from);
//   // std::set<int>::iterator el = idx.find(i);
//   std::pair<std::set<int>::iterator, bool> ret = idx.insert(i);
//   if (ret.second){ // already in set
//     return(dist[std::distance(idx.begin(), ret.first)]); // return distance at true index
//   } else { // was just added
//     dist.push_back(metric.distance(pts[from], pts[to])); // compute the and store new metric distance
//     return(dist.back());
//   }
// }
//
// template <typename METRIC_T>
// inline ANNdist DistCache<METRIC_T>::test2(const int from, const int to){
//   ANNdist pt_dist;
//   int min_idx = std::min(from, to), max_idx = std::max(from, to);
//   std::pair<int, int> edge_pair = std::make_pair(min_idx, max_idx);
//   if (dist_map.count(edge_pair) > 0){
//     pt_dist = dist_map[edge_pair];
//   }  else {
//     pt_dist = metric.distance(pts[from], pts[to]);
//     dist_map[edge_pair] = pt_dist; // save the distance
//   }
//   return(pt_dist);
// }
//
// template <typename METRIC_T>
// inline ANNdist DistCache<METRIC_T>::test3(const int from, const int to){
//   if (from == to){ return(0.0); }
//   const int i = from < to ? INDEX_TF(n, from, to) : INDEX_TF(n, to, from);
//   std::map<int, int>::iterator el = map_idx.lower_bound(i);
//   if (el != map_idx.end()){
//     // Rprintf("second: %d\n",el->second);
//     return(dist.at(el->second)); // return distance at true index
//   } else {
//     // map_idx.emplace_hint(el, i, map_idx.size()); // construct the index in constant time
//     map_idx.insert(el, std::make_pair(i, map_idx.size()));
//     dist.push_back(metric.distance(pts[from], pts[to])); // compute the and store new metric distance
//     return(dist.back());
//   }
// }
//
// template <typename METRIC_T>
// inline ANNdist DistCache<METRIC_T>::test4(const int from, const int to){
//   if (from == to){ return(0.0); }
//   const int i = from < to ? INDEX_TF(n, from, to) : INDEX_TF(n, to, from);
//   if (umap_idx.count(i) > 0){
//     // Rprintf("second: %d\n",el->second);
//     return(dist.at(umap_idx[i])); // return distance at true index
//   } else {
//     // map_idx.emplace_hint(el, i, map_idx.size()); // construct the index in constant time
//     umap_idx.insert(std::make_pair(i, map_idx.size()));
//     dist.push_back(metric.distance(pts[from], pts[to])); // compute the and store new metric distance
//     return(dist.back());
//   }
// }
//
// template <typename METRIC_T>
// inline ANNdist DistCache<METRIC_T>::test5(const int from, const int to){
//   const int min_idx = from < to ? from : to;
//   const int max_idx = from < to ? to : from;
//   if (E(min_idx, max_idx) == 0){
//     E(min_idx, max_idx) = metric.distance(pts[from], pts[to]);
//   }
//   return(E(min_idx, max_idx));
// }
//
// #include "ANN_util.h" // converting matrices to ANNpointarrays
// #include "profiler.h"
//
// // [[Rcpp::export]]
// NumericVector testDistCache(const NumericMatrix& x) {
//   std::chrono::time_point<std::chrono::high_resolution_clock> start;
//   std::chrono::time_point<std::chrono::high_resolution_clock> finish;
//   const int n = x.nrow();
//   const int d = x.ncol();
//   L_2 x_metric = L_2(d);
//   arma::mat dat;
//
//   ANNpointArray pts = matrixToANNpointArray(x);
//   DistCache< L_2 > my_cache = DistCache< L_2 >(n, pts, x_metric);
//   DistCache< L_2 > my_cache2 = DistCache< L_2 >(n, pts, x_metric);
//
//   Rcout << " -------- Testing set -------- " << std::endl;
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Creating new elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Accessing stored elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache2(i, j);
//       for (int j2 = 0; j2 < j; ++j2){
//         my_cache2(i, j2);
//       }
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Storing and accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   // Unordered_map testing
//   Rcout << " -------- Testing unordered map w/ direct pair -> dist -------- " << std::endl;
//   my_cache.dist_map = std::unordered_map<std::pair<int, int>, ANNdist, std::hash<std::pair<int, int>> >();
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test2(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Creating new elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test2(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   // BEGIN_PROFILE();
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache2.test2(i, j);
//       for (int j2 = 0; j2 < j; ++j2){
//         my_cache2.test2(i, j2);
//       }
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Storing and accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//   // END_PROFILE("Creating and accessing elements: ");
//
//
//   // ----- Map testing -----
//   Rcout << " -------- Testing map -------- " << std::endl;
//   my_cache.map_idx = std::map<int, int>();
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test3(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Creating new elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test3(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache2.test3(i, j);
//       for (int j2 = 0; j2 < j; ++j2){
//         my_cache2.test3(i, j2);
//       }
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Storing and accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   // ----- Unordered Map testing -----
//   Rcout << " -------- Testing unordered_map -------- " << std::endl;
//   my_cache.umap_idx = std::unordered_map<int, int>();
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test4(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Creating new elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test4(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache2.test4(i, j);
//       for (int j2 = 0; j2 < j; ++j2){
//         my_cache2.test4(i, j2);
//       }
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Storing and accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   // ----- Armadillo Sparse Matrix testing -----
//   Rcout << " -------- Testing Armadillo Sparse Matrix -------- " << std::endl;
//   my_cache.E = arma::sp_mat(my_cache.n, my_cache.n);
//   my_cache2.E = arma::sp_mat(my_cache.n, my_cache.n);
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test5(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Creating new elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache.test5(i, j);
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = i; j < n; ++j){
//       my_cache2.test5(i, j);
//       for (int j2 = 0; j2 < j; ++j2){
//         my_cache2.test5(i, j2);
//       }
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Storing and accessing elements: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//
//   return(NumericVector::create());
// }
//
//
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
//
// /*** R
// x <- matrix(rnorm(150))
// clustertree:::testDistCache(x)
// */
