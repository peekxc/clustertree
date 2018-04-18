//
// #include "ANN.h"
// #include "ANN_util.h"
// #include "RcppHeader.h"
// #include "FastPtrHash.h"
//
// #define TEST_MAX(x,y)		((x) > (y) ? (x) : (y))
//
// struct metric {
//   const int dim;
//   const ANNpointArray pts;
//   ANNpoint q, r;
//   metric(ANNpointArray _pts, const int d) : dim(d), pts(_pts){}
//   ANNdist distance1(const int i, const int j){
//     q = pts[i], r = pts[j];
//     ANNdist dist = 0.0;
//     for (int i = dim; i > 0; --i, ++q, ++r) {
//       dist += std::max(ANN_POW(*q - *r), std::max(*q, *r)); // compute length and adv coordinate
//     }
//     return(dist);
//   }
//   ANNdist distance2(const int i, const int j){
//     q = pts[i], r = pts[j];
//     ANNdist dist = 0.0;
//     for (int i = dim; i > 0; --i, ++q, ++r) {
//       dist += TEST_MAX(ANN_POW(*q - *r), TEST_MAX(*q, *r)); // compute length and adv coordinate
//     }
//     return(dist);
//   }
// };
//
// // [[Rcpp::export]]
// NumericMatrix timesTwo(const NumericMatrix& x) {
//   ANNpointArray pt_array = matrixToANNpointArray(x);
//
//   std::unordered_map<double**, int, FashPtrHash<double*>> ptr_idx_map = std::unordered_map<double**, int, FashPtrHash<double*>>();
//   std::chrono::time_point<std::chrono::high_resolution_clock> start;
//   std::chrono::time_point<std::chrono::high_resolution_clock> finish;
//
//   const int n = x.nrow();
//   metric m = metric(pt_array, x.ncol());
//   int c;
//
//
//   c = 0;
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = 0; j < n; ++j){
//       m.distance2(i, j);
//       c++;
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Distances 2: %f ms, %d\n", std::chrono::duration<double, std::milli>(finish - start).count(), c);
//
//   c = 0;
//   start = std::chrono::high_resolution_clock::now();
//   for (int i = 0; i < n; ++i){
//     for (int j = 0; j < n; ++j){
//       m.distance1(i, j);
//       c++;
//     }
//   }
//   finish = std::chrono::high_resolution_clock::now();
//   Rprintf("Distances 1: %f ms, %d\n", std::chrono::duration<double, std::milli>(finish - start).count(), c);
//
//
//
//
//   // const int n = x.nrow();
//   // start = std::chrono::high_resolution_clock::now();
//   // for (int i = 0; i < n; ++i){
//   //   ptr_idx_map.emplace(&pt_array[i], i);
//   //   // Rprintf("%p --> %d\n", &pt_array[i], i);
//   // }
//   // finish = std::chrono::high_resolution_clock::now();
//   // Rprintf("Emplacing: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//   //
//   // // ANNpoint& x_i = pt_array[5];
//   // // Rprintf("%p: %d\n", x_i, ptr_idx_map[&x_i]);
//   //
//   // ANNdist dist = 0.0;
//   // start = std::chrono::high_resolution_clock::now();
//   // for (int i = 0; i < n - 1; ++i){
//   //   ANNpoint& x_i_ref = pt_array[i];
//   //   ANNpoint& x_j_ref = pt_array[i+1];
//   //   ANNpoint x_i = pt_array[ptr_idx_map[&x_i_ref]];
//   //   ANNpoint x_j = pt_array[ptr_idx_map[&x_j_ref]];
//   //   dist += (*x_j - *x_i) * (*x_j - *x_i);
//   // }
//   // finish = std::chrono::high_resolution_clock::now();
//   // Rprintf("Accessing: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//   // Rcout << "Dist: " << dist << std::endl;
//   //
//   //
//   // dist = 0.0;
//   // start = std::chrono::high_resolution_clock::now();
//   // for (int i = 0; i < n - 1; ++i){
//   //   ANNpoint x_i = pt_array[i];
//   //   ANNpoint x_j = pt_array[i+1];
//   //   dist += (*x_j - *x_i) * (*x_j - *x_i);
//   // }
//   // finish = std::chrono::high_resolution_clock::now();
//   // Rprintf("Vector: %f ms\n", std::chrono::duration<double, std::milli>(finish - start).count());
//   // Rcout << "Dist: " << dist << std::endl;
//
//   return(x);
// }
//
//
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
//
// /*** R
// x <- cbind(rnorm(40000), rnorm(40000))
// x2 <- timesTwo(x)
// */
