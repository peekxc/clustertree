#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include <DT/structures/union_find.h>
#include <utilities.h>

// Computes the connection radius, i.e. the linkage criterion
inline double getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type);
void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind);
IntegerVector extractOrder(IntegerMatrix merge);
NumericMatrix primsRSL(const NumericVector r, const NumericVector r_k, const int n, const double alpha, const int type);
List mstToHclust(NumericMatrix mst, const int n);

#endif