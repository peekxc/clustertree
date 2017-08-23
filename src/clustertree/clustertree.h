#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include <DT/structures/union_find.h>
#include <DT/MST/dtb.h>
#include <ANN/ANN_util.h>
#include <utilities.h>
#include <metrics.h>
#include <clustertree/dtb_ct.h>

// Computes the connection radius, i.e. the linkage criterion
double getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type);
void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind);
IntegerVector extractOrder(IntegerMatrix merge);
NumericMatrix primsRSL(const NumericVector r, const NumericVector r_k, const int n, const double alpha, const int type);
List mstToHclust(NumericMatrix mst);
List dtbRSL(const NumericMatrix x, const NumericVector r_k, const double alpha, const int type);

#endif