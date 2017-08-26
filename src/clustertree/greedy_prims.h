#ifndef GREEDY_PRIMS_H
#define GREEDY_PRIMS_H

#include <ANN/ANNx.h>
#include <ANN/kd_tree/kd_tree.h>
#include <Rcpp.h>
using namespace Rcpp;

int findNearest(ANNpoint q, ANNkd_tree* qtree, double eps = 0);

#endif