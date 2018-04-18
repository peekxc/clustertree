#ifndef HCLUST_UTIL_H
#define HCLUST_UTIL_H

#include "RcppHeader.h"
#include "union_find.h"
#include "utilities.h"

// Recursively visit the merge matrix to extract an hclust sufficient ordering
void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind);
IntegerVector extractOrder(IntegerMatrix merge);
List mstToHclust(const IntegerMatrix& merge, const NumericVector& dist);


#endif