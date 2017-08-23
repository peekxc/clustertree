//----------------------------------------------------------------------
//                        Disjoint-set data structure
// File:                        union_find.h
//----------------------------------------------------------------------
// Copyright (c) 2017 Matt Piekenbrock.

// Class definition based off of data-structure described here:
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include <Rcpp.h>
using namespace Rcpp;

/* Rcpp UnionFind - implementation of the Disjoint Set data structure
 * Optimized for amortized constant-time using path compression and union by rank
 */
class UnionFind
{
private:
  IntegerVector parent, rank;

public:
  const unsigned int size;
  UnionFind(const int size);
  ~UnionFind();
  void Union(const int x, const int y);
  const int Find(const int x);

  // Retrieve the full integer vector representing component membership
  IntegerVector getCC();

  // Equality comparison operators to simplify finding distinct components
  bool operator==(UnionFind& other_cc);
  bool operator!=(UnionFind& other_cc);

  // Overload equality/assignment
  void merge(UnionFind &other_cc, LogicalVector& admitted);
}; // class UnionFind