//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.h
//----------------------------------------------------------------------
// Copyright (c) 2017 Matt Piekenbrock. 

// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include <Rcpp.h>
using namespace Rcpp;

class UnionFind
{
  Rcpp::IntegerVector parent;
  Rcpp::IntegerVector rank;
  
  public:
  UnionFind(const int size);
  ~UnionFind();
  void Union(const int x, const int y); 
  const int Find(const int x); 
  
}; // class UnionFind