//----------------------------------------------------------------------
//                        Disjoint-set data structure
// File:                        union_find.cpp
//----------------------------------------------------------------------
// Copyright (c) 2017 Matt Piekenbrock. All Rights Reserved.

// Class definition based off of data-structure described here:
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include "union_find.h"

UnionFind::UnionFind(const int size) : parent(size), rank(size), size(size)
{
  for (int i = 0; i < size; ++i)
  { parent[i] = i, rank[i] = 0; }
}

// Destructor not needed w/o dynamic allocation
UnionFind::~UnionFind() { }

void UnionFind::Union(const int x, const int y)
{
  const int xRoot = Find(x);
  const int yRoot = Find(y);

  // x and y are in the same set
  if (xRoot == yRoot)
   return;
  // Use rank to prioritize shorten the max length trees can grow
  else if (rank[xRoot] < rank[yRoot])
    parent[xRoot] = yRoot;
  else if (rank[xRoot] > rank[yRoot])
    parent[yRoot] = xRoot;
  // Arbitrarily assign parent tree
  else {
    parent[yRoot] = xRoot;
    rank[xRoot] = rank[xRoot] + 1;
  }
}

const int UnionFind::Find(const int x)
{
  // Use path compression to speed up subsequent finds
  if (parent[x] != x)
    parent[x] = Find(parent[x]);
  return(parent[x]);
}

// Return new integer vector representing the connected components
IntegerVector UnionFind::getCC(){
  IntegerVector cc = Rcpp::no_init(size);
  for (unsigned int i = 0; i < size; ++i){ cc[i] = Find(i); }
  return(cc);
}

// Simple method to print the CCs on one line
void UnionFind::printCC(){
  for (int j = 0; j < size; ++j){ Rcout << Find(j) << " "; }
  Rcout << std::endl;
}

// Comparison operators. Note that due to path compression, can't use const references
bool UnionFind::operator==(UnionFind& other_cc){
  bool matches = true;
  for (unsigned int i = 0; matches && i < size; ++i){
    matches = this->Find(i) == other_cc.Find(i);
  }
  return(matches);
};
bool UnionFind::operator!=(UnionFind& other_cc){
  return(!(*this == other_cc));
};
