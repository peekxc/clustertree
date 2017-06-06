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
