#ifndef NODE_BNB_H
#define NODE_BNB_H

#include "ANNx.h"

// Branch & Bound information
struct node_bnb {
  int id;
  ANNcoord max_radius, bound;
  ANNpoint centroid;
  // ANNpoint lo, hi;
  node_bnb(void);
};

#endif