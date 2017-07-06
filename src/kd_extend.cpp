#include <Rcpp.h>
using namespace Rcpp;

// Extensions to the core ANN KD tree node capabilities
#include "kd_tree.h"

// ---------------------------------------------------------
// --------------- Splitting node extensions ---------------
// ---------------------------------------------------------
void ANNkd_split::child_nodes(std::vector<ANNkd_node*>& nodes){
  ANNkd_ptr left_child = child[0], right_child = child[1];
  nodes.push_back(left_child);
  nodes.push_back(right_child);
  return;
}

// Internal splitting nodes contain 0 pts in the ANN library
void ANNkd_split::node_ids(std::vector<int>& ids){
  return;
}

// Recursively store the descendent node (ANNkd_ptr's). Goes left (points
// with negative distances relative to the cutting hpyerplane) then right.
void ANNkd_split::desc_nodes(std::vector<ANNkd_node*>& nodes){
  ANNkd_ptr left_child = child[0], right_child = child[1];

  // Store left child, descend left branch
  nodes.push_back(left_child);
  left_child->desc_nodes(nodes);

  // Store right child, descend right branch
  nodes.push_back(right_child);
  right_child->desc_nodes(nodes);
}


void ANNkd_split::desc_ids(std::vector<int>& ids){
  ANNkd_ptr left_child = child[0], right_child = child[1];
  // Splitting nodes don't contain points, so simply recurse
  left_child->desc_ids(ids);
  right_child->desc_ids(ids);
}

// ---------------------------------------------------------
// --------------- Leaf node extensions --------------------
// ---------------------------------------------------------

// Leaves don't contain children, so return
void ANNkd_leaf::child_nodes(std::vector<ANNkd_node*>& nodes){ return; }

// Push back all of the point ids in the bucket
void ANNkd_leaf::node_ids(std::vector<int>& ids){
  ids.insert(ids.end(), &bkt[0], &bkt[n_pts]); // 1 past the end should ok
}

// Leaves don't contain descendent nodes, so return
void ANNkd_leaf::desc_nodes(std::vector<ANNkd_node*>& nodes){ return; }

// Push back all of the point ids in the bucket
void ANNkd_leaf::desc_ids(std::vector<int>& ids){
  ids.insert(ids.end(), &bkt[0], &bkt[n_pts]); // 1 past the end should ok
}

// --- BD tree extensions ---
// #include "bd_tree.h"
// void ANNbd_shrink::child_ids(std::vector<int>& ids){
//   return;
// }
// void ANNbd_shrink::held_in_node(std::vector<int>& ids){
//   return;
// }
//
// ANNdist ANNbd_shrink::max_child_dist(){
//   return 0;
// }
// ANNdist ANNbd_shrink::max_desc_dist(){
//   return 0;
// }


// void ANNbd_shrink::held_in_node(std::vector<int>& ids){
//
// }

