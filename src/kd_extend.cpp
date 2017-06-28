#include <Rcpp.h>
using namespace Rcpp;

// Extensions to the core ANN KD tree node capabilities
#include "kd_tree.h"

// --- Splitting node extensions ---
void ANNkd_split::child_ids(std::vector<int>& ids){
  ANNkd_ptr left_child = child[0], right_child = child[1];
  left_child->child_ids(ids);
  right_child->child_ids(ids);
}

void ANNkd_split::held_in_node(std::vector<int>& ids){ return; }

ANNdist ANNkd_split::max_child_dist(){
  return 0;
}
ANNdist ANNkd_split::max_desc_dist(){
  return 0;
}

// --- Leaf node extensions ---
void ANNkd_leaf::child_ids(std::vector<int>& ids){
  for (int i = 0; i < this->n_pts; ++i){
    ids.push_back(this->bkt[i]);
  }
}
void ANNkd_leaf::held_in_node(std::vector<int>& ids){
  for (int i = 0; i < this->n_pts; ++i){
    ids.push_back(this->bkt[i]);
  }
}

ANNdist ANNkd_leaf::max_child_dist(){
  return 0;
}
ANNdist ANNkd_leaf::max_desc_dist(){
  return 0;
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
// Unimplemented
// void ANNbd_shrink::ann_dt_search(ANNdist box_dist){
//
// }

