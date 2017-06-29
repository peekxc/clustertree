#include <Rcpp.h>
using namespace Rcpp;

// Extensions to the core ANN KD tree node capabilities
#include "kd_tree.h"

// --- Splitting node extensions ---
void ANNkd_split::child_nodes(std::vector<ANNkd_node*>& nodes){
  ANNkd_ptr left_child = child[0], right_child = child[1];
  nodes.push_back(left_child);
  nodes.push_back(right_child);
  return;
}
void ANNkd_split::desc_nodes(std::vector<ANNkd_node*>& nodes){
  ANNkd_ptr left_child = child[0], right_child = child[1];

  // Store left child, descend left branch
  nodes.push_back(left_child);
  left_child->desc_nodes(nodes);

  // Store right child, descend right branch
  nodes.push_back(right_child);
  right_child->desc_nodes(nodes);
}
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
// ANNcoord* ANNkd_split::convex_subset(){
//   return(cd_bnds);
// }

// --- Leaf node extensions ---
void ANNkd_leaf::child_nodes(std::vector<ANNkd_node*>& nodes){
  return;
}
void ANNkd_leaf::desc_nodes(std::vector<ANNkd_node*>& nodes){
  return;
}
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
  std::vector<int>* pts_in_node = new std::vector<int>();
  this->held_in_node(*pts_in_node); // get the points held in this node
  double max_child_dist = std::numeric_limits<double>::infinity(), dist;

  // for(std::vector<int>::iterator pt = pts_in_node->begin(); pt != pts_in_node->end(); ++pt){
  //   dist = annDist(this->ro , centroid, *pt);
  //
  // }
  return 0;
}
ANNdist ANNkd_leaf::max_desc_dist(){
  return 0;
}

// ANNcoord* ANNkd_leaf::convex_subset(){
//   ANNcoord* rec = (ANNcoord*) calloc(2, sizeof(ANNcoord));
//
//   ANNorthRect bnd_box(this-> );
//   void annEnclRect(
//       ANNpointArray		pa,				// point array
//       ANNidxArray			pidx,			// point indices
//       int					n,				// number of points
//       int					dim,			// dimension
//       ANNorthRect			&bnds
//   ANNorthRect
//   // this->getStats(this->)
//
//   return(rec);
// }


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

