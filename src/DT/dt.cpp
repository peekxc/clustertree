#include <Rcpp.h>
using namespace Rcpp;

// Local header includes
#include "dt.h"
#include <ANN/kd_tree/kd_split.h> // ANNsplitRule

// Constructor initializes the bounds map (if pruning is enabled), and base case map
DualTree::DualTree(const bool prune, const int dim, Metric* m): use_pruning(prune), d(dim) {
  if (use_pruning){ bounds = new std::unordered_map<ANNkd_node*, const Bound& >(); }
  BC_check = new std::map< std::pair<int, int>, candidate_pair>();
  m_dist = m == NULL ? new L_2(dim) : m;
}

// Derivable setup function
void DualTree::setTrees(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR){
  // Check dimensionality, then assign trees if the same
  if (kd_treeR->theDim() != kd_treeQ->theDim() || kd_treeR->theDim() != d){ stop("Dimensionality of the query set does not match the reference set."); }
  // R_PRINTF("Setting up dual trees for %d dimensions\n", d)
  rtree = kd_treeR;
  qtree = kd_treeQ;
}

// Use specialized base cases based on whether the trees are identical
// (Identical implying they were built on the same data set, and with the same
// splitting criteria)
inline ANNdist DualTree::BaseCase(ANNkd_node* N_q_leaf, ANNkd_node* N_r_leaf){
  switch(qtree == rtree){
    case true:
      BaseCaseIdentity(N_q_leaf, N_r_leaf);
      break;
    case false:
      BaseCaseNonIdentity(N_q_leaf, N_r_leaf);
      break;
    }
}

// Generic function to comput ethe distance between a query and reference point
// Note that this function uses incremental distance calculations, stopping the
// calculation and returning Infinity if the distance is greater than both the
// given thresholds. This is used to dramatically reduce the number of floating
// point operations in higher-dimensional space
ANNdist DualTree::computeDistance(const int q_idx, const int r_idx, ANNdist eps1, ANNdist eps2){
  ANNcoord* qq = qtree->pts[q_idx];     // first coord of query point
  ANNcoord* pp = rtree->pts[r_idx];			// first coord of reference point
  ANNdist dist = 0;
  bool valid_dist = true;
  int d_i;
  for (d_i = 0; d_i < d; ++d_i) {
    dist += (*m_dist)(qq, pp, d_i);
    if(dist > eps1 && dist > eps2) { // check both since working with tree from identical data set
      break;
    }
  }
  valid_dist = d_i >= d && (ANN_ALLOW_SELF_MATCH || dist!=0); // ensure is valid distance
  return valid_dist ? dist : ANN_DIST_INF;
}

// Given two leaf nodes and an array of "best" distances so far
// void DualTree::compareLeaves(ANNkd_leaf* N_q_leaf, ANNkd_leaf* N_r_leaf,
//                              ANNdist* q_best_dist, // An array of the best
//                              ANNdist* r_best_dist
//                             ){
//   ANNdist dist;				// distance to data point
//   ANNdist min_dist_q, min_dist_r;			// distance to k-th closest point
//   ANNcoord* pp; // first coord of reference point
//   ANNcoord* qq; // first coord of query point
//   ANNcoord t;
//   int d_i;
//
//   ANNpoint p_q, p_r;
//   int q_idx, r_idx;
//   ANNdist min_dist_q, min_dist_r; // the minimum
//
//   for (int q_i = 0, q_idx = N_q_leaf->bkt[q_i]; q_i < N_q_leaf->n_pts; ++q_i, q_idx = N_q_leaf->bkt[q_i]){
//     min_dist_r = EL[q_idx].weight; // k-th smallest distance so far (query)
//     for (int r_i = 0, r_idx = N_r_leaf->bkt[r_i]; r_i < N_r_leaf->n_pts; ++r_i, r_idx = N_r_leaf->bkt[r_i]){
//
//
//
//       // Compute Base case, saving knn ids and distances along the way
//       if (!hasBeenChecked(q_idx, r_idx)) { ANN_PTS(2) // Has this pair been considered before?
//         min_dist_r = EL[r_idx].weight; // k-th smallest distance so far (reference)
//
//         qq = qtree->pts[q_idx];     // first coord of query point
//         pp = rtree->pts[r_idx];			// first coord of reference point
//         dist = 0;
//
//         // Incrementally compute distance. If at any dimension the distance exceeds the kth smallest
//         // distance so far, continue on to the next reference point
//         for(d_i = 0; d_i < d; d_i++) {
//           ANN_COORD(1)				// one more coordinate hit
//           ANN_FLOP(4)					// increment floating ops
//
//           t = *(qq++) - *(pp++);		// compute length and adv coordinate
//           // exceeds dist to k-th smallest?
//           dist = ANN_SUM(dist, ANN_POW(t));
//           if(dist > min_dist_q && dist > min_dist_r) { // check both since working with tree from identical data set
//             break;
//           }
//         }
//       }
//     } // end ref pts
//   } // end query pts
//
//
// }


void DualTree::printNode(ANNkd_leaf* N, int level){
  Rcout << "    ";
  for (int i = 0; i < level; i++)		// print indentation
    Rcout << "..";

  if (N == KD_TRIVIAL) {			// canonical trivial leaf node
    Rcout << "Leaf (trivial)\n";
  }
  else{
    Rcout << "Leaf n=" << N->n_pts << " <";
    for (int j = 0; j < N->n_pts; j++) {
      Rcout << N->bkt[j];
      if (j < N->n_pts-1) Rcout << ",";
    }
  #ifdef NDEBUG
      Rcout << "> (" << node_labels.at(N) << ")\n";
  #else
      Rcout << ">\n";
  #endif
  }
}

void DualTree::printNode(ANNkd_split* N, int level){
  IS_LEAF(N->child[ANN_HI]) ? printNode(AS_LEAF(N->child[ANN_HI]), level + 1) : printNode(AS_SPLIT(N->child[ANN_HI]), level + 1);
  Rcout << "    ";
  for (int i = 0; i < level; i++)		// print indentation
    Rcout << "..";
  Rcout << "Split cd=" << N->cut_dim << " cv=" << N->cut_val;
  Rcout << " lbnd=" << N->cd_bnds[ANN_LO];
  Rcout << " hbnd=" << N->cd_bnds[ANN_HI];
#ifdef NDEBUG
  Rcout << "> (" << node_labels.at(N) << ")\n";
#else
  Rcout << ">\n";
#endif
  IS_LEAF(N->child[ANN_LO]) ? printNode(AS_LEAF(N->child[ANN_LO]), level + 1) : printNode(AS_SPLIT(N->child[ANN_LO]), level + 1);
}

// TODO: Export better to R, figure out how to extract R output buffer instead of cout for windows users
void DualTree::PrintTree(ANNbool with_pts, bool ref_tree){
  ANNkd_tree* ctree = ref_tree ? rtree : qtree;
  Rcout << "ANN Version " << ANNversion << "\n";
  if (with_pts) {						// print point coordinates
    Rcout << "    Points:\n";
    for (int i = 0; i < ctree->n_pts; i++) {
      Rcout << "\t" << i << ": ";
      ANNpoint pt = ctree->pts[i];
      for (int j = 0; j < d; j++) {
        Rcout << pt[j];
        if (j < d-1) Rcout << " ";
      }
      Rcout << "\n";
    }
  }
  if (ctree->root == NULL)					// empty tree?
    Rcout << "    Null tree.\n";
  else {
    IS_LEAF(ctree->root) ? printNode(AS_LEAF(ctree->root), 0) : printNode(AS_SPLIT(ctree->root), 0);
  }
}

ANNkd_tree* DualTree::ConstructTree(ANNpointArray x, const int nrow, const int ncol, const int bkt_sz, ANNsplitRule split_rule){
  // Create kd tree, either a dual tree version with bounds or a regular ANN kd tree
  ANNkd_tree* kdTree;
  if (use_pruning){
    kdTree = new ANNkd_tree_dt(x, nrow, ncol, bkt_sz, ANN_KD_SUGGEST, bounds); // use bound base class
    R_INFO("Bounds computed at tree construction: " << bounds->size() << "\n")
    #ifdef NDEBUG
      node_labels = *new std::unordered_map<ANNkd_node*, char>();
      char letter='A';
      for (std::unordered_map<ANNkd_node*, const Bound& >::iterator it = bounds->begin(); it != bounds->end(); ++it){
        const Bound& c_bnd = it->second;
        node_labels.insert(std::make_pair(it->first, letter));
        Rcout << "Node: " << letter << ", Centroid: (";
        for (int i = 0; i < ncol; ++i){ Rcout << c_bnd.centroid[i] << ", "; }
        Rcout << "), Lambda: " << c_bnd.lambda << ", Rho: " << c_bnd.rho;
        // Rcout << " Bnd Box: ";
        //for (int i = 0; i < ncol; ++i){ Rcout <<  << ", "; }
        letter++;
      }
    #endif
  } else {
    kdTree = new ANNkd_tree(x, nrow, ncol, bkt_sz, ANN_KD_SUGGEST);
  }
  return kdTree;
}


// The minimum distance between any two points descendent of N_i and N_j
// ANNdist DualTree::min_dist(ANNkd_node* N_i, ANNkd_node* N_j){
//   if (N_i == N_j) return (ANNdist) 0.0; // Trivial case--nodes are the same
//   //Rcout << "Computing the minimum distance" << std::endl;
//   std::vector<int>* i_desc_ids = new std::vector<int>();
//   std::vector<int>* j_desc_ids = new std::vector<int>();
//   N_i->desc_ids(*i_desc_ids);
//   N_j->desc_ids(*j_desc_ids);
//
//   ANNdist min_dist = ANN_DIST_INF, dist;
//   for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
//     for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
//       ANNpoint p_i = this->rtree->pts[*i];
//       ANNpoint p_j = this->qtree->pts[*j];
//       dist = annDist(d, p_i, p_j);
//       if (dist < min_dist){
//         min_dist = dist;
//       }
//     }
//   }
//   return min_dist;
// }
//
// // The maximum distance between any two points descendent of N_i and N_j
// ANNdist DualTree::max_dist(ANNkd_node* N_i, ANNkd_node* N_j){
//   std::vector<int>* i_desc_ids = new std::vector<int>();
//   std::vector<int>* j_desc_ids = new std::vector<int>();
//   N_i->desc_ids(*i_desc_ids);
//   N_j->desc_ids(*j_desc_ids);
//
//   ANNdist max_dist = 0, dist;
//   for (std::vector<int>::iterator i = i_desc_ids->begin(); i != i_desc_ids->end(); ++i){
//     for (std::vector<int>::iterator j = j_desc_ids->begin(); j != j_desc_ids->end(); ++j){
//       ANNpoint p_i = this->rtree->pts[*i];
//       ANNpoint p_j = this->qtree->pts[*j];
//       if ((dist = annDist(d, p_i, p_j)) > max_dist){
//         max_dist = dist;
//       }
//     }
//   }
//   return max_dist;
// }
//
// // Get the maximum distance between the centroid of the convex subset (the box of the current node)
// // and the points within the children of the current node
// ANNdist DualTree::max_child_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
//
//   // If the maximum child distance has already been computed, return it
//   if (ni_bnd.rho != -1.0) return ni_bnd.rho;
//
//   // Else look at the points its holding
//   std::vector<int> point_ids = std::vector<int>();
//   N_i->node_ids(point_ids);
//
//   // Centroid computes (and stores) the bounding box, so access should be fine past here
//   ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;
//
//   // If it's a splitting node it contains no points. There are two options:
//   // 1) Search all of the descendents, computing the true maximum child distance
//   // 2) Approximate a looser upper bound using the distance between the precomputed centroid
//   // and the bounding rectangle orthogonal to the cutting dimension
//   if (point_ids.size() == 0){
//     ni_bnd.rho = annDist(d, ni_centroid, ni_bnd.bnd_box->hi); // <-- Looser upper bound
//     // Rcout << "== Computed child distance: " << ni_bnd.rho << std::endl;
//     // ni_bnd.rho = max_desc_dist(N_i, ni_bnd, ref_tree); // <-- Tighter but slower to compute bound
//     return ni_bnd.rho;
//   } else {
//     // If it's a leaf node and it contains no points, find the one farthest from the centroid
//     ANNdist max_dist = 0;
//     for (std::vector<int>::iterator pt_idx = point_ids.begin(); pt_idx != point_ids.end(); ++pt_idx){
//       ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
//       ANNdist dist = annDist(d, cpt, ni_centroid);
//       // Rcout << "Testing child distance: " << dist << std::endl;
//       if (dist > max_dist){ max_dist = dist; }
//     }
//     ni_bnd.rho = max_dist; // store the bound
//     return max_dist;
//   }
// }
//
// // Get the maximum distance between the centroid of the convex subset (the box of the current node)
// // and all of the descendent points within a given node
// ANNdist DualTree::max_desc_dist(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
//
//   // If the maximum desc. distance has already been computed, return it
//   if (ni_bnd.lambda != -1.0) return ni_bnd.lambda;
//
//   // Get the centroid
//   ANNpoint ni_centroid = (ni_bnd.centroid == NULL) ? (ANNpoint) this->centroid(N_i, ni_bnd, ref_tree) : ni_bnd.centroid;
//
//   // Get all descendent point indices
//   std::vector<int> desc_ids = std::vector<int>();
//   N_i->desc_ids(desc_ids);
//   ANNdist max_dist = 0;
//   bool check_point_knn = false;
//   for (std::vector<int>::iterator pt_idx = desc_ids.begin(); pt_idx != desc_ids.end(); ++pt_idx){
//     ANNpoint cpt = ref_tree ? rtree->pts[*pt_idx] : qtree->pts[*pt_idx];
//     ANNdist dist = annDist(d, cpt, ni_centroid);
//     if (dist > max_dist){ max_dist = dist; }
//     // if (knn->at(*pt_idx)->max_key() != ANN_DIST_INF){
//     //   check_point_knn = true;
//     //   Rcout << "*** " << N_i << " has finite knn points filled!" << std::endl;
//     // }
//   }
//   ni_bnd.lambda = max_dist; // Store the maximum distance by reference for later
//   return max_dist;
// }
//
//
// #include "kd_util.h"
// #include "ANN/ANNx.h"
// #include "ANN/ANNperf.h"
//
// ANNorthRect& DualTree::convex_subset(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
//   if (ni_bnd.bnd_box != NULL) return *ni_bnd.bnd_box;
//   // Get child ids
//   std::vector<int> desc_ids = std::vector<int>();
//   N_i->desc_ids(desc_ids);
//   ANNidxArray pidx = (ANNidxArray) &desc_ids[0];
//
//   // Compute the bounding rectangle
//   ANNorthRect* bnd_box = new ANNorthRect(d); // bounds are dynamically allocated
//   annEnclRect((ANNpointArray) ref_tree ? rtree->pts : qtree->pts,
//               (ANNidxArray) pidx,
//               (int) desc_ids.size(),
//               (int) d,
//               *bnd_box);
//   ni_bnd.bnd_box = bnd_box; // Store bounding rectangle
//   return *bnd_box;
// }
//
// ANNpoint DualTree::centroid(ANNkd_node* N_i, Bound& ni_bnd, bool ref_tree){
//   if (ni_bnd.centroid != NULL) return(ni_bnd.centroid);
//   ANNorthRect& bnds = ni_bnd.bnd_box == NULL ? convex_subset(N_i, ni_bnd, ref_tree) : *ni_bnd.bnd_box;
//   ANNpoint n_centroid = new ANNcoord[d];
//   for (int i = 0; i < d; ++i){
//     // Rcout << "Bounds (lo, hi): " << bnds.lo[i] << ", " << bnds.hi[i] << std::endl;
//     n_centroid[i] = ANNcoord((bnds.lo[i] + bnds.hi[i]) / 2.0);
//     // Rcout << "Centroid: " << centroid[i] << std::endl;
//   }
//   ni_bnd.centroid = n_centroid;
//   return(n_centroid);
// }
