#include <Rcpp.h>
using namespace Rcpp;

#include "dtb.h"

DualTreeBoruvka::DualTreeBoruvka(const bool prune, const int dim, const int n)
  : DualTree(prune, dim), CC(UnionFind(n)) {
  N_q_par = N_r_par = NULL;
  if (prune){
    R_INFO("Allocating KNN bound map\n")
    bnd_knn = new std::unordered_map<ANNkd_node*, BoundKNN& >();
  }
}

void DualTreeBoruvka::setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR) {
  (*this).DualTree::setup(kd_treeQ, kd_treeR);  // Call parent class setup to assign trees
  if (use_pruning){
    for (std::unordered_map<ANNkd_node*, const Bound& >::iterator bnd = bounds->begin(); bnd != bounds->end(); ++bnd){
      bnd_knn->insert(std::pair<ANNkd_node*, BoundKNN& >(bnd->first, (BoundKNN&) *new BoundKNN()));
    }
  }
}

inline ANNdist DualTreeBoruvka::BaseCase(ANNpoint p_q, ANNpoint p_r, const int q_idx, const int r_idx, ANNkd_node* N_q) {
  if (p_q == p_r) return 0;
  ANNdist dist = annDist(d, p_r, p_q);
  if (CC.Find(q_idx) != CC.Find(r_idx) && dist < D[CC.Find(q_idx)]){
    D[CC.Find(q_idx)] = dist;
    N[CC.Find(q_idx)] = double_edge((int) q_idx, (int) r_idx, (double) dist);
  }
  return dist;
}

// TODO
inline ANNdist DualTreeBoruvka::Score(ANNkd_node* N_q, ANNkd_node* N_r) {
  if (N_q == N_r) return 0;
  ANNdist min_dist_qr = min_dist(N_q, N_r); // minimum distance between two bounding rectangles
  ANNdist best_bound = B(N_q); // "best" lower bound
  if (min_dist_qr < best_bound){
    ANNkd_split* nq_spl = AS_SPLIT(N_q);
    ANNkd_split* nr_spl = AS_SPLIT(N_q);
    return ANN_DIST_INF; // Prune this branch
  }
  return min_dist_qr; // Recurse into this branch
}

void DualTreeBoruvka::DTB(NumericMatrix& x){


}

ANNdist DualTreeBoruvka::B(ANNkd_node* N_q){
  // If B has been computed before, return it, otherwise look at what needs computed
  BoundKNN& nq_knn_bnd = (BoundKNN&) (*bnd_knn)[N_q];
  if (nq_knn_bnd.B != ANN_DIST_INF){ return(nq_knn_bnd.B); }
  else { // Check to see if we can get better bound

    // Retrieve regular bound object
    const Bound& nq_bound = (const Bound&) (*bounds)[N_q];

    // There are 4 bounds to compute, although some apply only to leaf or split nodes
    ANNdist bound;
    if (IS_SPLIT(N_q)){
      ANNkd_split* nq_split = AS_SPLIT(N_q); // get split node
      ANNdist lchild_B = (*bnd_knn)[nq_split->child[ANN_LO]].B, rchild_B = (*bnd_knn)[nq_split->child[ANN_HI]].B;
      bound = std::min(
        std::min(lchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_LO]].lambda),
                 rchild_B + 2 * (nq_bound.lambda - (*bounds)[nq_split->child[ANN_HI]].lambda)),
                 N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B);
    } else {
      ANNkd_leaf* nq_leaf = AS_LEAF(N_q); // get leaf node
      ANNdist max_knn = nq_knn_bnd.maxKNN(nq_leaf->n_pts);
      bound = std::min(std::min(max_knn, nq_knn_bnd.min_knn + nq_bound.rho + nq_bound.lambda),
                       N_q_par == NULL ? ANN_DIST_INF : (*bnd_knn)[N_q_par].B);
    }

    // Final bound (lower)
    nq_knn_bnd.B = bound;
    return nq_knn_bnd.B;
  }
}

/*** R
timesTwo(42)
*/
