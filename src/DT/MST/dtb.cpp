#include <DT/MST/dtb.h> // Header file

DualTreeBoruvka::DualTreeBoruvka(const bool prune, const int dim): DualTree(prune, dim) {
    N_q_par = N_r_par = NULL;
}

// Now that the tree(s) are made, the knn-specific bound obejcts can be instantiated
void DualTreeBoruvka::setup(ANNkd_tree* kd_treeQ, ANNkd_tree* kd_treeR) {
  (*this).DualTree::setup(kd_treeQ, kd_treeR);  // Call parent class setup to assign trees
}
