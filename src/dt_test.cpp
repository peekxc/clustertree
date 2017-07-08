#include <Rcpp.h>
using namespace Rcpp;

#include "dt.h"

// void DualTree::test_cases(List& res){
//
//   // Print the reference tree
//   // rtree->Print((ANNbool) true, std::cout);
//
//   // Initialize KNN structures + bounds map
//   knn_initialize(5);
//
//   // Test all methods work for reference tree (recursively)
//   R_INFO(" ======== Beginning test cases for REFERENCE tree ======== ")
//     test_cases(res, rtree->root, 0, true);
//
//   // Test all methods work for query tree (recursively)
//   // Rcout << " ======== Beginning test cases for QUERY tree ======== " << std::endl;
//   // test_cases(res, qtree->root, 0, false);
// }
//
// int traversed = 0; // How many recursions have there been?
// void DualTree::test_cases(List& res, ANNkd_node* N, int depth, bool ref_tree) {
//   traversed++;
//
//   // Get child nodes
//   std::vector<ANNkd_node*> child_nodes = std::vector<ANNkd_node*>();
//   N->child_nodes(child_nodes);
//
//   // Print the point address to the current node for debugging
//   Rcout << "Node ptr: " << N << std::endl;
//
//   // Create bounds object for testing
//   if (bounds->find(N) == bounds->end()){ bounds->insert(std::pair<ANNkd_node*, Bound>(N, Bound())); }
//   Bound& N_bound = (*bounds)[N];
//
//   // Test bounding box calculation
//   ANNorthRect& bnd_box = convex_subset(N, N_bound, ref_tree);
//   if (N_bound.bnd_box != NULL){
//     Rcout << "Bounding Box: (";
//     for (int i = 0; i < d; ++i) { Rcout << N_bound.bnd_box->lo[i] << ", "; }
//     Rcout << "), (";
//     for (int i = 0; i < d; ++i) { Rcout << N_bound.bnd_box->hi[i] << ", "; }
//     Rcout << ")" << std::endl;
//   }
//
//   // Test retrieving the centroid
//   ANNpoint node_centroid = centroid(N, N_bound, ref_tree);
//   if (N_bound.centroid != NULL){
//     Rcout << "Centroid: ";
//     for (int i = 0; i < d; ++i) { Rcout << N_bound.centroid[i] << ", "; }
//     Rcout << std::endl;
//   }
//
//   // Test getting the max child distance
//   ANNdist n_max_cd = max_child_dist(N, N_bound, ref_tree);
//   Rcout << "Max child dist: " << n_max_cd << " (now " << N_bound.rho << ")" << std::endl;
//
//   // Test getting the max desc. distance
//   ANNdist n_max_dd = max_desc_dist(N, N_bound, ref_tree);
//   Rcout << "Max child dist: " << n_max_dd << " (now " << N_bound.lambda << ")" << std::endl;
//
//   // Recurse
//   // Rcout << "Children: ";
//   for (std::vector<ANNkd_node*>::iterator node = child_nodes.begin(); node != child_nodes.end(); ++node){
//     //Rcout << *node << ", ";
//     test_cases(res, (ANNkd_node*) *node, depth + 1, ref_tree);
//   }
//   // Rcout << std::endl;
//   //
//   //   // Test cases: starting at leaves, and working up the tree, test both static and recursive methods
//   //   Rcout << "Depth: " << depth << ", Traversed: " << traversed << std::endl;
//   //   R_FlushConsole();
//   //
//   //   // Test node id retrieval
//   //   std::vector<int> node_ids = std::vector<int>();
//   //   N->node_ids(node_ids);
//   //   Rcout << "Node ids: ";
//   //   for (int i = 0; i < node_ids.size(); ++i){ Rcout << node_ids.at(i) << ", "; }
//   //   Rcout << std::endl;
//
//   // Test Descendent (recursive) node retrieval
//   // std::vector<int> desc_ids = std::vector<int>();
//   // N->desc_ids(desc_ids);
//   // Rcout << "Desc. node ids: ";
//   // for (int i = 0; i < desc_ids.size(); ++i){
//   //   Rcout << desc_ids.at(i) << ", ";
//   // }
//   // Rcout << std::endl;
//
//
//   // std::copy(r_centroid, r_centroid + d, r_res.begin());
//   // std::copy(q_centroid, q_centroid + d, q_res.begin());
//   // res["r_centroid"] = r_res;
//   // res["q_centroid"] = q_res;
//   //
//   // // Test the various node extensions
//   // std::vector<int> desc_ids = std::vector<int>(), node_ids = std::vector<int>();
//   // std::vector<ANNkd_node*> desc_nodes = std::vector<ANNkd_node*>(), child_nodes = std::vector<ANNkd_node*>();
//   // rtree->root->node_ids(node_ids);
//   // rtree->root->child_nodes(child_nodes);
//   // rtree->root->desc_ids(desc_ids);
//   // rtree->root->desc_nodes(desc_nodes);
//   //
//   // res["node_ids"] = wrap(node_ids); // should be empty
//   // res["desc_ids"] = wrap(desc_ids); // should contain all ids, inorder
//
//
//   // Test bounds checking
//   // res["bounds"] = B(N);
//
//   // Return
//   return;
// }
//
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
//
// /*** R
// timesTwo(42)
// */
