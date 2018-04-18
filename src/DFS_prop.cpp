
#include "RcppHeader.h"
#include "ANN.h"
#include "ANN_util.h"
#include "kd_tree.h"
#include <unordered_map>

// struct parent_info{
//
// };
//
//
// void DFS(ANNkd_split* node){
//   DFS(node->child[ANN_LO]);
// }
// void DFS(ANNkd_leaf* node){
//   DFS(node->child[ANN_LO]);
// }


// NumericMatrix TreeProperties(const NumericMatrix& x) {
//   ANNpointArray pt_array = matrixToANNpointArray(x);
//
//   // Create the kd tree
//   NODE_INFO node_info; // container to store various properties related to tree construction
//   ANNkd_tree* kdTreeR = new ANNkd_tree(pt_array, x.nrow(), x.ncol(), 5, (ANNsplitRule) 5, node_info);
//
//   // Collect the parent information
//   std::unordered<ANNkd_node*, parent_info> info =   std::unordered<ANNkd_node*, parent_info>();
//
//   return(x);
// }

/*** R
DFS(42)
*/
