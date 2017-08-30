#include <Rcpp.h>
using namespace Rcpp;

#include "utilities.h"


int global_cid = 1;

// Framework for the Optimal Selection of Clusters
struct FOSC{

  void startDFS(){
    const int n = merge.nrow() + 1;
    int k = n - 2;
    IntegerVector split = merge(k, _);
  }

};


// Simple struct to keep track of branch information through the recursive calls
struct branch_info {
  std::list<int> ids; // list to track the cluster ids in branch (acculumative)
  IntegerVector points; // continuously update vector to track point ids
  branch_info(){
    ids = std::list<int>();
    points = IntegerVector();
  }
  void add2(int left, int right, const List& pt_map){
    ids.push_back(left);
    ids.push_back(right);
    IntegerVector left_ids = as<IntegerVector>(pt_map[left]);
    IntegerVector right_ids = as<IntegerVector>(pt_map[right]);
    std::copy(left_ids.begin(), left_ids.end(), points.end());
    std::copy(right_ids.begin(), right_ids.end(), points.end());
  }
};


// Given the map of branch ids to point ids, and a list a branch ids to extract,
// extract all of the point indices into a single integer vector
IntegerVector extract_ids(const List& pt_map, std::list<int>& cids){

  // Number of clusters
  const size_t n_cids = cids.size();

  // Figure out the length of the output vector
  size_t total_length = 0;
  for (std::list<int>::iterator cid = cids.begin(); cid != cids.end(); ++cid){
    total_length += as<IntegerVector>(pt_map[*cid]).size();
  }

  // Copy the points
  size_t index = 0;
  IntegerVector out = no_init(total_length);
  for (std::list<int>::iterator cid = cids.begin(); cid != cids.end(); ++cid){
    IntegerVector pt_ids = as<IntegerVector>(pt_map[*cid]);
    std::copy(pt_ids.begin(), pt_ids.end(), out.begin() + index);
    index += pt_ids.size();// update index
  }

  return out;
}

// Merge 3 integer vectors
IntegerVector merge3(const IntegerVector& b1, const IntegerVector& b2, const IntegerVector& b3){
  IntegerVector out = no_init(b1.size() + b2.size() + b3.size());
  std::copy(b1.begin(), b1.end(), out.begin());
  std::copy(b2.begin(), b2.end(), out.begin() + b1.size());
  std::copy(b3.begin(), b3.end(), out.begin() + b1.size() + b2.size());
  return out;
}

bool measureScore(const IntegerVector& l, const IntegerVector& p, const IntegerVector& r, Rcpp::Function f){

  // Disjoint cluster scores
  double l_score = f(l), r_score = f(r), v_score = f(p);

  // Score if they were all combined
  IntegerVector all_ids = merge3(l, p, r);
  double m_score = f(all_ids);

  // Return whether the merged score is greater than the child scores
  return m_score > l_score + r_score + v_score;
}


// F must be a scoring function that accepts an integer vector of ids
branch_info hcl_DFS(const int k, const IntegerMatrix& merge, const List& pt_map, Rcpp::Function f, std::list<int>& cl,
                    std::unordered_map<int>& score_map, bool propagate = true){

  // Base case
  if (all(split < 0).is_true()){
    int left = split[0], right = split[1];
    bool parentBetter = measureScore(pt_map[patch::to_string(left)],  // Point ids in left leaf
                                     pt_map[patch::to_string(k),      // Point ids in splitting branch
                                     pt_map[patch::to_string(right)], // Point ids in right leaf
                                     f);                              // Scoring function
    branch_info current_branch = branch_info();
    current_branch.add2(left, right, pt_map); // left and right point indices
    if (parentBetter){
      cl.push_back(k); // Keep the parent
    } else {
      cl.push_back(left); // Keep the left leaf
      cl.push_back(right); // Keep the right leaf
    }
    return current_branch;
  }
  // If one is a split and one is a leaf
  else if (split[0] < 0 || split[1] < 0) {
    const bool isLeaf = split[0] < 0;
    int leaf = isLeaf ? split[0] : split[1];
    int split_node = isLeaf ? split[1] : split[0];

    // Get the list of the best clusters for the split node
    branch_info current_cl = hcl_DFS(split_node, merge, p_map, f);

    // Gather the (extracted) point indices
    IntegerVector l_ids = pt_map[patch::to_string(leaf)];
    IntegerVector s_ids = current_cl.points;
    IntegerVector v_ids = pt_map[patch::to_string(k)];
    IntegerVector m
    // Compute the score. If propagating the score through the hierarchy, save the final
    // max score under the current branch node in the score map. Otherwise, 'refresh' the
    // score by computing the raw score for all of the points at every level.
    bool parentBetter = false;
    double total_score = 0;
    if (propagate && score_map.find(split_node) != score_map.end()){
      double s_score = score_map.at(split_node);
      double l_score = f(l_ids);
      double v_score = f(v_ids);
      double
      parentBetter = s_score + l_score + v_score;

      IntegerVector all_ids = merge3(pt_map[patch::to_string(leaf)],        // Point ids in the leaf
                                     pt_map[patch::to_string(k)],            // Point ids in splitting branch
                                           current_cl.points);
      double m_score = f(all_ids);
    } else {
      parentBetter = measureScore(pt_map[patch::to_string(leaf)],        // Point ids in the leaf
                                  pt_map[patch::to_string(k)],            // Point ids in splitting branch
                                  current_cl.points,                    // Point ids in the split node
                                  f);                                    // Scoring function

    }

    // Make the choice to keep the parent cluster id or the children
    if (parentBetter){
      cl.push_back(k); // Keep the parent (current) branch
      cl.remove(leaf); // Remove leaf from list of candidates
      for(std::list<int>::iterator cid = current_cl.ids.begin(); cid != current_cl.ids.end(); ++cid){
        cl.remove(*cid); // Remove all children ids from the current best clustering list
      }
      score_map.at(k) =
      Rcout << "keeping parent branch." << std::endl;
    } else {
      // Do nothing, children already in the best list

      Rcout << "keeping upper branches." << std::endl;
    }


  }
  // If both were split branches
  else {
    int left = split[0], right = split[1];
    IntegerVector left_pts = hcl_DFS(left, merge, p_map, f, score);
    IntegerVector right_pts = hcl_DFS(right, merge, p_map, f, score);
    IntegerVector parent_pts = pt_map[patch::to_string(k)] -1; // convert to 0-based
    IntegerVector all_ids = merge3(left_pts, right_pts, parent_ids);

    // Calculate cluster scores
    double left_score = f(left_pts);
    double right_score = f(right_pts);
    double virtual_score = f(parent_pts);
    double merge_score = f(all_ids);

    if (merge_score > left_score + right_score + virtual_score){
      cl[all_ids] = ++global_cid;
    } else{
      cl[parent_pts] = 0;
    }

  }
}

// Framework for Optimal Selection of Clusters (FOSC)
// Traverses a cluster tree hierarchy to compute a flat solution, maximizing the:
// - Unsupervised soln: the 'most stable' clusters following the give linkage criterion
// - SS soln w/ instance level Constraints: constraint-based w/ unsupervised tiebreaker
// - SS soln w/ mixed objective function: maximizes J = α JU + (1 − α) JSS (normalized to the unit interval)
// Each recursive call returns the stability and constraint scores of the current branch before normalization.
// The normalization that occurs w.r.t total stability is stored in the "score" attribute of each branch.
// [[Rcpp::export]]
NumericVector fosc(const List& st, std::list<int>& sc) // instance-level constraints
{
  // Base case: at a leaf
  if (){
    List cl = cl_tree[cid];
    sc.push_back(stoi(cid)); // assume the leaf will be a salient cluster until proven otherwise
    return(NumericVector::create((double) cl["stability"],  // Leaf total stability == regular stability score
                                 (double) alpha < 1 ? cl["vscore"] : 0));
  } else {
    // Non-base case: at a merge of clusters, determine which to keep
    List cl = cl_tree[cid];

    // Get child stability/constraint scores
    NumericVector scores, stability_scores = NumericVector(), constraint_scores = NumericVector();
    IntegerVector child_ids = cl_hierarchy[cid];
    for (int i = 0, clen = child_ids.length(); i < clen; ++i){
      int child_id = child_ids.at(i);
      scores = fosc(cl_tree, patch::to_string(child_id), sc, cl_hierarchy, alpha, n_constraints, constraints);
      stability_scores.push_back(scores.at(0)); // stability score for child
      constraint_scores.push_back(scores.at(1)); // constraint score for child
    }

    // Compare and update stability & constraint scores
    double split_stability = (double) sum(stability_scores);  // stability if the branches were disjoint
    double merge_stability = (double) cl["stability"]; // stability if the branches were merged
    double split_constraint_score = (double) sum(constraint_scores) + (double) computeVirtualNode(cl["contains"], constraints)/2*n_constraints;
    double merge_constraint_score = (double) cl["vscore"];

    // Compute total scores
    double split_score = alpha * (split_stability/max_stability) + (1 - alpha) * split_constraint_score;
    double merge_score = alpha * (merge_stability/max_stability) + (1 - alpha) * merge_constraint_score;
    bool merge_children = merge_score > split_score; // Whether to merge the child branches or keep them as disjoint clusters
    cl["score"] = merge_children ? merge_score : split_score;
    cl["vscore"] = merge_children ? merge_constraint_score : split_constraint_score;

    // Prune children and add parent (cid) if need be
    if (merge_children && cid != "0") {
      IntegerVector children = all_children(cl_hierarchy, stoi(cid)); // use all_children to prune subtrees
      for (int i = 0, clen = children.length(); i < clen; ++i){ sc.remove(children.at(i)); }
      sc.push_back(stoi(cid));
    }

    // Save scores for traversal up and for later
    cl_tree[cid] = cl;

    // Return this sub trees score
    return(NumericVector::create((double) cl["score"], alpha < 1 ? (double) cl["vscore"] : 0));
  }
}

/*** R
sl <- clustertree::clustertree(iris[, 1:4], k = 5L)
sl_si <- clustertree:::simplified_hclust(sl$hc, 2)

*/
