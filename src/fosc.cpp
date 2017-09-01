#include <Rcpp.h>
using namespace Rcpp;

#include "utilities.h"


int global_cid = 1;

// struct node {
//   int id;
// };
//
// struct split_node : node {
//
// };
//
// struct leaf_node : node {
//   IntegerVector idx;
// };

#include <Rcpp.h>
using namespace Rcpp;

// //' @export SimpleHclust
// class SimpleHclust {
// public:
//   SimpleHclust(int x): x_(x) {}
//   int getValue() { return x_; }
//   void addValue(int y) { x_ += y; }
//   void merge(const SimpleHclust& rhs) { x_ += rhs.x_; }
// private:
//   int x_;
// };
//
//
// RCPP_EXPOSED_CLASS(SimpleHclust)
//   RCPP_MODULE(mod_simplehclust) {
//
//     class_<SimpleHclust>("SimpleHclust")
//
//     .constructor<int>("sets initial value")
//
//     .method("getValue", &SimpleHclust::getValue, "Returns the value")
//     .method("addValue", &SimpleHclust::addValue, "Adds a value")
//     .method("merge", &SimpleHclust::merge, "Merges another SimpleHclust into this object")
//     ;
//   }


// Framework for the Optimal Selection of Clusters
// struct FOSC{
//   split_node*
//   void startDFS(){
//     const int n = merge.nrow() + 1;
//     int k = n - 2;
//     IntegerVector split = merge(k, _);
//   }
//
//
//
// };
//
void printVec(const IntegerVector& vec){
  for (IntegerVector::const_iterator it = vec.begin(); it != vec.end(); ++it){
    Rcout << *it << ", ";
  }
  Rcout << std::endl;
}
// Simple struct to keep track of branch information through the recursive calls
struct branch_info {
  std::list<int> ids; // list to track the cluster ids in branch (acculumative)
  IntegerVector points; // continuously update vector to track point ids
  branch_info(){
    ids = std::list<int>();
    points = IntegerVector();
  }

  // Add a cluster, and it's associated point ids, to the current branch information
  void add(int cid, const List& pt_map){
    ids.push_back(cid);
    IntegerVector point_ids = as<IntegerVector>(pt_map[patch::to_string(cid)]);
    for (IntegerVector::iterator pid = point_ids.begin(); pid != point_ids.end(); ++pid){
      points.push_back(*pid);
    }
    // printVec(point_ids);
    //std::copy(point_ids.begin(), point_ids.end(), points.begin() + points.size());
  }

  // Merge a branch cluster and point indices into this one
  void merge(branch_info& other_branch){
    ids.insert(ids.end(), other_branch.ids.begin(), other_branch.ids.end());
    std::copy(other_branch.points.begin(), other_branch.points.end(), points.end());
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
  if (b1.size() > 0) std::copy(b1.begin(), b1.end(), out.begin());
  if (b2.size() > 0) std::copy(b2.begin(), b2.end(), out.begin() + b1.size());
  if (b3.size() > 0) std::copy(b3.begin(), b3.end(), out.begin() + b1.size() + b2.size());
  return out;
}

bool measureScore(const IntegerVector& l, const IntegerVector& p, const IntegerVector& r, Rcpp::Function f){

  // Disjoint cluster scores
  double l_score = as<double>(f(l)), r_score = as<double>(f(r)), v_score = as<double>(f(p));

  // Score if they were all combined
  IntegerVector all_ids = merge3(l, p, r);
  double m_score = as<double>(f(all_ids));

  // Return whether the merged score is greater than the child scores
  return m_score > l_score + r_score + v_score;
}


// F must be a scoring function that accepts an integer vector of ids
branch_info& hcl_DFS(const int k, const IntegerMatrix& merge, const List& pt_map, Rcpp::Function f, std::list<int>& cl,
                    std::unordered_map<int, double>& score_map, bool propagate = true){

  // The current points
  IntegerVector split = merge(k, _);
  Rcout << "split: " << split[0] << ", " << split[1] << std::endl;

  // Gather the (extracted) point indices
  int left = split[0], right = split[1];
  IntegerVector l_ids, r_ids;

  // Two leaf case
  if (all(split < 0).is_true()){
    l_ids = pt_map[patch::to_string(left)];
    r_ids = pt_map[patch::to_string(right)];
  } else if (left < 0 || right < 0) {
    const bool isLeaf = left < 0;
    int leaf = isLeaf ? left : right;
    int split_node = isLeaf ? right : left;
    l_ids = pt_map[patch::to_string(left)];
    branch_info current_cl = hcl_DFS(split_node-1, merge, pt_map, f, cl, score_map, propagate);
  } else {

  }




  IntegerVector v_ids = pt_map[patch::to_string(k+1)];
  IntegerVector m_ids = merge3(l_ids, v_ids, r_ids);

  // Get the list of the best clusters for both split nodes
  branch_info cl_l = hcl_DFS(left-1, merge, pt_map, f, cl, score_map, propagate);
  branch_info cl_r = hcl_DFS(right-1, merge, pt_map, f, cl, score_map, propagate);


  // Create a new set of branch information
  branch_info current_branch = branch_info();
  current_branch.add(left, pt_map); // left and right point indices
  current_branch.add(right, pt_map); // left and right point indices
  if (parentBetter){
    cl.push_back(k+1); // Keep the parent
  } else {
    cl.push_back(left); // Keep the left leaf (already 1-based)
    cl.push_back(right); // Keep the right leaf (already 1-based)
  }
  return current_branch;



  // Gather the (extracted) point indices
  IntegerVector l_ids = cl_l.points;
  IntegerVector s_ids = cl_r.points;

  // Calculate cluster scores
  double l_score = as<double>(f(l_ids)); // Left branch score
  double r_score = as<double>(f(r_ids)); // Right branch score
  double v_score = as<double>(f(v_ids)); // Virtual score (poitns inside the parent branch)
  double m_score = as<double>(f(m_ids)); // Merged score (all points)

  // Determine whether parent is better or not; compute final score
  bool parentBetter = m_score > r_score + l_score + v_score;
  double total_score = parentBetter ? m_score : r_score + l_score + v_score;


  // If one is a split and one is a leaf
  else if (split[0] < 0 || split[1] < 0) {
    const bool isLeaf = split[0] < 0;
    int leaf = isLeaf ? split[0] : split[1];
    int split_node = isLeaf ? split[1] : split[0];

    // Get the list of the best clusters for the split node
    branch_info current_cl = hcl_DFS(split_node-1, merge, pt_map, f, cl, score_map, propagate);

    // Compute the score. If propagating the score through the hierarchy, save the final
    // max score under the current branch node in the score map. Otherwise, 'refresh' the
    // score by computing the raw score for all of the points at every level.
    bool parentBetter = false;
    double total_score = 0;
    if (propagate && score_map.find(split_node) != score_map.end()){
      double s_score = score_map.at(split_node);
      double l_score = as<double>(f(l_ids));
      double v_score = as<double>(f(v_ids));
      double m_score = as<double>(f(m_ids));
      parentBetter = m_score > s_score + l_score + v_score;
      total_score = parentBetter ? m_score : s_score + l_score + v_score;
    } else {
      parentBetter = m_score > s_score + l_score + v_score;
      total_score = parentBetter ? m_score : s_score + l_score + v_score;
    }

    // Make the choice to keep the parent cluster id or the children
    if (parentBetter){
      cl.push_back(k-1); // Keep the parent (current) branch
      cl.remove(leaf); // Remove leaf from list of candidates
      for(std::list<int>::iterator cid = current_cl.ids.begin(); cid != current_cl.ids.end(); ++cid){
        cl.remove(*cid); // Remove all children ids from the current best clustering list
      }
      score_map[k-1] = total_score;
    } else {
      // Do nothing, children already in the best list
      score_map[k-1] = total_score; // Always update score map
    }

    // Update the current branch
    current_cl.add(leaf, pt_map);
    return current_cl;
  }
  // If both were split branches
  else {
    int left = split[0], right = split[1];


    IntegerVector v_ids = pt_map[patch::to_string(k)];
    IntegerVector m_ids = merge3(l_ids, v_ids, s_ids);

    // Calculate cluster scores
    double l_score = as<double>(f(l_ids));
    double r_score = as<double>(f(s_ids));
    double v_score = as<double>(f(v_ids));
    double m_score = as<double>(f(m_ids));

    bool parentBetter = m_score > l_score + r_score + v_score;
    double total_score = parentBetter ? m_score : l_score + r_score + v_score;

    // Make the choice to keep the parent cluster id or the children
    if (parentBetter){
      cl.push_back(k); // Keep the parent (current) branch
      for(std::list<int>::iterator cid = cl_l.ids.begin(); cid != cl_l.ids.end(); ++cid){
        cl.remove(*cid); // Remove all children ids from the current best clustering list
      }
      for(std::list<int>::iterator cid = cl_r.ids.begin(); cid != cl_r.ids.end(); ++cid){
        cl.remove(*cid); // Remove all children ids from the current best clustering list
      }
      score_map[k] = total_score;
      Rcout << "keeping parent branch." << std::endl;
    } else {
      // Do nothing, children already in the best list
      score_map[k] = total_score; // Always update score map
      Rcout << "keeping upper branches." << std::endl;
    }
    cl_l.merge(cl_r);
    return cl_l;
  }
}


// Framework for the Optimal Selection of Clusters (FOSC)
IntegerVector fosc(const List& hc, Rcpp::Function f) {
  const IntegerMatrix& merge = hc["merge"];
  const NumericVector& height = hc["height"];
  const List& pt_map = hc["idx"];

  const int n = merge.nrow() + 1;
  std::list<int> cl = std::list<int>(); // cluster output list
  std::unordered_map<int, double> score_map = std::unordered_map<int, double>();
  hcl_DFS(n - 2, merge, pt_map, f, cl, score_map);
  return(wrap(cl));
}

/*** R

library("clustertree")
sl <- clustertree::clustertree(iris[, 1:4], k = 5L)
sl_si <- clustertree:::simplified_hclust(sl$hc, 2)

sf <- function(cl) { if(length(cl) == 0) return(0) else return(max(cl)) }
# what <- clustertree:::fosc(sl_si, sf)
*/
