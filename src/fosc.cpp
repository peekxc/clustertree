#include "RcppHeader.h"
#include "utilities.h"
#include "union_find.h"

int global_cid = 1;

// MERGE := a valid merge matrix
// F1 := Function to apply
// S1 := Stability score for left branch
// S2 := Stability score for right branch
// DEFINE FOSC(CID, MERGE, F1):
//  IF ( IS_LEAF(CID) ):
//    RETURN F1(CID)
//  ELSE:
//    LEFT_SCORE = FOSC(MERGE[CID]->LEFT, F1)
//    RIGHT_SCORE = FOSC(MERGE[CID]->RIGHT, F1)
//    VIRTUAL_SCORE = F1( BRANCH_IDS )
//  MERGE_SCORE = F1(  )
//  SPLIT_SCORE = AGGREGATE( LEFT_SCORE, RIGHT_SCORE, VIRTUAL_SCORE )
//  RETURN ( MAX(SCORE)S)
//
#include <Rcpp.h>
using namespace Rcpp;

void printVec(const IntegerVector& vec){
  for (IntegerVector::const_iterator it = vec.begin(); it != vec.end(); ++it){ Rcout << *it << ", "; }
  Rcout << std::endl;
}

// a predicate implemented as a class:
struct is_in {
  std::list<int> ids;
  bool operator() (const int& value) { return (std::find(ids.begin(), ids.end(), value) != ids.end()); }
};

// Merge 3 integer vectors
IntegerVector merge3(const IntegerVector& b1, const IntegerVector& b2, const IntegerVector& b3){
  IntegerVector out = no_init(b1.size() + b2.size() + b3.size());
  if (b1.size() > 0) std::copy(b1.begin(), b1.end(), out.begin());
  if (b2.size() > 0) std::copy(b2.begin(), b2.end(), out.begin() + b1.size());
  if (b3.size() > 0) std::copy(b3.begin(), b3.end(), out.begin() + b1.size() + b2.size());
  return out;
}


// Simple struct to keep track of branch information through the recursive calls
struct branch {
  double score; // The branchs quality 'score'
  std::list<int> ids; // list to track the cluster ids in branch (acculumative)
  IntegerVector points; // continuously update vector to track point ids
  branch(){
    score = 0.0;
    ids = std::list<int>();
    points = IntegerVector();
  }
  branch(IntegerVector& pt_ids, double _score){
    ids = std::list<int>();
    points = pt_ids;
    score = _score;
  }

  branch(std::list<int>& cids, IntegerVector& pt_ids, double _score){
    ids = cids;
    points = pt_ids;
    score = _score;
  }
  ~branch(){ }
  branch(const branch& rhs) : ids(rhs.ids), points(rhs.points){}
  // branch& operator=(const branch& rhs){
  //   ids = rhs.ids;
  //   points = rhs.points;
  //   return *this;
  // }

  // Add the cluster ids to the branch
  void add_clusters(std::list<int>& cids){
    for(std::list<int>::iterator cid = cids.begin(); cid != cids.end(); ++cid){
      this->ids.push_back(*cid);
    }
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

  // Add the actual point indices directly to the pont array
  void add(IntegerVector pt_ids){
    for (IntegerVector::iterator pid = pt_ids.begin(); pid != pt_ids.end(); ++pid){
      points.push_back(*pid);
    }
  }

  // Merge a branch cluster and point indices into this one
  void merge(branch& other_branch){
    ids.insert(ids.end(), other_branch.ids.begin(), other_branch.ids.end());
    std::copy(other_branch.points.begin(), other_branch.points.end(), points.end());
  }
};

// Given a index, traverse the merge matrix, unioning leaves until arriving at that index.
// Once at the given index, break out of the loop return the point indices part of the
// local branch as if it were cut there.
// [[Rcpp::export]]
IntegerVector getMergeIndices(const IntegerMatrix& merge, const int m_i){
  if (m_i < 0) stop("'getMergeIndices' expects a positive index.");
  if (m_i > merge.nrow() + 1) stop("'getMergeIndices' expects an index less than the size of the data set.");
  const int n = merge.nrow() + 1;
  UnionFind ds = UnionFind(n);
  NumericVector merge_height = NumericVector(n);

  // Component tracker
  std::vector<int> c_idx = std::vector<int>(n-1, 0);

  // Traverse merge
  int i;
  for (i = 0; i < n -1; ++i){
    int lm = merge(i, 0), rm = merge(i, 1);
    IntegerVector m = IntegerVector::create(lm, rm);
    if (all(m < 0).is_true()){
      int left_id = (-lm)-1, right_id = (-rm)-1;
      ds.Union(left_id, right_id); // Union the leaves
      c_idx[i] = ds.Find(left_id); // record their component
    } else if (any(m < 0).is_true()) {
      int leaf = lm < 0 ? (-lm)-1 : (-rm)-1;
      int comp = lm < 0 ? rm - 1 : lm - 1;
      ds.Union(leaf, c_idx[comp]); // Merge leaf with the component
      c_idx[i] = ds.Find(leaf);
    } else {
      int left_id = c_idx[lm - 1], right_id = c_idx[rm - 1];
      ds.Union(ds.Find(left_id), ds.Find(right_id));
      c_idx[i] = ds.Find(left_id);
    }
    if (i == m_i - 1){ break; }
  }
  const int cid = c_idx[i];
  IntegerVector cl = ds.getCC();
  IntegerVector idx = which_cpp(cl, cid);
  return(idx + 1);
}

// f1 :=
// f2 := merged function
struct fosc {
  const IntegerMatrix& merge;
  Rcpp::Function f1;
  Rcpp::Function f2;
  std::list<int> cluster_ids;
  List idx_map;
  NumericVector b_scores;
  IntegerVector cluster_out;
  fosc(const IntegerMatrix& m, Rcpp::Function f_1, Rcpp::Function f_2, List& pt_map)
    : merge(m), f1(f_1), f2(f_2) {
    cluster_ids = std::list<int>(); // keep 1-based
    idx_map = pt_map; // 1-based
    b_scores = NumericVector(m.nrow());
    cluster_out = IntegerVector(m.nrow() + 1);
  }

  // void removeClusters(){
  //
  // }

  // Retrieve the point indices of the given component index. If the id given
  // is negative, retrieve the singleton as a positive index. All assumed 1-based.
  IntegerVector getPointIdx(const int id){
    std::string key_str = patch::to_string(id);
    if ((bool) idx_map.containsElementNamed(key_str.c_str())){
      return as<IntegerVector>(idx_map[key_str]);
    } else if (id < 0){
      int singleton = abs(id);
      return IntegerVector::create(singleton);
    } else { return IntegerVector(); }
  }

  // Expects k to be 0-based
  branch recurse(int k){

    // If k < 0, it's a leaf
    if (k < 0) {
      // Get the current set of point ids
      IntegerVector point_ids = getPointIdx(k + 1);

      // Create a new branch for the leaf
      branch leaf = branch();
      leaf.ids.push_back(k + 1); // should equal leaf index
      leaf.points = point_ids;
      leaf.score = as<double>(f1(leaf.points));

      // Add the leaf to the cluster id list
      cluster_ids.push_back(k + 1); // 1-based
      return leaf;
    }
    else {
      const int left = merge(k, 0) - 1, right = merge(k, 1) - 1; // 0-based

      // Recurse Left & Right
      branch left_br = recurse(left), right_br = recurse(right); // 0-based
      Rcout << "Left: " << left + 1 << ", Right: " << right + 1 << std::endl;

      // Virtual branch
      IntegerVector current_ids = getPointIdx(k + 1); // 0-based pass
      branch v_br = branch(current_ids, current_ids.size() == 0 ? 0 : as<double>(f1(current_ids)));

      // Create the merged branch, moving back up
      branch merged_branch = branch();
      IntegerVector m_pts = merge3(left_br.points, right_br.points, v_br.points);
      double m_score = m_pts.size() == 0 ? 0 : as<double>(f1(m_pts));
      double s_score = m_pts.size() == 0 ? 0 : as<double>(f2(left_br.points, right_br.points, v_br.points));
      // Rcout << "Merged Score: " << m_score << std::endl;
      // Rcout << "Split Score: %.12f" << s_score << std::endl;
      Rprintf("Merged Score: %.14f\n", m_score);
      Rprintf("Split Score: %.14f\n", s_score);
      Rcout << "Merge == Split: " << int(m_score == s_score) << ", (" << fabs(m_score - s_score) << ")" <<  std::endl;
      Rcout << "Precision: " << std::numeric_limits<double>::epsilon() << std::endl;
      Rcout << "Under precision: " << int(fabs(m_score - s_score) < std::numeric_limits<double>::epsilon()) << std::endl;

      merged_branch.ids = std::list<int>();
      if (m_score >= s_score || fabs(m_score - s_score) < std::numeric_limits<double>::epsilon()) {
        merged_branch.ids.push_back(k + 1); // 1-based

        // Update clusters
        Rcout << "Use upper: (pushing " << k + 1 << ")" << std::endl;
        cluster_ids.push_back(k + 1); // 1-based

        Rcout << "Current List: ";
        for (std::list<int>::iterator it = cluster_ids.begin(); it != cluster_ids.end(); ++it){
          Rcout << *it << ", ";
        }
        Rcout << std::endl;


        for (std::list<int>::iterator it = left_br.ids.begin(); it != left_br.ids.end(); ++it){
          Rcout << "(L) Removing " << (*it)  << std::endl;
          cluster_ids.remove(*it);
        }
        for (std::list<int>::iterator it = right_br.ids.begin(); it != right_br.ids.end(); ++it){
          Rcout << "(R) Removing " << (*it) << std::endl;
          cluster_ids.remove(*it);
        }

        Rcout << "Current List: ";
        for (std::list<int>::iterator it = cluster_ids.begin(); it != cluster_ids.end(); ++it){
          Rcout << *it << ", ";
        }
        Rcout << std::endl;

      } else {
        // Add the cluster indices from the left and right branches to the merged branch
        Rcout << "Use lower" << std::endl;
        for (std::list<int>::iterator it = left_br.ids.begin(); it != left_br.ids.end(); ++it){ merged_branch.ids.push_back(*it); }
        for (std::list<int>::iterator it = right_br.ids.begin(); it != right_br.ids.end(); ++it){ merged_branch.ids.push_back(*it); }
      }
      // Always pick the higher score moving up
      merged_branch.score = std::max(m_score, s_score); // Set stability score for this branch
      b_scores[k] = merged_branch.score;

      // Save the points
      merged_branch.points = m_pts;
      Rcout << "Current List: ";
      for (std::list<int>::iterator it = cluster_ids.begin(); it != cluster_ids.end(); ++it){
        Rcout << *it << ", ";
      }
      Rcout << std::endl;
      //Rcpp::stop("depth here");
      return(merged_branch);
    }
  }

};

// [[Rcpp::export]]
List FOSC_int(const IntegerMatrix& merge, Rcpp::Function f_1, Rcpp::Function f_2){
  List null_map = List::create();
  const int n = merge.nrow() + 1;
  fosc hierarchy = fosc(merge, f_1, f_2, null_map);
  branch root = hierarchy.recurse(n - 2);
  Rcout << "Final Score: " << root.score << std::endl;
  List res = List::create(_["cl"] = wrap(hierarchy.cluster_ids),
                          _["cids"] = wrap(root.ids),
                          _["pids"] = wrap(root.points),
                          _["scores"] = hierarchy.b_scores,
                          _["cl_out"] = hierarchy.cluster_out);
  return(res);
}



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
// branch* hcl_DFS(const int k, const IntegerMatrix& merge, const List& pt_map,
//                      Rcpp::Function f1, Rcpp::Function f2, // Scoring Functions
//                      std::list<int>& cl,
//                     std::unordered_map<int, double>& score_map, bool propagate = true){
//
//   // The current points
//   IntegerVector split = merge(k, _);
//
//   // Gather the (extracted) point indices
//   int left = split[0], right = split[1];
//   IntegerVector v_ids, m_ids; // left, right, virtual, and merged ids
//   //std::list<int> l_cids, r_cids;
//   branch* l_branch, *r_branch;
//
//   // Cases
//   if (all(split < 0).is_true()){
//
//     // Add point indices and cluster ids
//     l_branch = new branch(std::list<int>(1, left), pt_map[])
//     l_branch->add(left, pt_map);
//     l_branch->add_clusters();
//
//     // Right branch
//     r_branch->add(right, pt_map);
//     r_branch->add_clusters(std::list<int>(1, right));
//
//     Rcout << "Adding clusters: " << left << ", " << right << std::endl;
//     cl.push_back(left); // Assume the leaf is in the final result
//     cl.push_back(right); // Assume the leaf is in the final result
//
//   } else if (left < 0 || right < 0) {
//     const bool isLeaf = left < 0;
//     int leaf = isLeaf ? left : right;
//     int split_node = isLeaf ? right : left;
//     cl.push_back(leaf); // Assume the leaf is in the final result
//     Rcout << "Adding clusters: " << leaf << std::endl;
//
//     // Recurse down the split node
//     branch* current_cl = hcl_DFS(split_node-1, merge, pt_map, f1, f2, cl, score_map, propagate);
//
//     // Assign left and right ids
//     l_branch->points = isLeaf ? pt_map[patch::to_string(left)] : current_cl->points;
//     r_branch->points = isLeaf ? current_cl->points : pt_map[patch::to_string(right)];
//
//     // Assign cluster ids
//     l_branch->ids = isLeaf ? std::list<int>(1, left) : current_cl->ids;
//     r_branch->ids = isLeaf ? current_cl->ids : std::list<int>(1, right);
//
//     // No longer need that branch
//     delete current_cl;
//   } else {
//     branch* current_cl1 = hcl_DFS(left-1, merge, pt_map, f1, f2, cl, score_map, propagate);
//     branch* current_cl2 = hcl_DFS(right-1, merge, pt_map, f1, f2, cl, score_map, propagate);
//     l_ids = current_cl1->points;
//     r_ids = current_cl2->points;
//     l_cids = current_cl1->ids;
//     r_cids = current_cl2->ids;
//     delete current_cl1;
//     delete current_cl2;
//   }
//
//   // Create a new set of branch information
//   branch* cb = new branch(); // current branch
//   if (k != merge.nrow() - 1){ // Don't evaluate root as a possibility
//
//     // Compute virtual point ids (if any) and total merged ids
//     v_ids = IntegerVector();
//     if ((bool) pt_map.containsElementNamed(patch::to_string(k+1).c_str())){
//       v_ids = pt_map[patch::to_string(k+1)];
//     }
//     m_ids = merge3(l_ids, v_ids, r_ids);
//
//     // Populate the current branch information
//     cb->add(m_ids); // point indices
//     cb->add_clusters(l_cids); // cluster ids in the left branch
//     cb->add_clusters(r_cids); // cluster ids in the right branch
//
//     // Calculate cluster scores
//     double ch_score = as<double>(f2(l_ids, r_ids, v_ids)); // Children (split) score
//     double m_score = as<double>(f1(m_ids)); // Merged score (all points)
//
//     Rcout << "here3" << std::endl;
//
//     // Determine whether parent is better or not; compute final score
//     bool parentBetter = m_score >= std::max(cb->score ch_score;
//
//     Rcout << "split: " << split[0] << ", " << split[1] << " ";
//     Rcout << "Parent score: " << m_score << " < vs. " << ch_score << ">" << std::endl;
//     Rcout << "Parent better ?" << (parentBetter ? "True" : "False") << std::endl;
//     Rcout << "Cluster ids (component " << k+1 << "): ";
//     for (std::list<int>::iterator cid = cb->ids.begin(); cid != cb->ids.end(); ++cid){
//       Rcout << *cid << ", ";
//     }
//     Rcout << std::endl;
//     // Test condition: If the current merge is better, push back the current branches merged component
//     // id. Otherwise, add the left and right branch ids to the current clustering list
//     if (parentBetter){
//       cl.push_back(k+1); // Keep the parent
//       Rcout << "Adding cluster: " << k+1 << std::endl;
//
//       // Remove all children ids from the current best clustering list
//       for(std::list<int>::iterator cid = cb->ids.begin(); cid != cb->ids.end(); ++cid){
//         Rcout << "1. Removing: " << *cid << std::endl;
//         cl.remove(*cid); // Remove from the final clustering result any child cluster ids
//       }
//     } else {
//     // The parent is not better;
//       // Remove all children ids from the current best clustering list
//       // for(std::list<int>::iterator cid = cb->ids.begin(); cid != cb->ids.end(); ++cid){
//       //   if (*cid != left && *cid != right){
//       //     Rcout << "2. Removing: " << *cid << std::endl;
//       //     cl.remove(*cid); // if there's a child cluster id not in the current split, remove it
//       //   }
//       // }
//     }
//     cb->score = parentBetter ? m_score : ch_score;
//     return cb; // Return the newly formed current branch
//   } else {
//     // We've arrived back at the root; cleanup
//     delete cb;
//   }
//   return NULL;
// }
//

// Framework for the Optimal Selection of Clusters (FOSC)
// IntegerVector fosc(const List& hc, Rcpp::Function f1, Rcpp::Function f2) {
//   const IntegerMatrix& merge = hc["merge"];
//   const NumericVector& height = hc["height"];
//   const List& pt_map = hc["idx"];
//
//   const int n = merge.nrow() + 1;
//   std::list<int> cl = std::list<int>(); // cluster output list
//
//   // Create a map for each node in the tree (internal and leaf)
//   std::unordered_map<int, branch> b_map = std::unordered_map<int, branch>();
//   hcl_DFS(n - 2, merge, pt_map, f1, f2, cl, score_map);
//   return(wrap(cl));
// }

/*** R

idx <- sample(1:nrow(iris), size = 25)
sl <- hclust(dist(iris[idx, 1:4]), method = "single")
X <- as.integer(iris[idx, 5])

## Truth
plot(dendextend::color_branches(as.dendrogram(sl), clusters = X))

## ARI
# f1 <- function(x) { mclust::adjustedRandIndex(rep(1L, length(x)), X[x]) }
# f2 <- function(x_i, x_j, x_k) { mean(f1(x_i) + f1(x_j) + f1(x_k)) }

## Cluster Purity
n <- length(X)
f1 <- function(x) {
  if (length(x) == 0) return(0)
  else { Ct <- table(X[x]); t_j <- Ct[which.max(Ct)]; return(as.numeric(t_j/n)) }
}
f2 <- function(x_i, x_j, x_k) { mean(f1(x_i) + f1(x_j) + f1(x_k)) }

what <- clustertree:::FOSC_int(sl$merge, f1, f2)

cl_out <- rep(0L, length(X))
cid <- 1L
for (comp_id in what$cl){
  if(comp_id > 0){
    leaf_idx <- clustertree:::getMergeIndices(merge = sl$merge, comp_id)
    cl_out[which(sl$order %in% leaf_idx)] <- cid
    cid <- cid + 1L
  }
}

## FOSC w/ Cluster Purity
plot(dendextend::color_branches(as.dendrogram(sl), clusters = cl_out))
#lapply(what$cl, function(cid) clustertree:::getMergeIndices(merge = sl$merge,  cid))

#sapply(sl$height, function(h) cutree(sl, h = h))[what$cl, ]
#cl <- match(what$cl_out, sort(unique(what$cl_out)))


# Truth
plot(dendextend::color_branches(as.dendrogram(sl), clusters = X))

#
cl_out <- rep(0L, length(X))
cl_out[sl$order][118] <- 1L
cl_out[sl$order][132] <- 1L




# library("clustertree")
# sl <- clustertree::clustertree(iris[, 1:4], k = 5L)
# sl_si <- clustertree:::simplified_hclust(sl$hc, 2)
#
# sf_1 <- function(cl) { if(length(cl) == 0) return(0) else return(max(cl)) }
# sf_2 <- function(l, r, v) { cl <- c(l, r, v); if(length(cl) == 0) return(0) else return(max(cl)) }
#
#
# ARI <- function(X, Y){
#   n <- length(X)
#   cross <- unname(as.matrix(table(X, Y)))
#   { a <- apply(cross, 1, sum); b <- apply(cross, 2, sum) }
#   index <- sum(choose(cross, 2))
#   exp_idx <- (sum(choose(a, 2)) * sum(choose(b, 2)))/choose(n, 2)
#   max_idx <- 0.5 * (sum(choose(a, 2)) + sum(choose(b, 2)))
#   return ((index - exp_idx) / (max_idx - exp_idx))
# }
# s1 <- compiler::cmpfun(ARI)
#
# truth_cl <- hdbscan(moons, minPts = 5L)$cluster
# sf_1 <- function(cl) {
#   if(length(cl) == 0)
#     return(0)
#   else {
#     s1(truth_cl[cl], rep(1L, length(cl)))
#   }
# }
#
# sf_2 <- function(l, r, v) {
#   cl <- c(l, r, v)
#   if(length(cl) == 0) return(0)
#   else {
#     un_cl <- c(rep(1L, length(l)), rep(2L, length(r)), rep(3L, length(v)))
#     return(s1(truth_cl[cl], un_cl))
#   }
# }
#
# hc <- hclust(dist(moons), method = "single")
# hc$idx <- as.list(structure(1:nrow(moons), names = as.character(-(1:nrow(moons)))))
# what <- clustertree:::fosc(hc, sf_1, sf_2)
#
# cl <- rep(0L, nrow(iris))
# cl[sl_si$idx$`-7`] <- 1L
# cl[sl_si$idx$`-1`] <- 2L
# cl[sl_si$idx$`5`] <- 3L
# plot(iris[, 1:2], col = cut(sl_si, k = 2)+1L)
#
#
# t1 <- clustertree::clustertree(moons, k = 5L)
# t1$hc$idx <- as.list(structure(1:nrow(moons), names = as.character(-(1:nrow(moons)))))
# clustertree:::fosc(t1$hc, sf_1, sf_2)

*/
