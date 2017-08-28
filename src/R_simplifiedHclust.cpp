#include <Rcpp.h>
using namespace Rcpp;

#include "unordered_map"
#include "union_find.h" // disjoint set data structure
#include "hclust_util.h" // mstToHclust, various hclust tools
#include "utilities.h" // various helpful utilities for working with Rcpp

// Simple struct to hold simplified cluster info
struct cl_info{
  IntegerVector contains;
  NumericVector eps;
  double eps_death, eps_birth;
  bool processed;
  int n_children;
  IntegerVector children;
  cl_info(){
    contains = IntegerVector();
    eps = NumericVector();
    processed = false;
    eps_death = eps_birth = n_children = 0;
  };
};

// Populate the from, to, and height vectors needed to make a hierarchy
int genSimMerges(int cid, std::unordered_map<int, cl_info>& cl,
                 IntegerVector& from, // from indices
                 IntegerVector& to, // to indices
                 IntegerVector& level_idx, // record the comp index from and to merged into
                 NumericVector& height, // Height of the tree
                 std::map<std::pair<int, int>, int>& parent_map, // map of child pairs to parent id
                 List& cl_hierarchy,
                 int level  // In the case of ties in height, prioritize by level to ensure mappings work out
){
  if (cl[cid].n_children == 0){
    return -cid;
  } else {
    // Get the left and right child cluster ids
    int left_id = cl[cid].children.at(0), right_id = cl[cid].children.at(1);
    cl_hierarchy[patch::to_string(cid)] = IntegerVector::create(left_id, right_id);

    // Store the parent, indexed by the children! This is necessary for later.
    parent_map.insert(std::make_pair(std::minmax(left_id, right_id), cid));

    // Recursively go left and right
    int left_pt = genSimMerges(left_id, cl, from, to, level_idx, height, parent_map, cl_hierarchy, level+1);
    int right_pt = genSimMerges(right_id, cl, from, to, level_idx, height, parent_map, cl_hierarchy, level+1);
    double cheight = std::max(cl[std::abs(left_pt)].eps_birth, cl[std::abs(right_pt)].eps_birth);

    from.push_back(left_pt);
    to.push_back(right_pt);
    level_idx.push_back(level);

    // Trick
    // while(any(height == cheight).is_true()){
    //   cheight += std::numeric_limits<double>::epsilon();
    // }
    height.push_back(cheight);


    // Return arbitrary point index; let the disjoint set handle the tracing of the indices later on
    return cid;
  }
}

// [[Rcpp::export]]
IntegerVector cut_simplified_hclust(List hcl, IntegerVector cl_in, const int big_n){
  IntegerMatrix merge = hcl["merge"];
  const int n = merge.nrow() +1;
  List idx = hcl["idx"];

  IntegerVector cl_out = Rcpp::no_init(big_n);
  std::vector<int> cont = std::vector<int>(n-1, 0);
  int left_cid =0, right_cid = 0;

  for (int i=0; i < n-1; ++i){
    int lm = merge(i, 0), rm = merge(i, 1);
    IntegerVector m = IntegerVector::create(lm, rm);

    if (all(m < 0).is_true()){
      left_cid = cl_in[(-lm)-1], right_cid = cl_in[(-rm)-1];
      IntegerVector left_ids = idx[patch::to_string(lm)];
      IntegerVector right_ids = idx[patch::to_string(rm)];
      cl_out[left_ids - 1] = left_cid;
      cl_out[right_ids - 1] = right_cid;
      cont[i] = left_cid == right_cid ? left_cid : 0;
    } else if (any(m < 0).is_true()){
      int leaf = m.at(0) < 0 ? m.at(0) : m.at(1);
      int comp = m.at(0) < 0 ? m.at(1) : m.at(0);
      int leaf_cid = cl_in[(-leaf)-1];
      IntegerVector leaf_ids = idx[patch::to_string(leaf)];
      cl_out[leaf_ids -1] = leaf_cid;
      IntegerVector comp_ids = idx[patch::to_string(comp)];
      left_cid = right_cid = leaf_cid;
      if (cont[comp-1] == leaf_cid){
        cl_out[comp_ids-1] = leaf_cid;
        cont[i] = leaf_cid;
      } else {
        cl_out[comp_ids-1] = 0;
        cont[i] = 0; // Mark this component as noise
      }
    } else {
      left_cid = cont[m.at(0) - 1], right_cid = cont[m.at(1) - 1];
      IntegerVector left_ids = idx[patch::to_string(m.at(0))];
      IntegerVector right_ids = idx[patch::to_string(m.at(1))];
      cl_out[left_ids-1] = left_cid;
      cl_out[right_ids-1] = right_cid;
      if (left_cid == right_cid){
        cont[i] = left_cid;
      } else {
        cont[i] = 0;
      }
    }
  }
  // Last step: compare most recent left and right clusters to determine what root labels
  IntegerVector root_ids = idx["0"];
  //Rcout << "root ids: " << root_ids.at(0) << std::endl;
  if (left_cid == right_cid){
    cl_out[root_ids - 1] = left_cid;
  } else {
    cl_out[root_ids - 1] = 0;
  }

  return(cl_out);
}


// Given an hclust object and a minimum cluster size, traverse the tree divisely to create a
// simplified hclust object, where each leaf contains at least min_sz points
// [[Rcpp::export]]
List simplified_hclust(List hcl, const int min_sz) {
  // Extract hclust info
  NumericMatrix merge = hcl["merge"];
  NumericVector eps_dist = hcl["height"];
  IntegerVector pt_order = hcl["order"];
  int n = merge.nrow() + 1, k;

  //  Auxilary information not saved
  std::vector<int> cl_tracker = std::vector<int>(n-1 , 0); // cluster component each step
  std::vector<int> member_sizes = std::vector<int>(n-1, 0); // Size each step

  // Primary information needed
  std::unordered_map<int, cl_info> cl = std::unordered_map<int, cl_info>();
  cl.reserve(n);

  // First pass: Hclust object are intrinsically agglomerative structures. Splitting the hierarchy
  // divisively requires knwoledge of the sizes of each branch. So agglomerate up the hierarchy, recording member sizes.
  // This enables a dynamic programming strategy to improve performance below.
  for (k = 0; k < n-1; ++k){
    int lm = merge(k, 0), rm = merge(k, 1);
    IntegerVector m = IntegerVector::create(lm, rm);
    if (all(m < 0).is_true()){
      member_sizes[k] = 2;
    } else if (any(m < 0).is_true()) {
      int pos_merge = (lm < 0 ? rm : lm), merge_size = member_sizes[pos_merge - 1];
      member_sizes[k] = merge_size + 1;
    } else {
      // Record Member Sizes
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];
      member_sizes[k] = merge_size1 + merge_size2;
    }
  }

  // Initialize root (unknown size, might be 0, so don't initialize length)
  int root_id = 0, global_cid = 0;
  cl[root_id].eps_birth = eps_dist.at(eps_dist.length()-1);

  // Second pass: Divisively split the hierarchy, recording the epsilon and point index values as needed
  for (k = n-2; k >= 0; --k){
    // Current Merge
    int lm = merge(k, 0), rm = merge(k, 1), cid = cl_tracker.at(k);
    IntegerVector m = IntegerVector::create(lm, rm);

    // Get reference to the current cluster info
    cl_info& info = cl[cid]; // construct one using default constructor if doesn't exist

    // Trivial case: split into singletons, record eps, contains, and ensure eps_death is minimal
    if (all(m < 0).is_true()){
      info.contains.push_back(-lm), info.contains.push_back(-rm);
      double noise_eps = info.processed ? info.eps_death : eps_dist.at(k);
      info.eps.push_back(noise_eps), info.eps.push_back(noise_eps);
      if (!info.processed && eps_dist.at(k) < info.eps_death){
        info.eps_death = eps_dist.at(k);
      }
    } else if (any(m < 0).is_true()) {
      // Record new point info and mark the non-singleton with the cluster id
      info.contains.push_back(-(lm < 0 ? lm : rm));
      info.eps.push_back(info.processed ? info.eps_death : eps_dist.at(k));
      cl_tracker.at((lm < 0 ? rm : lm) - 1) = cid;
    } else {
      int merge_size1 = member_sizes[lm-1], merge_size2 = member_sizes[rm-1];

      // The cluster-simplification step: if two daughter nodes have 'runt sizes' at least as
      // big as the threshold given, consider this as a true split.
      if (merge_size1 >= min_sz && merge_size2 >= min_sz){
        // Record death of current cluster
        info.eps_death = eps_dist.at(k);
        info.processed = true;

        // Mark the lower merge steps as new clusters
        info.children = IntegerVector::create(global_cid+1, global_cid+2);
        int l_index = global_cid+1, r_index = global_cid+2;
        cl_tracker.at(lm - 1) = ++global_cid, cl_tracker.at(rm - 1) = ++global_cid;

        // Record the distance the new clusters appeared and initialize containers
        cl_info& left_node = cl[l_index], &right_node = cl[r_index];
        left_node.eps_birth = eps_dist.at(k), right_node.eps_birth = eps_dist.at(k);
        left_node.eps_death = eps_dist.at(lm - 1), right_node.eps_birth = eps_dist.at(rm - 1);
        left_node.processed = false, right_node.processed  = false;
        info.n_children = merge_size1 + merge_size2;
      } else {
        // Inherit cluster identity
        cl_tracker.at(lm - 1) = cid,  cl_tracker.at(rm - 1) = cid;
      }
    }
  }

  // Point ids of original observations
  List point_idx = List();

  // The tree has not been divisely split into a simplified hierarchy, but now again,
  // to create an 'hclust' object, the hierarchical structure must be built agglomeratively.
  IntegerVector from = IntegerVector(), to = IntegerVector(), level_idx = IntegerVector();
  NumericVector height = NumericVector();
  std::map<std::pair<int, int>, int> parent_map = std::map<std::pair<int, int>, int>(); // need to record the parents of pairwise children
  List cl_hierarchy = List();
  genSimMerges(0, cl, from, to, level_idx, height, parent_map, cl_hierarchy, 0);

  // Create a map from original cluster ids to normalized ones
  IntegerVector leaf_from = as<IntegerVector>(from[from < 0]), leaf_to = as<IntegerVector>(to[to < 0]);
  IntegerVector all_ids = unique(combine(leaf_from, leaf_to)); // original leaf ids
  std::transform(all_ids.begin(), all_ids.end(), all_ids.begin(), ::abs);
  std::unordered_map<int, int> singleton_map = std::unordered_map<int, int>();
  int i = 1;
  for (IntegerVector::iterator it = all_ids.begin(); it != all_ids.end(); ++it, ++i){
    singleton_map[*it] = i;
    point_idx[patch::to_string(-singleton_map[*it])] = cl[*it].contains;
  }

  // Reorder the indices by height
  IntegerVector merge_order = order_(height) - 1; // 0-based indices
  NumericVector o_height = height[merge_order];
  IntegerVector o_from = from[merge_order], o_to = to[merge_order];

  // Get index of the true (original) component indices formed at every merge level
  const int merge_sz = o_height.size();
  std::vector<int> comp_idx = std::vector<int>(merge_sz);
  std::unordered_map<int, int> comp_idx2 = std::unordered_map<int, int>(merge_sz);
  IntegerMatrix new_merge = IntegerMatrix(merge_sz, 2);
  for(i = 0; i < merge_sz; ++i){
    IntegerVector pt_ids = IntegerVector::create(o_from.at(i), o_to.at(i));
    if (all(pt_ids < 0).is_true()){
      int left = singleton_map[-pt_ids.at(0)];
      int right = singleton_map[-pt_ids.at(1)];
      new_merge(i, _) = NumericVector::create(-left, -right);
      int orig_comp = parent_map.at(std::minmax(-pt_ids.at(0), -pt_ids.at(1)));
      comp_idx[i] = orig_comp;
      comp_idx2[orig_comp] = i;
      point_idx[patch::to_string(i)] = cl[orig_comp].contains;
    } else if (any(pt_ids < 0).is_true()){
      int leaf = pt_ids.at(0) < 0 ? -pt_ids.at(0) : -pt_ids.at(1);
      int comp = pt_ids.at(0) < 0 ? pt_ids.at(1) : pt_ids.at(0);
      int orig_comp = parent_map.at(std::minmax(leaf, comp));
      comp_idx[i] = orig_comp;
      comp_idx2[orig_comp] = i;
      new_merge(i, _) = NumericVector::create(-singleton_map[leaf], comp_idx2[comp]+1);
      point_idx[patch::to_string(i)] = cl[orig_comp].contains;
    } else {
      int orig_comp = parent_map.at(std::minmax(pt_ids.at(0), pt_ids.at(1)));
      comp_idx[i] = orig_comp;
      comp_idx2[orig_comp] = i;
      new_merge(i, _) = NumericVector::create(comp_idx2[pt_ids.at(0)]+1, comp_idx2[pt_ids.at(1)]+1);
      point_idx[patch::to_string(i)] = cl[orig_comp].contains;
    }
  }

  List simple_hclust = List::create(
    _["merge"] = new_merge,
    _["height"] = o_height,
    _["order"] = extractOrder(new_merge),
    _["labels"] = R_NilValue
    //_["mst"] = mst
  );
  simple_hclust["idx"] = point_idx;

  // Prepend "simplified_hclust" to the S3 class name and return
  std::vector<std::string> class_names = std::vector<std::string>(2);
  class_names[0] = "simplified_hclust"; // Prefer simplified S3 methods
  class_names[1] = "hclust"; // It still IS-A valid hclust object
  simple_hclust.attr("class") = wrap(class_names);
  return(simple_hclust);
}