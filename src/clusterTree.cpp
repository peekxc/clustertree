#include <Rcpp.h>
using namespace Rcpp;

// Includes
#include "union_find.h"
#include "utilities.h"
#include <unordered_map>
#include <stack>

static int global_cid = 0;
static int counter = 0; // global counter

// Computes the connection radius, i.e. the linkage criterion, required to
// connect edges. If the there is some saliency criteria preventing connectivity
// to improve the robustness of the estimator (such as a minimum kNN distance criteria),
// returns the maximum floating point representation. Note that infinity cannot be used
// as some of the MST-based estimators use infinity to distinguish which vertices to test.
inline double getConnectionRadius(double dist_ij, double radius_i, double radius_j, double alpha, const int type) {

  // Only admit edges with finite weight if the neighborhood radii allow
  // Note that RSL will always form a complete hierarchy, so returning the numerical
  // limits maximum isn't necessary.
  double R;
  switch(type){
    // Robust Single Linkage from 2010 paper
    case 0:
      return std::max(dist_ij / alpha, std::max(radius_i, radius_j));
      break;
    // kNN graph from Algorithm 2 from Luxburgs 2014 paper
    case 1:
      R = alpha * std::max(radius_i, radius_j);
      return dist_ij <= R ? R : std::numeric_limits<double>::max();
      break;
    // mutual kNN graph from Algorithm 2 from Luxburgs 2014 paper
    case 2:
      R = alpha * std::min(radius_i, radius_j);
      return dist_ij <= R ? R : std::numeric_limits<double>::max();
      break;
    default:
      Rcpp::stop("Not a valid neighborhood query type");
    }
  return std::numeric_limits<double>::max();
}

// Recursively visit the merge matrix to extract an hclust sufficient ordering
void visit(const IntegerMatrix& merge, IntegerVector& order, int i, int j, int& ind) {
  // base case
  if (merge(i, j) < 0) {
    order.at(ind++) = -merge(i, j);
  }
  else {
    visit(merge, order, merge(i, j) - 1, 0, ind);
    visit(merge, order, merge(i, j) - 1, 1, ind);
  }
}

// Build a visually compatible ordering
IntegerVector extractOrder(IntegerMatrix merge){
  IntegerVector order = IntegerVector(merge.nrow()+1);
  int ind = 0;
  visit(merge, order, merge.nrow() - 1, 0, ind);
  visit(merge, order, merge.nrow() - 1, 1, ind);
  return(order);
}

// Regular Kruskal's MST
// [[Rcpp::export]]
NumericMatrix kruskalsMST(const NumericVector dist_x){
  std::string message = "kruskalsMST expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }

  // Set up resulting MST
  const int n = as<int>(dist_x.attr("Size"));
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get sorted radii
  NumericVector sorted_x = Rcpp::clone(dist_x).sort(false);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order");
  IntegerVector r_order = as<IntegerVector>(order(dist_x)) - 1;

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  int i = 0, crow = 0;
  for (NumericVector::const_iterator dist_ij = sorted_x.begin(); dist_ij != sorted_x.end(); ++dist_ij, ++i) {
    // Retrieve index of x_i and x_j
    int x_i = INDEX_TO(i, n), x_j = INDEX_FROM(i, n, x_i);
    if (components.Find(x_i) != components.Find(x_j)){
      mst(crow++, _) = NumericVector::create(x_i, x_j, *dist_ij);
      components.Union(x_i, x_j);
    }
  }
  return(mst);
}

/* primsMST
 * Compute MST using variant of Prim's based on a priority-first search for a dense graph
 */
// [[Rcpp::export]]
NumericMatrix primsMST(const NumericVector dist_x){
  std::string message = "primsMST expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }

  // Set up resulting MST
  const int n = as<int>(dist_x.attr("Size"));
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));

  double min, cedge_weight = 0.0;
  int c_i = 0, min_id = n - 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen

    // Compare all the new edge weights w/ the "current best" edge weights
    for (int t_i = 0; t_i < n; ++t_i) {
      if (t_i == c_i) continue;
      int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right

      // MST step, make sure node isn't already in the spanning tree
      if (v_selected[t_i] < 0) {

        // Updates edges on the fringe with lower weights if detected
        cedge_weight = dist_x[d_i];
        if (cedge_weight < fringe[t_i].weight) {
          fringe[t_i].weight = cedge_weight;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        // If edge 'on the fringe' is less than any of the current (head) nodes weight,
        // change the head-facing node to the edge of the fringe
        if (fringe[t_i].weight < min) {
          min = fringe[t_i].weight, min_id = t_i;
        }
      }
    }

    mst(n_edges, _) = NumericVector::create(min_id, c_i, min);
    v_selected[c_i] = 1;
    c_i = min_id;
  }
  return(mst);
}


/*
* Compute MST using variant of Prim's, constrained by the radius of the Balls around each x_i.
* Requires several array-type or indicator variables, namely:
* v_selected := array of indicators of spanning tree membership (-1 implies non-membership)
* c_i := index of current (head) node (relative to r_k)
* t_i := index of node to test against (relative to r_k)
* d_i := index of distance from current node to test node (relative to r)
*/
// [[Rcpp::export]]
NumericMatrix primsRSL(const NumericVector r, const NumericVector r_k, const int n, const double alpha, const int type){
  // Set up resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Data structures for prims
  std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
  std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));

  double min = std::numeric_limits<double>::infinity(), priority = 0.0;
  int c_i = 0, min_id = n - 1;
  for (int n_edges = 0; n_edges < n - 1; n_edges++) {
    if (n_edges % 1000 == 0) Rcpp::checkUserInterrupt();
    min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen

    // Compare all the new edge weights w/ the "current best" edge weights
    for (int t_i = 0; t_i < n; ++t_i) {
      // Graph subset step: ensure both c_i and t_i have neighborhood radii at least as big as the current radius
      if (t_i == c_i) continue;
      int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right

      // MST step, make sure node isn't already in the spanning tree
      if (v_selected[t_i] < 0) {

        // Generic RSL step
        priority = getConnectionRadius(r[d_i], r_k[c_i], r_k[t_i], alpha, type);
        if (priority < fringe[t_i].weight) { // using max above implicitly ensures c_i and t_i are connected
          fringe[t_i].weight = priority;
          fringe[t_i].to = c_i; // t_i indexes the 'from' node
        }

        // An edge 'on the fringe' might be less than any of the current nodes weights
        if (fringe[t_i].weight < min) {
          min = fringe[t_i].weight, min_id = t_i;
        }
      }
    }
    // Rcout << "Adding edge: (" << min_id << ", " << c_i << ") [" << min << "]" << std::endl;
    mst(n_edges, _) = NumericVector::create(min_id, c_i, type == 0 ? min * alpha : min);
    v_selected[c_i] = 1;
    c_i = min_id;
  }
  return(mst);
}


/* mstToHclust
 * Given a minimum spanning tree of the columnar form (<from>, <to>, <height>), create a valid hclust object
 * using a disjoint-set structure to track components
 * Notes: expects 0-based mst indices, and that all edge indices are 0 <= i < n
 */
// [[Rcpp::export(name = "mstToHclust")]]
List mstToHclust(NumericMatrix mst, const int n){

  // Check indices to make sure visiting the ordering is ok!
  NumericMatrix::Column from = mst.column(0);
  NumericMatrix::Column to = mst.column(1);
  if (all(from >= 0 & from <= n-1).is_false() || all(to >= 0 & to <= n-1).is_false()){
    Rcpp::stop("Invalid spanning tree indices passed.");
  }

  // Extract merge heights and associated order of such heights
  NumericVector dist = mst.column(2);
  IntegerVector height_order = order_(dist) - 1;

  // Set up components
  UnionFind components = UnionFind(n);
  IntegerVector comp_index = rep(0, n);
  std::vector<int> assigned = std::vector<int>(n, 0); // 0 to indicate node has not been assigned

  // Go through each row of MST, tracking component indices
  IntegerMatrix merge = IntegerMatrix(n - 1, 2);
  for (int i = 0; i < n - 1; ++i) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();
    NumericVector crow = mst.row(height_order.at(i));
    int from = crow.at(0), to = crow.at(1);
    int from_comp = components.Find(from), to_comp = components.Find(to);
    components.Union(from, to);

    // CASE 1: Two singletons a merging
    if (assigned.at(from) == 0 && assigned.at(to) == 0) { // Rcout << "1";
      merge.row(i) = IntegerVector::create(-(from + 1), -(to + 1));
      comp_index.at(components.Find(from)) = i;
    }
    // CASE 2: One singleton and one component are merging
    else if (assigned.at(from) == 0 || assigned.at(to) == 0) { // Rcout << "2";
      int leaf = assigned.at(from) == 0 ? from : to;
      int comp = assigned.at(from) == 0 ? to_comp : from_comp;
      merge.row(i) = IntegerVector::create(-(leaf+1), comp_index.at(comp)+1);
      comp_index.at(components.Find(leaf)) = i;
    }
    // CASE 3: Two components
    else {
      merge.row(i) = IntegerVector::create(comp_index.at(from_comp)+1, comp_index.at(to_comp)+1);
      comp_index.at(components.Find(from)) = i;
    }
    assigned.at(from) = assigned.at(to) = 1;
  }

  // Extractor merge order and return
  List res = List::create(
    _["merge"] = merge,
    _["height"] = dist[height_order],
    _["order"] = extractOrder(merge),
    _["labels"] = R_NilValue
    //_["mst"] = mst
  );
  res.attr("class") = "hclust";
  return(res);
}


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


// void normalizeIndices(std::unordered_map<int, cl_info>& cl, NumericMatrix& merge_heights){
//   NumericVector height = merge_heights.column(2);
//   NumericMatrix merge = merge_heights(_, Range(0,1));
//   IntegerVector merge_order = order_(height);
//
//   for (int i = 0; i < merge_heights.nrow(); ++i){
//     NumericMatrix::Row m_idx = merge.row(i);
//
//   }
//
// }

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
  int root_id = 0;
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

  // return(List::create(_["from"] = o_from, _["to"] = o_to, _["level_idx"] = level_idx, _["height"] = o_height,
  //                     _["comp_idx"] = comp_idx, _["cl_hierarchy"] = cl_hierarchy,
  //                     _["new_merge"] = new_merge, _["comp_idx2"] = comp_idx2));
  // // Normalize the indices
  // IntegerVector all_ids = Rcpp::unique(combine(from, to)); // original leaf ids
  // IntegerVector norm_from = match(from, all_ids) - 1; // 0-based
  // IntegerVector norm_to = match(to, all_ids) - 1; // 0-based
  //
  // // Use the mst to generate an hclust object
  // const int n_mst = all_ids.size();
  // NumericMatrix mst = no_init_matrix(n_mst - 1, 3);
  // mst(_, 0) = as<NumericVector>(norm_from);
  // mst(_, 1) = as<NumericVector>(norm_to);
  // mst(_, 2) = as<NumericVector>(height);
  // List simple_hclust = mstToHclust(mst, n_mst);


  // // Step 1 to point simplification: Record the singleton points in a new list
  // List point_idx = List(all_ids.size());
  // int i = 1;
  // for (IntegerVector::iterator id = all_ids.begin(); id != all_ids.end(); ++id, ++i){
  //   point_idx[patch::to_string(i)] = cl[*id].contains; // save the points corresponding to each cluster
  // }
  //
  // // Step 2 to point simplification: One more agglomerative pass...
  // // Need to map the new non-singleton normalized component indices to the points they contain
  // IntegerMatrix norm_merge = simple_hclust["merge"];
  // const int nrow_merge = norm_merge.nrow();
  // IntegerMatrix orig_merge = IntegerMatrix(nrow_merge, 2);
  // std::vector<int> comp_map = std::vector<int>(nrow_merge); // map from merge step to the original parent idx
  // for (int row_i = 0; row_i < nrow_merge; ++row_i){
  //   IntegerVector pt_ids = norm_merge(row_i, _);
  //   if (all(pt_ids < 0).is_true()){
  //     int orig_left = all_ids.at((-pt_ids.at(0))-1), orig_right = all_ids.at((-pt_ids.at(1))-1);
  //     Rcout << "Case 0: (" << orig_left << ", " << orig_right << ")\n";
  //     orig_merge(row_i, _) = IntegerVector::create(-orig_left, -orig_right);
  //     std::pair<int, int> key = std::minmax(orig_left, orig_right);
  //     if (parent_map.find(key) != parent_map.end()){
  //       int orig_parent_idx = parent_map.at(key);
  //       comp_map[row_i] = orig_parent_idx;
  //     }
  //   } else if (any(pt_ids < 0).is_true()){
  //     int orig_leaf = pt_ids.at(0) < 0 ? all_ids.at((-pt_ids.at(0))-1) : all_ids.at((-pt_ids.at(1))-1);
  //     int orig_comp = pt_ids.at(0) < 0 ? comp_map.at(pt_ids.at(1) - 1) : comp_map.at(pt_ids.at(0) - 1);
  //     Rcout << "Case 1: (" << orig_leaf << ", " << orig_comp << ")\n";
  //     orig_merge(row_i, _) = IntegerVector::create(-orig_leaf, orig_comp);
  //     std::pair<int, int> key = std::minmax(orig_leaf, orig_comp);
  //     if (parent_map.find(key) != parent_map.end()){
  //       int orig_parent_idx = parent_map.at(key);
  //       comp_map[row_i] = orig_parent_idx;
  //     }
  //   } else {
  //     int orig_left = comp_map.at(pt_ids.at(0) - 1);
  //     int orig_right = comp_map.at(pt_ids.at(1) - 1);
  //     orig_merge(row_i, _) = IntegerVector::create(orig_left, orig_right);
  //     Rcout << "Case 2: (" << orig_left << ", " << orig_right << ")\n";
  //     std::pair<int, int> key = std::minmax(orig_left, orig_right);
  //     if (parent_map.find(key) != parent_map.end()){
  //       int orig_parent_idx = parent_map.at(key);
  //       comp_map[row_i] = orig_parent_idx;
  //     }
  //   }
  // }
  //
  // for (std::map<std::pair<int, int>, int>::iterator it = parent_map.begin(); it != parent_map.end(); ++it){
  //   Rcout << "(" << it->first.first << ", " << it->first.second << ") ==> " << it->second << std::endl;
  // }
  //

  // Add the original point indices to the hclust object, and prepend "simplified_hclust" to
  // the S3 class name
  // simple_hclust["leaf_idx"] = point_idx;
  // std::vector<std::string> class_names = std::vector<std::string>(2);
  // class_names[0] = "simplified_hclust"; // Prefer simplified S3 methods
  // class_names[1] = "hclust"; // It still IS-A valid hclust object
  // simple_hclust.attr("class") = wrap(class_names);
  // return List::create(simple_hclust, _["parents"] = wrap(comp_map), _["leaves"] = all_ids, _["orig_merge"] = orig_merge,
  //                     _["cl_hierarchy"] = cl_hierarchy,
  //                     _["mst"] = mst);
}

// leaves <- unname(unlist(sapply(sapply(sapply(ls(what$cl_hierarchy), function(key) { cand <- what$cl_hierarchy[key]; Filter(f = function(x) ! x %in% ls(what$cl_hierarchy), cand)}), unlist), unname)))

// [[Rcpp::export]]
List clusterTree(const NumericVector dist_x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095,
                 const int type = 0) {
  std::string message = "clusterTree expects a 'dist' object.";
  if (!dist_x.hasAttribute("class") || as<std::string>(dist_x.attr("class")) != "dist") { stop(message); }
  if (!dist_x.hasAttribute("method")) { stop(message); }
  if (!dist_x.hasAttribute("Size")){ stop(message); }
  if (as<std::string>(dist_x.attr("method")) != "euclidean") { warning("RSL expects euclidean distances.");}

  // Number of data points
  const int n = as<int>(dist_x.attr("Size"));

  // Get ordered radii, use R order function to get consistent ordering
  NumericVector r = Rcpp::clone(dist_x);

  // Run the MST with set parameters
  NumericMatrix mst = primsRSL(r, r_k, n, alpha, type);

  // Convert to HCLUST object
  List res = mstToHclust(mst, n);
  return (res);
}


