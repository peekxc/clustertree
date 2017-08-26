#include "greedy_prims.h"

int findNearest(ANNpoint q, ANNkd_tree* qtree, double eps){

  int ANNkdDim = qtree->dim;				// dimension of space
  double ANNkdMaxErr = ANN_POW(1.0 + eps);;			// max tolerable squared error
  ANNpointArray	pts = qtree->pts;				// the points

  std::stack<ANNdist> box_dists = std::stack<ANNdist>(); // manage the stack manually!
  std::stack<ANNcoord> cut_diffs = std::stack<ANNcoord>(); // manage the stack manually!
  std::stack<ANNkd_node*> nodes = std::stack<ANNkd_node*>(); // manage the stack manually!

  // Initial prep
  ANNdist box_dist = annBoxDistance(q, qtree->bnd_box_lo, qtree->bnd_box_hi, qtree->dim);
  box_dists.push(box_dist); // first box_dist
  nodes.push(qtree->root); // first node to search (root)
  ANNmin_k knn_queue = ANNmin_k(1);	// create set to update for k = 1 nearest neighbors

  // Use hash map to track which nodes have been visited
  std::unordered_map<ANNkd_node*, bool> visited = std::unordered_map<ANNkd_node*, bool>();

  // Start the iteration!
  while(!nodes.empty()){
    ANNkd_node* cnode = nodes.top();
    ANNdist cbox_dist = box_dists.top();

    // "Base case"
    if(IS_LEAF(cnode)){
      ANNkd_leaf* cleaf = AS_LEAF(cnode);
      ANNdist dist;				// distance to data point
      ANNcoord* pp;				// data coordinate pointer
      ANNcoord* qq;				// query coordinate pointer
      ANNdist min_dist;			// distance to current closest point
      ANNcoord t;
      int d;

      min_dist = knn_queue.max_key(); // smallest distance so far
      for (int i = 0; i < n; i++) {	// check points in bucket
        pp = pts[cleaf->bkt[i]];			// first coord of next data point
        qq = q;					// first coord of query point
        dist = 0;

        for(d = 0; d < ANNkdDim; d++) {
          t = *(qq++) - *(pp++);		// compute length and adv coordinate
          // exceeds dist to k-th smallest?
          if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
            break;
          }
        }

        if (d >= ANNkdDim && (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
          // add it to the list
          knn_queue.insert(dist, cleaf->bkt[i]);
          min_dist = knn_queue.max_key();
        }
      }

      // Go back up
      nodes.pop();
      cbox_dist.pop();
    }
    else if(IS_SPLIT(cnode)){
      ANNkd_split* cspl_node = AS_SPLIT(cnode);
      ANNcoord cut_diff = q[cspl_node->cut_dim] - cspl_node->cut_val; // distance to cutting plane
      if (cut_diff < 0) { // left of cutting plane
        if (!visited[cnode]){
          visited[cnode] = true; // mark this split node as visited
          nodes.push(cspl_node->child[ANN_LO]);
          box_dists.push(cbox_dist);
          continue; // "recurse"
        } else {
          ANNcoord box_diff = cspl_node->cd_bnds[ANN_LO] - q[cspl_node->cut_dim];
          if (box_diff < 0)				// within bounds - ignore
            box_diff = 0;
          // distance to further box
          box_dist = (ANNdist) ANN_SUM(box_dist,
                      ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));
          if (box_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push(cspl_node->child[ANN_HI]);
            box_dists.push(box_dist);
            continue;
          }
        }
      }
      else { // right of cutting plane
        if (!visited[cnode]){
          visited[cnode] = true; // mark this split node as visited
          nodes.push(cspl_node->child[ANN_HI]); // visit closer child first
          box_dists.push(cbox_dist);
          continue; // "recurse"
        } else {
          ANNcoord box_diff = q[cspl_node->cut_dim] - cspl_node->cd_bnds[ANN_HI];
          if (box_diff < 0)				// within bounds - ignore
            box_diff = 0;
          // distance to further box
          box_dist = (ANNdist) ANN_SUM(box_dist,
                      ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

          // visit further child if close enough
          if (box_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push(cspl_node->child[ANN_LO]);
            box_dists.push(box_dist);
            continue;
          }
        }
      }
    } // else if(IS_SPLIT(cnode)){
  } // while(!nodes.empty())
}

// void ann_greedy_prims(const NumericMatrix& query_pts,
//                        ANNkd_tree* qtree,
//                        double				eps)			// the error bound
// {
//   const int n = query_pts.nrow();
//
//   // Set up resulting MST
//   NumericMatrix mst = NumericMatrix(n - 1, 3);
//
//   // Data structures for prims
//   std::vector<int> v_selected = std::vector<int>(n, -1); // -1 to indicate node is not in MST
//   std::vector<edge> fringe = std::vector<edge>(n, edge(-1, std::numeric_limits<double>::infinity()));
//
//   int				ANNkdDim = qtree->dim;				// dimension of space
//   ANNpoint		ANNkdQ;					// query point
//   double			ANNkdMaxErr;			// max tolerable squared error
//   ANNpointArray	pts = qtree->pts;				// the points
//   ANNmin_k		*ANNkdPointMK;			// set of k closest points
//
//   double min = std::numeric_limits<double>::infinity(), priority = 0.0;
//   int c_i = 0, min_id = n - 1;
//
//   ANNpoint q = qtree->pts[0]; // first query point
//   ANNdist box_dist = annBoxDistance(q, qtree->bnd_box_lo, qtree->bnd_box_hi, qtree->dim);
//
//   std::stack<ANNdist> box_dists = std::stack<ANNdist>(); // manage the stack manually!
//   std::stack<ANNcoord> cut_diffs = std::stack<ANNcoord>(); // manage the stack manually!
//   std::stack<ANNkd_node*> nodes = std::stack<ANNkd_node*>();
//
//   box_dists.push(box_dist); // first box_dist
//   nodes.push(qtree->root);
//
//   ANNmin_k knn_queue = ANNmin_k(1);	// create set to update for k = 1 nearest neighbors
//
//   std::unordered_map<ANNkd_node*, bool> visited = std::unordered_map<ANNkd_node*, bool>();
//
//   for (int n_edges = 0; n_edges < n - 1; n_edges++) {
//     if (n_edges % 100 == 0) Rcpp::checkUserInterrupt();
//     min = std::numeric_limits<double>::infinity(); // Reset so new edge is always chosen
//
//     ANNkdMaxErr = ANN_POW(1.0 + eps);
//     // The search
//
//
//       // Compare all the new edge weights w/ the "current best" edge weights
//       for (int t_i = 0; t_i < n; ++t_i) {
//         // Graph subset step: ensure both c_i and t_i have neighborhood radii at least as big as the current radius
//         if (t_i == c_i) continue;
//         int d_i = t_i > c_i ? INDEX_TF(n, c_i, t_i) : INDEX_TF(n, t_i, c_i); // bigger index always on the right
//
//         // MST step, make sure node isn't already in the spanning tree
//         if (v_selected[t_i] < 0) {
//
//           // Generic RSL step
//           priority = getConnectionRadius(r[d_i], r_k[c_i], r_k[t_i], alpha, type);
//           if (priority < fringe[t_i].weight) { // using max above implicitly ensures c_i and t_i are connected
//             // Rcout << "Updating fringe " << t_i << "w/ radii (" << r_k[c_i] << ", " << r_k[t_i] << ")" << std::endl;
//             // Rcout << "current: " << c_i << ", to: " << t_i << std::endl;
//             fringe[t_i].weight = priority;
//             fringe[t_i].to = c_i; // t_i indexes the 'from' node
//           }
//
//           // An edge 'on the fringe' might be less than any of the current nodes weights
//           if (fringe[t_i].weight < min) {
//             min = fringe[t_i].weight, min_id = t_i;
//           }
//         }
//       }
//       // Rcout << "Adding edge: (" << min_id << ", " << c_i << ") [" << min << "]" << std::endl;
//       mst(n_edges, _) = NumericVector::create(min_id, c_i, min);
//       v_selected[c_i] = 1;
//       c_i = min_id;
//     }
//     return(mst);
//   }
