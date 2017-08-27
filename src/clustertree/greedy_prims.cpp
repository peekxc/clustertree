#include "greedy_prims.h"

simple_edge findNearest(const int q_idx, ANNkd_tree* qtree, std::vector<int> ignore, double eps){

  //const int n = qtree->nPoints();
  int ANNkdDim = qtree->dim;				// dimension of space
  double ANNkdMaxErr = ANN_POW(1.0 + eps);;			// max tolerable squared error
  ANNpointArray	pts = qtree->pts;				// the points
  ANNpoint q = pts[q_idx]; // query pt;

  std::stack<ANNdist> box_dists = std::stack<ANNdist>(); // manage the stack manually!
  std::stack<ANNcoord> cut_diffs = std::stack<ANNcoord>(); // manage the stack manually!
  std::stack<ANNkd_node*> nodes = std::stack<ANNkd_node*>(); // manage the stack manually!

  // Initial prep
  ANNdist box_dist = annBoxDistance(q, qtree->bnd_box_lo, qtree->bnd_box_hi, qtree->dim);
  ANNmin_k knn_queue = ANNmin_k(2);	// create set to update for k = 2 nearest neighbors

  // Use hash map to track which nodes have been visited
  std::unordered_map<ANNkd_node*, bool> visited = std::unordered_map<ANNkd_node*, bool>();

  // Start with the root
  box_dists.push(box_dist); // first box_dist
  nodes.push(qtree->root); // first node to search (root)

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
      for (int i = 0; i < cleaf->n_pts; i++) {	// check points in bucket
        pp = pts[cleaf->bkt[i]];			// first coord of next data point
        qq = q;					// first coord of query point
        dist = 0;

        // If the point is initialized in the ignore vector, don't calculate distance
        if (ignore[cleaf->bkt[i]] >= 0 || cleaf->bkt[i] == q_idx){ // don't allow self match
          dist = std::numeric_limits<ANNdist>::max();
          d = 0; // mark the distance as invalid
        } else {
          for(d = 0; d < ANNkdDim; d++) {
            t = *(qq++) - *(pp++);		// compute length and adv coordinate
            // exceeds dist to k-th smallest?
            if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
              break;
            }
          }
        }

        if (d >= ANNkdDim) { // allow 0-dist to allow duplicate points, just not the same ids
          // add it to the list
          knn_queue.insert(dist, cleaf->bkt[i]);
          min_dist = knn_queue.max_key();
        }
      }

      // Go back up
      nodes.pop();
      box_dists.pop();
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
          cbox_dist = (ANNdist) ANN_SUM(cbox_dist,
                      ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

          // This split node has now been visited, and all calculation needed have been performed
          // So pop it off the stack and evaluate whether to keep going
          nodes.pop();
          box_dists.pop();

          // Check to to see if should go right
          if (cbox_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push(cspl_node->child[ANN_HI]);
            box_dists.push(cbox_dist);
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
          cbox_dist = (ANNdist) ANN_SUM(cbox_dist,
                      ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

          // This split node has now been visited, and all calculation needed have been performed
          // So pop it off the stack and evaluate whether to keep going
          nodes.pop();
          box_dists.pop();

          // visit further child if close enough
          if (cbox_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push(cspl_node->child[ANN_LO]);
            box_dists.push(cbox_dist);
            continue;
          }
        }
      }
    } // else if(IS_SPLIT(cnode)){
  } // while(!nodes.empty())

  simple_edge shortest_edge = simple_edge();
  //int same_id = knn_queue.ith_smallest_info(0) == q_idx;
  shortest_edge.to = knn_queue.ith_smallest_info(0); // Should never be the same id
  shortest_edge.weight = knn_queue.ith_smallest_key(0); // Should never be the same id
  return shortest_edge;
}

void knn_stack(const int q_idx, // the query point (index)
               const int k, 		// number of near neighbors to return
               ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
               ANNdistArray		dd,				// the approximate nearest neighbor
               ANNkd_tree* qtree, // The kd tree used to make the points
               double eps){

  const int n = qtree->nPoints();
  int ANNkdDim = qtree->dim;				// dimension of space
  double ANNkdMaxErr = ANN_POW(1.0 + eps);;			// max tolerable squared error
  ANNpointArray	pts = qtree->pts;				// the points
  ANNpoint q = pts[q_idx]; // query pt;

  std::vector<ANNdist> box_dists = std::vector<ANNdist>();
  std::vector<ANNcoord> cut_diffs = std::vector<ANNcoord>();
  std::vector<ANNkd_node*> nodes = std::vector<ANNkd_node*>();
  const int expected_depth = floor(log2(n));
  box_dists.reserve(expected_depth);
  cut_diffs.reserve(expected_depth);
  nodes.reserve(expected_depth);
  // std::stack<ANNdist> box_dists = std::stack<ANNdist>(); // manage the stack manually!
  // std::stack<ANNcoord> cut_diffs = std::stack<ANNcoord>(); // manage the stack manually!
  // std::stack<ANNkd_node*> nodes = std::stack<ANNkd_node*>(); // manage the stack manually!

  // Initial prep
  ANNdist box_dist = annBoxDistance(q, qtree->bnd_box_lo, qtree->bnd_box_hi, qtree->dim);
  ANNmin_k knn_queue = ANNmin_k(k);	// create set to update for k = 2 nearest neighbors

  // Use hash map to track which nodes have been visited
  std::unordered_map<ANNkd_node*, bool> visited = std::unordered_map<ANNkd_node*, bool>();

  // Start with the root
  box_dists.push_back(box_dist); // first box_dist
  nodes.push_back(qtree->root); // first node to search (root)

  // Start the iteration!
  while(!nodes.empty()){
    ANNkd_node* cnode = nodes.back();
    ANNdist cbox_dist = box_dists.back();
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
      for (int i = 0; i < cleaf->n_pts; i++) {	// check points in bucket
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

        if (d >= ANNkdDim) { // allow 0-dist to allow duplicate points, just not the same ids
          // add it to the list
          knn_queue.insert(dist, cleaf->bkt[i]);
          min_dist = knn_queue.max_key();
        }
      }

      // Go back up
      nodes.pop_back();
      box_dists.pop_back();
    }
    else if(IS_SPLIT(cnode)){
      const ANNkd_split& cspl_node = *AS_SPLIT(cnode);
      ANNcoord cut_diff = q[cspl_node.cut_dim] - cspl_node.cut_val; // distance to cutting plane
      if (cut_diff < 0) { // left of cutting plane
        if (!visited[cnode]){
          visited[cnode] = true; // mark this split node as visited
          nodes.push_back(cspl_node.child[ANN_LO]);
          box_dists.push_back(cbox_dist);
          continue; // "recurse"
        } else {
          ANNcoord box_diff = cspl_node.cd_bnds[ANN_LO] - q[cspl_node.cut_dim];
          if (box_diff < 0)				// within bounds - ignore
            box_diff = 0;
          // distance to further box
          cbox_dist = (ANNdist) ANN_SUM(cbox_dist,
                       ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

          // This split node has now been visited, and all calculation needed have been performed
          // So pop it off the stack and evaluate whether to keep going
          nodes.pop_back();
          box_dists.pop_back();

          // Check to to see if should go right
          if (cbox_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push_back(cspl_node.child[ANN_HI]);
            box_dists.push_back(cbox_dist);
            continue;
          }
        }
      }
      else { // right of cutting plane
        if (!visited[cnode]){
          visited[cnode] = true; // mark this split node as visited
          nodes.push_back(cspl_node.child[ANN_HI]); // visit closer child first
          box_dists.push_back(cbox_dist);
          continue; // "recurse"
        } else {
          ANNcoord box_diff = q[cspl_node.cut_dim] - cspl_node.cd_bnds[ANN_HI];
          if (box_diff < 0)				// within bounds - ignore
            box_diff = 0;
          // distance to further box
          cbox_dist = (ANNdist) ANN_SUM(cbox_dist,
                       ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

          // This split node has now been visited, and all calculation needed have been performed
          // So pop it off the stack and evaluate whether to keep going
          nodes.pop_back();
          box_dists.pop_back();

          // visit further child if close enough
          if (cbox_dist * ANNkdMaxErr < knn_queue.max_key()){
            nodes.push_back(cspl_node.child[ANN_LO]);
            box_dists.push_back(cbox_dist);
            continue;
          }
        }
      }
    } // else if(IS_SPLIT(cnode)){
  } // while(!nodes.empty())

  for (int i = 0; i < k; i++) {		// extract the k-th closest points
    dd[i] = knn_queue.ith_smallest_key(i);
    nn_idx[i] = knn_queue.ith_smallest_info(i);
  }
  return;
}
