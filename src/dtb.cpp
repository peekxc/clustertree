#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>     // std::cout, std::ostream, std::ios
#include <fstream>      // std::filebuf

// ANN library
#include "ANN/ANN.h"
#include "kd_search.h" // kd-search declarations
#include "bd_tree.h"	 // bd-tree declarations

// KD tree R-accessible functions
#include "R_kdtree.h"

// ANN extensions
#include "DT_Abstract.h"

// Static globals to simplify recursion
int				    kd_dim;				    // dimension of space
ANNpoint		  kd_query_pt;			// query point
double			  kd_max_sqerror;		// max tolerable squared error
ANNpointArray	kd_pts;				    // the points
ANNmin_k		  *k_min_pts;			  // set of k closest points

// This function assumes the caller is the 'reference' tree, with the query tree
// attached via pointer. This sets up the dual tree traversal on the augmented split and leaf nodes
void ANNkd_tree::ann_dt_kSearch(
    ANNpoint		  q,				    // query point
    int				    k,				    // number of near neighbors to return
    ANNidxArray	  nn_idx,			  // nearest neighbor array (modified)
    ANNdistArray  dd,				    // dist to near neighbors (modified)
    ANNkd_tree* query_tree,  // query tree for dual tree traversal
    double			  eps           // error bound
){
  kd_dim = dim;						// copy arguments to static equivs
  kd_query_pt = q;
  kd_pts = pts;
  ANNptsVisited = 0;					// initialize count of points visited

  if (k > n_pts) {					// too many near neighbors?
    annError("Requesting more near neighbors than data points", ANNabort);
  }

  kd_max_sqerror = ANN_POW(1.0 + eps);
  ANN_FLOP(2)							// increment floating op count

  k_min_pts = new ANNmin_k(k);		// create set for closest k points

  // search starting at the root
  //query_tree->root->getStats()
  // pidx[0];
  root->ann_dt_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));

  for (int i = 0; i < k; i++) {		// extract the k-th closest points
    dd[i] = k_min_pts->ith_smallest_key(i);
    nn_idx[i] = k_min_pts->ith_smallest_info(i);
  }
  delete k_min_pts;				// deallocate closest point set
}

//----------------------------------------------------------------------
//	kd_split::ann_dt_search - search a splitting node (dual tree version)
//----------------------------------------------------------------------

void ANNkd_split::ann_dt_search(ANNdist box_dist)		// dual tree search
{
  // check dist calc term condition
  if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

  // First, it stores the integer cutting dimension (from 0 to d - 1) indicating the coordinate axis that
  // is orthogonal to the cutting hyperplane.
  // this->cut_dim;

  // Second it stores the cutting value where this plane intersects this axis.
  // this->cut_val;

  // It also contains the two pointers to the left and right children
  // (containing points lying to the low and high side of the cutting plane, respectively).
  // this->child;

  // In order to implement incremental distance computation, it also stores two associated pieces of
  // information with the cell. Consider the two hyperplanes bounding the cell that are parallel to the
  // cutting plane. It stores the values at which these two hyperplanes intersect the cutting dimension's axis.
  // this->cd_bnds;

  // distance to cutting plane
  ANNcoord cut_diff = kd_query_pt[cut_dim] - cut_val; // kd_query_pt := static copy of query point (ANNpoint)

  if (cut_diff < 0) {					// left of cutting plane
    child[ANN_LO]->ann_dt_search(box_dist);// visit closer child first

    ANNcoord box_diff = cd_bnds[ANN_LO] - kd_query_pt[cut_dim];
    if (box_diff < 0)				// within bounds - ignore
      box_diff = 0;
    // distance to further box
    box_dist = (ANNdist) ANN_SUM(box_dist,
                ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

    // visit further child if close enough
    if (box_dist * kd_max_sqerror < k_min_pts->max_key())
      child[ANN_HI]->ann_dt_search(box_dist);

  }
  else {								// right of cutting plane
    child[ANN_HI]->ann_dt_search(box_dist);// visit closer child first

    ANNcoord box_diff = kd_query_pt[cut_dim] - cd_bnds[ANN_HI];
    if (box_diff < 0)				// within bounds - ignore
      box_diff = 0;
    // distance to further box
    box_dist = (ANNdist) ANN_SUM(box_dist,
                ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

    // visit further child if close enough
    if (box_dist * kd_max_sqerror < k_min_pts->max_key())
      child[ANN_LO]->ann_dt_search(box_dist);

  }
  ANN_FLOP(10)		// increment floating ops
  ANN_SPL(1)			// one more splitting node visited
}


void ANNkd_leaf::ann_dt_search(ANNdist box_dist)
{
  ANNdist dist;				// distance to data point
  ANNcoord* pp;				// data coordinate pointer
  ANNcoord* qq;				// query coordinate pointer
  ANNdist min_dist;			// distance to k-th closest point
  ANNcoord t;
  int d;

  min_dist = k_min_pts->max_key(); // k-th smallest distance so far

  for (int i = 0; i < n_pts; i++) {	// check points in bucket

    pp = kd_pts[bkt[i]];			// first coord of next data point
    qq = kd_query_pt;					// first coord of query point
    dist = 0;

    for(d = 0; d < kd_dim; d++) {
      ANN_COORD(1)				// one more coordinate hit
      ANN_FLOP(4)					// increment floating ops

      t = *(qq++) - *(pp++);		// compute length and adv coordinate
      // exceeds dist to k-th smallest?
      if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
        break;
      }
    }

    if (d >= kd_dim &&					// among the k best?
        (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
      // add it to the list
      k_min_pts->insert(dist, bkt[i]);
      min_dist = k_min_pts->max_key();
    }
  }
  ANN_LEAF(1)							  // one more leaf node visited
  ANN_PTS(n_pts)						// increment points visited
  ANNptsVisited += n_pts;		// increment number of points visited
}
// double knnBaseCase(ANNkd_split N_q, ANNkd_split N_r){
//   N_q.
// }
//

// double Score(ANNcoord* N_q, ANNcoord* N_r){
//   N_q
// }

// // Base Depth-First Search framework for the Dual Tree traversals
// void DFS(DT_Abstract* qtree, DT_Abstract* rtree){
//   // if (Score){
//   //
//   // }
// }


// [[Rcpp::export]]
List DTB(NumericMatrix x, int k = 5) {

  // Copy data over to ANN point array
  int nrow = x.nrow(), ncol = x.ncol();

  ANNpointArray dataPts = annAllocPts(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      (dataPts[i])[j] = x(i, j);
    }
  }

  // Create kd tree
  ANNkd_tree* kdTree = new ANNkd_tree(dataPts, nrow, ncol, 30, ANN_KD_SUGGEST);

  // Note: the search also returns the point itself (as the first hit)!
  // So we have to look for k+1 points.
  ANNdistArray dists = new ANNdist[k+1];
  ANNidxArray nnIdx = new ANNidx[k+1];

  // Distance matrix of kNN distances
  NumericMatrix d(x.nrow(), k);

  // Id matrix of knn indices
  IntegerMatrix id(x.nrow(), k);

  // Loop through the query points and perform the search
  for (int i=0; i < x.nrow(); i++) {
    if (!(i % 100)) Rcpp::checkUserInterrupt();

    // Retrieve the query point
    NumericVector query_pt = x.row(i);
    ANNpoint queryPt = &query_pt.at(0);

    // Search the kd tree
    kdTree->ann_dt_kSearch(queryPt, k+1, nnIdx, dists, kdTree,  0);

    // Remove self matches
    IntegerVector ids = IntegerVector(nnIdx, nnIdx+k+1);
    LogicalVector take = ids != i;
    ids = ids[take];
    id(i, _) = ids + 1;

    // Convert to regular, non-squared euclidean distances
    NumericVector ndists = NumericVector(dists, dists+k+1)[take];
    d(i, _) = sqrt(ndists);
  }
  List res = List::create(_["dist"] = d, _["id"] = id);

  DT_Abstract abstract_dt = DT_Abstract(kdTree);
  IntegerVector child_ids = abstract_dt.child_nodes();
  res["child_ids"] = child_ids;
  res["ids"] = abstract_dt.getIDXArray();
  res["box_dists"] = abstract_dt.convex_subset();
  // std::filebuf fb;
  // fb.open ("test_kdtree.txt",std::ios::out);
  // std::ostream os(&fb);
  // kdTreeR->Dump(ANNtrue, os);
  // fb.close();

  return(res);
}

/*** R
  test_set <- matrix(c(13, 291,
                       57, 145,
                       115, 232,
                       86, 27,
                       145, 28,
                       262, 203,
                       320, 261,
                       379, 174,
                       261, 71,
                       325, 57), byrow=T, ncol=2)
  test_set[, 2] <- 321 - test_set[, 2]
  clustertree:::DTB(test_set, k = 3L)
*/
