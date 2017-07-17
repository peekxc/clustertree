#include "ANNdt.h"
#include "utilities.h"

ANNkd_tree_dt::ANNkd_tree_dt(					// construct from point array
  ANNpointArray		pa,				// point array (with at least n pts)
  int					n,				// number of points
  int					dd,				// dimension
  int					bs,				// bucket size
  ANNsplitRule		split, // splitting method
  std::unordered_map<ANNkd_node*, const Bound& >* bounds) // bounds
{
  R_PRINTF("Dimension of trees: %d\n", dd);
  SkeletonTree(n, dd, bs);			// set up the basic stuff
  pts = pa;							// where the points are
  if (n == 0) return;					// no points--no sweat

  ANNorthRect bnd_box(dd);			// bounding box for points
  annEnclRect(pa, pidx, n, dd, bnd_box);// construct bounding rectangle
  R_PRINTF("Root box: [(%f, %f), (%f, %f)]\n", bnd_box.lo[0], bnd_box.hi[0], bnd_box.lo[1], bnd_box.hi[1]);

  // Compute the centroid of the current rectangle
  ANNpoint centroid = new ANNcoord[dd];
  for (int i = 0; i < dd; ++i){ centroid[i] = ANNcoord((bnd_box.lo[i] + bnd_box.hi[i]) / 2.0); }
  R_PRINTF("Centroid: (%f, %f)\n", centroid[0], centroid[1]);

  // copy to tree structure
  bnd_box_lo = annCopyPt(dd, bnd_box.lo);
  bnd_box_hi = annCopyPt(dd, bnd_box.hi);

  switch (split) {					// build by rule
  case ANN_KD_STD:					// standard kd-splitting rule
    root = rkd_tree_pr(pa, pidx, n, dd, bs, bnd_box, kd_split, bounds, centroid);
    break;
  case ANN_KD_MIDPT:					// midpoint split
    root = rkd_tree_pr(pa, pidx, n, dd, bs, bnd_box, midpt_split, bounds, centroid);
    break;
  case ANN_KD_FAIR:					// fair split
    root = rkd_tree_pr(pa, pidx, n, dd, bs, bnd_box, fair_split, bounds, centroid);
    break;
  case ANN_KD_SUGGEST:				// best (in our opinion)
  case ANN_KD_SL_MIDPT:				// sliding midpoint split
    root = rkd_tree_pr(pa, pidx, n, dd, bs, bnd_box, sl_midpt_split, bounds, centroid);
    break;
  case ANN_KD_SL_FAIR:				// sliding fair split
    root = rkd_tree_pr(pa, pidx, n, dd, bs, bnd_box, sl_fair_split, bounds, centroid);
    break;
  default:
    annError("Illegal splitting method", ANNabort);
  }
  R_INFO("Root created: " << root << "\n")

  Bound& bnd = *new Bound();
  bnd.bnd_box = new ANNorthRect(dim, bnd_box.lo, bnd_box.hi); // Copy bounding box
  bnd.centroid = annCopyPt(dim, centroid);
  ANNdist tmp = ANN_DIST_INF;
  for (int i = 0; i < n; ++i){
    ANNdist centroid_dist = annDist(dim, (ANNpoint) centroid, (ANNpoint) pa[pidx[i]]);
    tmp = centroid_dist < tmp ? centroid_dist : tmp;
  };

  // Update max child and max desc. distance with upper bounds
  bnd.lambda = annDist(dim, bnd.bnd_box->hi, bnd.centroid); // upper bound
  bnd.rho = bnd.lambda; // also upper bound

  // Insert root as final key
  bounds->insert(std::pair< ANNkd_ptr, const Bound& >(root, bnd));

}

// ANNkd_tree_dt::ANNkd_tree_dt

// Recursive procedure to actually create the tree
ANNkd_ptr rkd_tree_pr(				// recursive construction of kd-tree
    ANNpointArray		pa,				// point array
    ANNidxArray			pidx,			// point indices to store in subtree
    int					n,				// number of points
    int					dim,			// dimension of space
    int					bsp,			// bucket space
    ANNorthRect			&bnd_box,		// bounding box for current node
    ANNkd_splitter		splitter, // splitting routine
    std::unordered_map<ANNkd_node*, const Bound& >* bounds, // bounds map
    ANNpoint centroid) // centroid of rectangle for current node
{
  if (n <= bsp) {						// n small, make a leaf node
    if (n == 0)						// empty leaf node
      return KD_TRIVIAL;			// return (canonical) empty leaf
    else	{    // construct the node and return
      ANNkd_node* new_leaf = new ANNkd_leaf(n, pidx);
      Bound& leaf_bnd = *new Bound();
      leaf_bnd.bnd_box = new ANNorthRect(dim, bnd_box.lo, bnd_box.hi); // Copy bounding box
      for (int i = 0; i < dim; ++i){ centroid[i] = ANNcoord((bnd_box.lo[i] + bnd_box.hi[i]) / 2.0); }
      leaf_bnd.centroid = annCopyPt(dim, centroid);
      // leaf_bnd.centroid = annCopyPt(dim, centroid);
      ANNdist tmp = ANN_DIST_INF;
      for (int i = 0; i < n; ++i){
        ANNdist centroid_dist = annDist(dim, (ANNpoint) centroid, (ANNpoint) pa[pidx[i]]);
        tmp = centroid_dist < tmp ? centroid_dist : tmp;
      };
      leaf_bnd.lambda = tmp; // annDist(dim, leaf_bnd.bnd_box->hi, leaf_bnd.centroid); // upper bound
      leaf_bnd.rho = tmp; // maximum child distance
      bounds->insert(std::pair< ANNkd_ptr, const Bound& >(new_leaf, leaf_bnd));
      return new_leaf;
    }
  }
  else {								// n large, make a splitting node
    int cd;							// cutting dimension
    ANNcoord cv;					// cutting value
    int n_lo;						// number on low side of cut
    ANNkd_node *lo, *hi;			// low and high children

    // invoke splitting procedure
    (*splitter)(pa, pidx, bnd_box, n, dim, cd, cv, n_lo);

    // Create the bounds object for this node
    Bound& node_bnds = *new Bound();

    // Copy bounding box
    node_bnds.bnd_box = new ANNorthRect(dim, bnd_box.lo, bnd_box.hi); // Copy bounding box
    for (int i = 0; i < dim; ++i){ centroid[i] = ANNcoord((bnd_box.lo[i] + bnd_box.hi[i]) / 2.0); }
    node_bnds.centroid = annCopyPt(dim, centroid);

    ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
    ANNcoord hv = bnd_box.hi[cd];

    bnd_box.hi[cd] = cv;			// modify bounds for left subtree
    lo = rkd_tree_pr(					// build left subtree
      pa, pidx, n_lo,			// ...from pidx[0..n_lo-1]
      dim, bsp, bnd_box, splitter,
      bounds, centroid); // propagate bounds map and centroid
    bnd_box.hi[cd] = hv;			// restore bounds

    bnd_box.lo[cd] = cv;			// modify bounds for right subtree
    hi = rkd_tree_pr(					// build right subtree
      pa, pidx + n_lo, n-n_lo,// ...from pidx[n_lo..n-1]
      dim, bsp, bnd_box, splitter,
      bounds, centroid); // propagate bounds map and centroid
    bnd_box.lo[cd] = lv;			// restore bounds

    // create the splitting node
    ANNkd_split *ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);

    // Update max child and max desc. distance with upper bounds
    ANNdist hi_dist = annDist(dim, node_bnds.bnd_box->hi, node_bnds.centroid);
    ANNdist lo_dist = annDist(dim, node_bnds.bnd_box->lo, node_bnds.centroid);
    node_bnds.lambda = hi_dist > lo_dist ? hi_dist : lo_dist; // upper bound
    node_bnds.rho = node_bnds.lambda; // use same bound

    // Save bounding information in the hash map
    bounds->insert(std::pair< ANNkd_ptr, const Bound& >(ptr, node_bnds));
    return ptr; // return pointer to this node
  }
}
