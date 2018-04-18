#include "node_bnb.h"
#include "ANN_util.h"
#include "kd_tree.h"
#include "kd_util.h"
#include "kd_split.h"

#include "RcppHeader.h"

static int g_id = 0;

// List::create(_["id"] = g_id++, _["max_radius"] = max_radius, _["centroid"] = centroid, _["bound"] = ANN_DIST_INF,
//              _["lb"] = ptToVec(bnd_box.lo, dim), _["ub"] = ptToVec(bnd_box.hi, dim), _["idx"] = idxArrayToVec(pidx, n));


ANNkd_ptr rkd_tree_dt(				// recursive construction of kd-tree
    ANNpointArray		pa,				// point array
    ANNidxArray			pidx,			// point indices to store in subtree
    int					n,				// number of points
    int					dim,			// dimension of space
    int					bsp,			// bucket space
    ANNorthRect			&bnd_box,		// bounding box for current node
    ANNkd_splitter		splitter, // splitting routine
    NODE_INFO& node_info) // node info
{
  // Compute and store the maximum radius around the bounding box
  node_bnb c_node_info = node_bnb();
  c_node_info.id = g_id++;
  c_node_info.lo = annCopyPt(dim, bnd_box.lo); // for debugging
  c_node_info.hi = annCopyPt(dim, bnd_box.hi); // for debugging
  c_node_info.centroid = annAllocPt(dim, 0);
  for (int d_i = 0; d_i < dim; ++d_i){
    double delta_di = bnd_box.hi[d_i] - bnd_box.lo[d_i]; // box length in dimension d_i
    c_node_info.centroid[d_i] = bnd_box.lo[d_i] + delta_di/2; // centroid center is lies along the halfspace
  }
  // max radius is half of the squared euclidean radius of the smallest circumsphere enclosing the box
  c_node_info.max_radius = sqrt(annDist(dim, bnd_box.lo, bnd_box.hi))/2;

  // If is a leaf node
  if (n <= bsp) {						// n small, make a leaf node
    if (n == 0)						// empty leaf node
      return KD_TRIVIAL;			// return (canonical) empty leaf
    else {// construct the node and return
      ANNkd_leaf* new_leaf = new ANNkd_leaf(n, pidx);
      new_leaf->id = c_node_info.id;
      try { node_info.at(new_leaf->id) = c_node_info; }
      catch (std::out_of_range ex) {
        node_info.resize(new_leaf->id+1);
        node_info.at(new_leaf->id) = c_node_info;
      }
      // node_info[c_node_info.id] = c_node_info;
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

    ANNcoord lv = bnd_box.lo[cd];	// save bounds for cutting dimension
    ANNcoord hv = bnd_box.hi[cd];

    bnd_box.hi[cd] = cv;			// modify bounds for left subtree
    lo = rkd_tree_dt(					// build left subtree
      pa, pidx, n_lo,			// ...from pidx[0..n_lo-1]
      dim, bsp, bnd_box, splitter, node_info);
    bnd_box.hi[cd] = hv;			// restore bounds

    bnd_box.lo[cd] = cv;			// modify bounds for right subtree
    hi = rkd_tree_dt(					// build right subtree
      pa, pidx + n_lo, n-n_lo,// ...from pidx[n_lo..n-1]
      dim, bsp, bnd_box, splitter, node_info);
    bnd_box.lo[cd] = lv;			// restore bounds

    // create the splitting node
    ANNkd_split* ptr = new ANNkd_split(cd, cv, lv, hv, lo, hi);

    // The children should also be stored in the tree; track the parent
    ptr->id = c_node_info.id;
    c_node_info.parent = nullptr; // temporary
    try { node_info.at(ptr->id) = c_node_info; }
    catch (std::out_of_range ex) {
      node_info.resize(ptr->id+1);
      node_info.at(ptr->id) = c_node_info;
    }
    // node_info[lo->id].parent = ptr;
    // node_info[hi->id].parent = ptr;
    return ptr;						// return pointer to this node
  }
}

ANNkd_tree::ANNkd_tree(					// construct from point array
  ANNpointArray		pa,				// point array (with at least n pts)
  int					n,				// number of points
  int					dd,				// dimension
  int					bs,				// bucket size
  ANNsplitRule split, // splitting method
  NODE_INFO& node_info // Map to store various things per-node at tree construction time
){
  g_id = 0;
  SkeletonTree(n, dd, bs);			// set up the basic stuff
  pts = pa;							// where the points are
  if (n == 0) return;					// no points--no sweat

  ANNorthRect bnd_box(dd);			// bounding box for points
  annEnclRect(pa, pidx, n, dd, bnd_box);// construct bounding rectangle
  // copy to tree structure
  bnd_box_lo = annCopyPt(dd, bnd_box.lo);
  bnd_box_hi = annCopyPt(dd, bnd_box.hi);

  switch (split) {					// build by rule
  case ANN_KD_STD:					// standard kd-splitting rule
    root = rkd_tree_dt(pa, pidx, n, dd, bs, bnd_box, kd_split, node_info);
    break;
  case ANN_KD_MIDPT:					// midpoint split
    root = rkd_tree_dt(pa, pidx, n, dd, bs, bnd_box, midpt_split, node_info);
    break;
  case ANN_KD_FAIR:					// fair split
    root = rkd_tree_dt(pa, pidx, n, dd, bs, bnd_box, fair_split, node_info);
    break;
  case ANN_KD_SUGGEST:				// best (in our opinion)
  case ANN_KD_SL_MIDPT:				// sliding midpoint split
    root = rkd_tree_dt(pa, pidx, n, dd, bs, bnd_box, sl_midpt_split, node_info);
    break;
  case ANN_KD_SL_FAIR:				// sliding fair split
    root = rkd_tree_dt(pa, pidx, n, dd, bs, bnd_box, sl_fair_split, node_info);
    break;
  default:
    annError("Illegal splitting method", ANNabort);
  }
}

