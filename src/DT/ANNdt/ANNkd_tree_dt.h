#ifndef ANN_KD_TREE_DT_H
#define ANN_KD_TREE_DT_H

#include "Bound.h" // Bounds object + ANN includes
#include <unordered_map> // Hash table maps node pointers to their corresponding bounds
#include <ANN/kd_tree/kd_tree.h> // Base class
#include <ANN/kd_tree/kd_split.h>	// kd-tree splitting rules
#include <ANN/kd_tree/kd_util.h> // Utilities + Need for bounding box computations
#include <utilities.h> // R_PRINTF
#include <metrics.h> // Metrics

// No searching is done via this tree directly! This inherited relationship only exists to provide an
// augmented procedure for building a kd tree suitable for dual tree traversals. Namely, the exact procedure in which
// the tree is built is the same, depending on the splitting scheme, however the bounds traditionally computed
// and deallocated in the normal tree construction are saved here to allow bounding-computations later on.
class DLL_API ANNkd_tree_dt : public ANNkd_tree {
  friend class DualTree;
public:
  ANNkd_tree_dt(){}; // default constructor
  ANNkd_tree_dt(							// build from point array
    ANNpointArray		pa,				// point array (with at least n pts)
    int					n,				// number of points
    int					dd,				// dimension
    Metric& m,            // The metric distance to build the tree on
    std::unordered_map<ANNkd_node*, const Bound& >& bounds, // Map of the bounds to compute recursively in the tree construction
    int					bs = 1,				// bucket size
    ANNsplitRule		split = ANN_KD_SUGGEST // splitting method
    );
  ~ANNkd_tree_dt();						// tree destructor
};

//----------------------------------------------------------------------
//		External entry points
//----------------------------------------------------------------------
ANNkd_ptr rkd_tree_pr(				// recursive construction of kd-tree
    ANNpointArray		pa,				// point array (unaltered)
    ANNidxArray			pidx,			// point indices to store in subtree
    int					n,				// number of points
    int					dim,			// dimension of space
    Metric& m,            // The metric distance to build the tree on
    int					bsp,			// bucket space
    ANNorthRect			&bnd_box,		// bounding box for current node
    ANNkd_splitter		splitter, 	// splitting routine
    std::unordered_map<ANNkd_node*, const Bound& >& bounds // additional argument to store bound objects
);

#endif