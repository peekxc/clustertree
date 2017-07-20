#include "ANN/ANN.h"
// #include "dt_knn.h"

// struct DualTreeKNN; //forward declaration

// Forwarder to forward any dual tree subclass's score function to the appropriate context
// template <typename DT>
// void score_forwarder(void* context, ANNkd_node* lhs, ANNkd_node* rhs) {
//   static_cast<DT*>(context)->Score(lhs, rhs);
// }

// Trivial struct to use in the sorting network. Allows quick prioritization of recursive
// branches to follow.
struct NodeScore{
  ANNkd_node* lhs, *rhs;
  ANNdist score;
  NodeScore(ANNkd_node* l, ANNkd_node* r, ANNdist lr_score) : lhs(l), rhs(r) { score = lr_score; }
};

// Simple 4-part sorting network
// This sorting network will sort a fixed 4-length array of node scores very quickly,
// however it does NOT do bounds checking.
static inline void sort4_sn(NodeScore* d){
#define min(x, y) (x.score<y.score?x:y)
#define max(x, y) (x.score<y.score?y:x)
#define SWAP(x,y) { \
  const NodeScore a = min(d[x], d[y]);            \
  const NodeScore b = max(d[x], d[y]);           \
  d[x] = a; d[y]= b; }
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(0, 2);
  SWAP(1, 3);
  SWAP(1, 2);
#undef SWAP
#undef min
#undef max
}

// Simple 4-part sorting network
static inline void sort4_sn_stable(NodeScore* d){
#define min(x, y) (x.score<y.score?x:y)
#define max(x, y) (x.score<y.score?y:x)
#define SWAP(x,y) {                               \
  const NodeScore a = min(d[x], d[y]);            \
  const NodeScore b = max(d[x], d[y]);            \
  d[x] = a; d[y]= b; }
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(0, 2);
  SWAP(1, 3);
  SWAP(1, 2);
  SWAP(0, 1);
  SWAP(2, 3);
  SWAP(0, 2);
  SWAP(1, 3);
  SWAP(1, 2);
#undef SWAP
#undef min
#undef max
}


// Simple 3-part sorting network
static inline void sort3_sn(NodeScore* d){
#define min(x, y) (x.score<y.score?x:y)
#define max(x, y) (x.score<y.score?y:x)
#define SWAP(x,y) {                               \
  const NodeScore a = min(d[x], d[y]);            \
  const NodeScore b = max(d[x], d[y]);            \
  d[x] = a; d[y]= b; }
  SWAP(1, 2);
  SWAP(0, 2);
  SWAP(0, 1);
#undef SWAP
#undef min
#undef max
}

// Trivial 2-part sorting network just to keep syntax clean
static inline void sort2_sn(NodeScore* d){
#define min(x, y) (x.score<y.score?x:y)
#define max(x, y) (x.score<y.score?y:x)
#define SWAP(x,y) {                               \
  const NodeScore a = min(d[x], d[y]);            \
  const NodeScore b = max(d[x], d[y]);            \
  d[x] = a; d[y]= b; }
  SWAP(0, 1);
#undef SWAP
#undef min
#undef max
}
