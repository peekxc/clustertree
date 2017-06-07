#include <Rcpp.h>
using namespace Rcpp;

// // Includes
// #include "union_find.h"
// #include "utilities.h"
//
// // Retrieve the full integer vector representing component membership (by reference)
// void CC(UnionFind& uf, IntegerVector& cc, std::vector<int>& admitted){
//   for (std::vector<int>::iterator it = admitted.begin(); it < admitted.end(); ++it){
//     cc[*it] = uf.Find(*it);
//   }
// }
//
// // Retrieve the full integer vector representing component membership
// IntegerVector getCC(UnionFind& uf, std::vector<bool>& admitted){
//   IntegerVector cc = Rcpp::no_init(uf.size);
//   for (int i = 0; i < uf.size; ++i){ cc[i] = admitted[i] ? uf.Find(i) : 0; }
//   return(cc);
// }
//
// void CC_noise(UnionFind& uf, std::vector<bool>& admitted){
//   // for (std::vector<int>::iterator it = admitted.begin(); it < admitted.end(); ++it){
//   //   cc[*it] = uf.Find(*it);
//   // }
//   for (int i = 0; i < uf.size; ++i){
//     Rcout << (admitted[i] ? uf.Find(i) : 0) << " ";
//   }
// }
//
// void printCC(UnionFind& uf, std::vector<bool>& admitted){
//   // for (std::vector<int>::iterator it = admitted.begin(); it < admitted.end(); ++it){
//   //   cc[*it] = uf.Find(*it);
//   // }
//   for (int i = 0; i < uf.size; ++i){
//     Rcout << (admitted[i] ? uf.Find(i) : 0) << " ";
//   }
// }
//
// // Retrieve the full integer vector representing component membership (by reference)
// void CC(UnionFind& uf, IntegerVector& cc, std::vector<bool>& admitted){
//   for (int i = 0; i < uf.size; ++i){
//     cc[i] = admitted[i] ? uf.Find(i) : 0;
//   }
// }
//
//
// /* -- CC_changed --
// * This function controls for the detection of a change between two sets of connected
// * components, restricted to only look for changes between points that have been admitted into
// * G_r.
// */
// bool CC_changed(IntegerVector current_CC, IntegerVector previous_CC, std::vector<int>& admitted){
//   if (admitted.size() < 2) return(false);
//   for (std::vector<int>::iterator it = admitted.begin(); it != admitted.end(); ++it){
//     if (current_CC[*it] != previous_CC[*it]) {// && contains(admitted, current_CC[*it])){
//       return(true);
//     }
//   }
//   return(false);
// }
//
// // Returns the indices of admitted points that were newly included
// IntegerVector which_added(IntegerVector current_CC, IntegerVector previous_CC,
//                           std::vector<int>& admitted,
//                           std::vector<bool>& included){
//   IntegerVector changed = IntegerVector();
//   for (std::vector<int>::iterator it = admitted.begin(); it != admitted.end(); ++it){
//     if (current_CC[*it] != previous_CC[*it] && included[*it] == 0) {// && contains(admitted, current_CC[*it])){
//       changed.push_back(*it);
//     }
//   }
//   return(changed);
// }
//
// IntegerMatrix mergeOrder(const NumericVector x){
//   const int n = as<int>(x.attr("Size"));
//   Function order = Function("order");
//   IntegerVector r_order = as<IntegerVector>(order(x)) - 1;
//   IntegerMatrix merge_order = IntegerMatrix(x.length(), 2);
//   int i = 0;
//   for(NumericVector::const_iterator r = x.begin(); r != x.end(); ++r, ++i){
//     int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
//     merge_order(i, _) = IntegerVector::create(from, to);
//   }
//   return(merge_order);
// }
//
// List checkKruskals(const NumericVector x, const NumericVector r_k, const int k, const double alpha = 1.414213562373095){
//
//   // Number of data points
//   const int n = as<int>(x.attr("Size"));
//
//   // Get sorted radii
//   NumericVector radii = Rcpp::clone(x).sort(false) * alpha;
//   NumericVector l2 = Rcpp::clone(x).sort(false) ;
//
//   // Get order of original; use R order function to get consistent ordering
//   Function order = Function("order");
//   IntegerVector r_order = as<IntegerVector>(order(x)) - 1;
//
//   // Set up connected components (initialized as singletons)
//   UnionFind components = UnionFind(n);
//
//   // Set up vector representing inclusion and admittance criteria
//   // admittance := only consider x_i s.t. \{ x_i : r_k(x_i) \leq r \}
//   // included := 0 if x_i is a singleton, 1 if x_i is connected to some component
//   std::vector<bool> admitted = std::vector<bool>(n, false);
//   std::vector<int> admitted2 = std::vector<int>();
//   std::vector<bool> included = std::vector<bool>(n, false);
//
//   // Admission tracker; rather than check admissability of every point at every radii,
//   // sort by increasing radius and make one comparison each loop
//   int c_rk = 0;
//   NumericVector sorted_r_k = clone(r_k).sort(false);
//   IntegerVector order_r_k = as<IntegerVector>(order(r_k)) - 1;
//
//   // Track the 0-based sorted lambda index distinct CCs occurred at
//   IntegerVector index_conjoined = IntegerVector(n, 0);
//
//   int i = 0, c = 0;
//   // Data structures needed to build the 'hclust' object
//   IntegerVector component_index = IntegerVector(n);
//
//   // Set up vector whose values, indexed by component, represent last merge step the component was a part of
//   IntegerMatrix merge = IntegerMatrix(n - 1, 2);
//   NumericVector height = NumericVector(n - 1);
//   IntegerMatrix merge_order = IntegerMatrix(radii.length(), 2);
//
//   List cc_list = List(n, IntegerVector());
//
//   IntegerVector previous_CC = seq(0, n - 1);
//   IntegerVector current_CC = clone(previous_CC);
//
//   IntegerVector n_admitted = IntegerVector(n, 0);
//   NumericVector R = NumericVector(n);
//
//   // Iterate r from 0 \to \infty
//   std::list<double_edge> to_connect = std::list<double_edge>();
//   for(NumericVector::iterator r = radii.begin(), dist_ij = l2.begin(); r != radii.end(); ++r, ++i, ++dist_ij){
//     int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
//     //int from_comp = components.Find(from), to_comp = components.Find(to);
//     merge_order(i, _) = IntegerVector::create(from, to);
//
//     // Sets admission status of points
//     while(c_rk < n && sorted_r_k.at(c_rk) <= (*dist_ij)) {
//       admitted.at(int(order_r_k.at(c_rk))) = true;
//       admitted2.push_back(int(order_r_k.at(c_rk)));
//       //Rcout << i << ": Set " << int(order_r_k.at(c_rk)) << " admitted == TRUE (@radius == " << *dist_ij << ")" << std::endl;
//       //Rcout << "HEAD: " << sorted_r_k.at(c_rk) << ", " << sorted_r_k.at(c_rk+1) << ", " << sorted_r_k.at(c_rk+2) << std::endl;
//       c_rk++;
//       // if (i == 33 || i == 34 || i == 35){
//       //   Rcout << "N admitted: " << admitted2.size() << std::endl;
//       // }
//       // Reconnect edges with distances < the current radius that have now been admitted to G_r
//       for (std::list<double_edge>::iterator it = to_connect.begin(); it != to_connect.end();){
//         double_edge hidden_edge = *it;
//         if (hidden_edge.from == 64 && hidden_edge.to == 53){
//           Rcout << "admitted index: " << c_rk << std::endl;
//           Rcout << "index: " << i << std::endl;
//           Rcout << "from: " << hidden_edge.from << " to: " << hidden_edge.to << std::endl;
//           Rcout << "admitted from: " << admitted[hidden_edge.from] << " to: " << admitted[hidden_edge.to] << std::endl;
//         }
//         if (admitted[hidden_edge.from] && admitted[hidden_edge.to]){
//           components.Union(hidden_edge.from, hidden_edge.to);
//           to_connect.erase(it++);
//         } else { it++; }
//       }
//
//     }
//
//     // if (i == 81){
//     //   Rcout << "from: " << from << " to: " << to << std::endl;
//     //   Rcout << "admitted from: " << admitted[from] << " to: " << admitted[to] << std::endl;
//     // }
//
//     // Union points which have been admitted, (potentially) mutating the current set of connected components
//     if (admitted.at(from) && admitted.at(to)){
//       //Rcout << "Unioning: " << from << ", " << to << std::endl;
//       components.Union(from, to);
//     } else
//       // If one or more of the points is not current in G_r, need to save it to a list to check for admissability
//       // later as the growing radius approaches varying neighorhood sizes
//     {
//       if (i == 81){
//         Rcout << "Pushing back " << from << ", " << to << std::endl;
//       }
//       to_connect.push_back(double_edge(from, to, *r));
//     }
//
//
//
//     // Update current connected components
//     //CC(components, current_CC, admitted);
//     CC(components, current_CC, admitted2);
//
//     if (i == 123 || i == 124 || i == 125){
//       Rcout << i << ": ";
//       printCC(components, admitted);
//       Rcout << std::endl;
//       Rcout << "from: " << from << ", to: " << to << std::endl;
//       // Rcout << "Admitted: ";
//       // for (int ii = 0; ii < n; ++ii){
//       //   Rcout << admitted.at(ii);
//       // }
//       // Rcout << std::endl << "radius: " << *r << std::endl;
//       // Rcout << "From comp: " << from_comp << ", To comp: " << to_comp << std::endl;
//       // Rcout << "Admitted? " << admitted.at(from) << " and " << admitted.at(to) << std::endl;
//       Rcout << "Component change? " << CC_changed(current_CC, previous_CC, admitted2) << std::endl;
//     }
//
//     // If the component has changed, updated with the earliest index that happened
//     if (CC_changed(current_CC, previous_CC, admitted2)){
//       n_admitted.at(c) = c_rk;
//       cc_list.at(c) = getCC(components, admitted); //clone(current_CC);
//       index_conjoined.at(c) = i;
//       R.at(c) = (*dist_ij);
//       c++;
//     }
//
//
//     // Construct a graph G_r with nodes \{ x_i : r_k(x_i) \leq r \}.
//     // if (admitted.at(from) && admitted.at(to)){
//
//     // Include edge (x_i, x_j) if \lVert x_i - x_j \rVert \leq \alpha r.
//     // If the edge inclusion results in a component change, record the change
//     // if (CC_changed(current_CC, previous_CC, admitted2)){
//     //   Rcout << i << ": Including edge: " << from << " -- " << to << " | ";
//     //   printCC(components, admitted);
//     //   Rcout << std::endl;
//     //   index_conjoined.at(c) = i;
//     //   height.at(c) = (*r);
//     //   R.at(c) = (*dist_ij);
//     //
//     //   // Save current connected components
//     //   cc_list.at(c) = getCC(components, admitted); //clone(current_CC);
//     //
//     //   IntegerVector noise = which_added(current_CC, previous_CC, admitted2, included);
//     //
//     //   for (IntegerVector::iterator new_pt = noise.begin(); new_pt != noise.end(); ++new_pt){
//     //     int cFrom = *new_pt; // cFrom := the newly included point that used to be not included in G_r
//     //     int cTo = current_CC[cFrom]; // cTo := the component cFrom belongs too
//     //
//     //     // Hclust class requirements: when building the merge matrix, positive merge indices represent agglomerations,
//     //     // and negative indicate singletons. This requires more indexing vectors to enable tracking component
//     //     // membership changes
//     //     if (!included.at(cFrom) && !included.at(cTo)) {
//     //       merge(c++, _) = IntegerVector::create(-(cFrom + 1), -(cTo + 1));
//     //       included.at(cFrom) = included.at(cTo) = true;
//     //     } else if (!included.at(cFrom) || !included.at(cTo)){
//     //       int leaf = (!included.at(cFrom)) ? cFrom : cTo;
//     //       int comp = (!included.at(cFrom)) ? to_comp : from_comp;
//     //       merge(c++, _) = IntegerVector::create(-(leaf + 1), component_index.at(comp));
//     //       included.at(leaf) = true;
//     //     } else {
//     //       merge(c++, _) = IntegerVector::create(component_index.at(from_comp), component_index.at(to_comp));
//     //     }
//     //     component_index.at(from_comp) = component_index.at(to_comp) = c;
//     //   }
//
//     //}
//     previous_CC = clone(current_CC);
//     // components.Union(from, to);
//   }
//   //}
//   // IntegerVector final_component = IntegerVector(n, 0);
//   // for (int i = 0; i < n; ++i){
//   //   final_component.at(i) = components.Find(i);
//   // }
//
//
//   // Extractor merge order and return
//   List res = List::create(
//     //_["merge"] = merge,
//     //_["height"] = height,
//     //_["order"] = extractOrder(merge),
//     _["labels"] = R_NilValue,
//     _["i"] = index_conjoined,
//     _["mo"] = merge_order,
//     _["cc"] = cc_list,
//     _["n_admitted"] = n_admitted,
//     _["dist_ij"] = R
//   );
//   res.attr("class") = "hclust";
//
//   return(res);
// }