// Algorithm 2: kNN Graph
// r_k := vector of k - 1 nearest neighbors
// knn_indices := id's of knn correspond to r_k
List kruskalsKNN(const NumericVector dist_x,
                 const NumericVector r_k,
                 const IntegerVector knn_indices,
                 const int n, const int k, const double alpha){

  Rcpp::Rcout << "Beginning kruskals" << std::endl;
  // Set up resulting MST
  NumericMatrix mst = NumericMatrix(n - 1, 3);

  // Get sorted radii
  NumericVector lambda = Rcpp::clone(r_k).sort(false);

  // Get order of original; use R order function to get consistent ordering
  Function order = Function("order"), duplicated = Function("duplicated");
  IntegerVector rk_order = as<IntegerVector>(order(r_k)) - 1;
  LogicalVector admitted = LogicalVector(n, false);

  IntegerVector x_order = as<IntegerVector>(order(dist_x)) - 1;
  NumericVector inc_dist = Rcpp::clone(dist_x).sort(false);

  // Create disjoint-set data structure to track components
  UnionFind components = UnionFind(n);
  UnionFind prev_components = UnionFind(n);
  UnionFind dist_components = UnionFind(n);

  bool connect_next = false;
  int i = 0, cr_k = 0, crow = 0, px_i = INDEX_TO(x_order.at(0), n), px_j = INDEX_FROM(x_order.at(0), n, px_i);
  for(NumericVector::iterator dist_ij = inc_dist.begin(); dist_ij != inc_dist.end(); ++dist_ij, ++i){
    int x_j = INDEX_TO(x_order.at(i), n), x_i = INDEX_FROM(x_order.at(i), n, x_j);
    if (dist_components.Find(x_i) != dist_components.Find(x_j)) dist_components.Union(x_i, x_j);
    while(cr_k < n - 1 && r_k.at(rk_order.at(cr_k)) <= (*dist_ij)) { admitted.at(int(rk_order.at(cr_k++))) = true; }
    if (connect_next){
      prev_components.Union(px_i, px_j);
      connect_next = false;
    }
    if ((*dist_ij)/alpha <= lambda.at(cr_k) && admitted.at(x_i) && admitted.at(x_j)){
      //IntegerVector CCs = components.getCC();
      components.merge(dist_components, admitted);
      if (components != prev_components)
      {
        mst(crow, _) = NumericVector::create(x_i, x_j, *dist_ij);
        crow++;
        components.Union(x_i, x_j);
        px_i = x_i, px_j = x_j;
        connect_next = true;
        if (crow == n) break;
      }
    }
  }

  // // Create disjoint-set data structure to track components
  // UnionFind components = UnionFind(n);
  // i = 0;
  // int crow = 0;
  // for (NumericVector::const_iterator r = lambda.begin(); r != lambda.end(); ++r, ++i) {
  //   if (i % 100 == 0) Rcpp::checkUserInterrupt();
  //
  //   // Point x_i becomes admitted into the graph G_r
  //   admitted.at(r_order.at(i)) = true;
  //
  //   // Retrieve index of x_i and x_j
  //   // if (i < 2){
  //   //   int from = INDEX_FROM_KNN(i+1, k);
  //   //   Rcout << "rorder_i: " << r_order.at(i) << std::endl;
  //   //   Rcout << "index from: " << ((int) from) << std::endl;
  //   //   Rcout << "to: " << knn_indices.at( r_order.at(i)) - 1 << std::endl;
  //   //   Rcout << "to_i: " <<  r_order.at(knn_indices.at( r_order.at(i)) - 1) << std::endl;
  //   // }
  //   int x_i = r_order.at(i); //(int) INDEX_FROM_KNN(i+1, k);
  //   int x_j = int(knn_indices.at(x_i) - 1);
  //   // int to = INDEX_TO(r_order.at(i), n), from = INDEX_FROM(r_order.at(i), n, to);
  //
  //
  //   if (admitted.at(x_i) && admitted.at(x_j)) {
  //     if (components.Find(x_i) != components.Find(x_j)){
  //       mst(crow, _) = NumericVector::create(x_i, x_j, *r);
  //       crow++;
  //       components.Union(x_i, x_j);
  //       if (crow == n) break;
  //     }
  //   }
  // }

  return(List::create(_["mst"] = mst, _["admitted"] = admitted));
}