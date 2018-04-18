#include "metrics.h"

// SEXP chooseMetric_int(std::string metric_name, List config) {
//   if (metric_name == "euclidean")     {
//     L_2* metric = new L_2();
//     if (!config.containsElementNamed("d")){
//       Rcpp::stop(std::string("Invalid parameters given for: ")+metric_name);
//     }
//     metric->init(config);
//     Rcpp::XPtr< L_2 > p((L_2*) metric, true);
//     return wrap(p);
//   }
//   else if (metric_name == "manhattan"){
//     L_1* metric = new L_1();
//
//     if (!config.containsElementNamed("d")){
//       Rcpp::stop(std::string("Invalid parameters given for: ")+metric_name);
//     }
//
//     metric->init(config);
//     Rcpp::XPtr< L_1 > p((L_1*) metric, true);
//     return wrap(p);
//   }
//   else if (metric_name == "maximum")  {
//     L_inf* metric = new L_inf();
//
//     if (!config.containsElementNamed("d")){
//       Rcpp::stop(std::string("Invalid parameters given for: ")+metric_name);
//     }
//
//     metric->init(config);
//     Rcpp::XPtr< L_inf > p((L_inf*) metric, true);
//     return wrap(p);
//   }
//   else if (metric_name == "minkowski")  {
//     L_p* metric = new L_p();
//
//     if (!config.containsElementNamed("d") || !config.containsElementNamed("p")){
//       Rcpp::stop(std::string("Invalid parameters given for: ")+metric_name);
//     }
//
//     metric->init(config);
//     Rcpp::XPtr< L_p > p((L_p*) metric, true);
//     return wrap(p);
//   }
//   return 0;
// }
//
// Metric& getMetric(SEXP metric_ptr){
//   Rcpp::XPtr< Metric > metric_xptr((SEXP) metric_ptr);
//   Metric& metric = *metric_xptr;
//   return metric;
// }

// Double-check slicing isn't occurring
// void testMetric(SEXP metric_exp){
//   Metric& metric = getMetric(metric_exp);
//   double x[2] = {1, 5};
//   double y[2] = {4, 6};
//   Rcout << "dist = " << metric((ANNpoint) &x, (ANNpoint) &y) << std::endl;
// }


