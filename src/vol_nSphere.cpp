#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double vol_nSphere(const int n, const double R = 1) {
  switch(n){
    case 0:
      return 1.0;
    case 1:
      return 2.0 * R;
    case 2:
      return PI * std::pow(R, 2);
    case 3:
      return (4.0 * PI) / 3.0 * std::pow(R, 3);
    case 4:
      return std::pow(PI, 2) / 2.0 * std::pow(R, 4);
    case 5:
      return 8.0 * std::pow(PI, 2) / 15.0 * std::pow(R, 5);
    case 6:
      return std::pow(PI, 3) / 6.0 * std::pow(R, 6);
    case 7:
      return 16.0 * std::pow(PI, 3) / 105.0 * std::pow(R, 7);
    case 8:
      return std::pow(PI, 4) / 24.0 * std::pow(R, 8);
    default:
      return std::pow((double) PI, (n/2.0)) / tgamma((n/2.0) + 1.0) * std::pow(R, n);
  }
}
