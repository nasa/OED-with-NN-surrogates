#include <Rcpp.h>
using namespace Rcpp;

#include <limits>

// [[Rcpp::export]]
int closest_sample(const NumericVector x, const NumericMatrix& constituent_samples){
  unsigned int n = constituent_samples.nrow(), i = 0, idx = 0;
  double d = std::numeric_limits<double>::max(), tmp_d;
  for (i = 0; i < n; ++i){
    const Rcpp::NumericVector& v1 = constituent_samples.row(i);
    tmp_d = sum(pow(v1-x, 2.0));
    if (tmp_d < d){
      d = tmp_d;
      idx = i; 
    }
  }
  return(idx);
}

// NumericVector cmi(const NumericMatrix& X, const NumericMatrix& Y, const NumericMatrix& Z, const int k, const double alpha){
//   
// }

// [[Rcpp::export]]
NumericVector utility(const NumericMatrix& d, const List& B) {
  // 
  // gen_fiber_p <- d[, 1:2]
  // NumericMatrix fiber_prop = d(_, Rcpp::Range(2)); 
  // gen_matrix_p <- d[, 3:4]
  // gen_fvf <- t(sapply(d[, 5], function(vf) rep(vf, nrow(lam_layups))))
  // gen_layup <- truncate_between(d[, 6], 0 , 1)
}



/*** R


*/