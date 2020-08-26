#include <Rcpp.h>
using namespace Rcpp;

double dot_col(const NumericMatrix& P, int j, const NumericVector& x) {
  
  int size = x.size();
  double cp = 0;
  
  for (int i = 0; i < size; i++) {
    cp += P(i, j) * x[i];
  }
  
  return cp;
} 

// [[Rcpp::export]]
NumericMatrix rtmvnorm_gibbs_cpp(const NumericMatrix& P,
                                 const NumericVector& sd,
                                 const NumericVector& lower,
                                 const NumericVector& upper,
                                 const LogicalVector& fixed,
                                 const IntegerVector& to_return,
                                 NumericVector& x,
                                 int n_sim,
                                 int burn_in) {
  
  // Rcout << sum(fixed) << std::endl;
  // to_return <- rep(-1, d)
  // to_return[ind] <- seq_along(ind) - 1L
  
  
  NumericMatrix res(n_sim, Rcpp::sum(to_return >= 0));
  
  int n_total = n_sim + burn_in;
  int d = sd.size();
  
  for (int k = 0; k < n_total; k++) {
    for (int j = 0; j < d; j++) {
      if (!fixed[j]) {
        double mu_j = dot_col(P, j, x);
        double Fa = ::Rf_pnorm5(lower[j], mu_j, sd[j], 1, 0);
        double Fb = ::Rf_pnorm5(upper[j], mu_j, sd[j], 1, 0);
        double U_ab = ::Rf_runif(Fa, Fb);
        x[j] = ::Rf_qnorm5(U_ab, mu_j, sd[j], 1, 0);
        if (k >= burn_in && to_return[j] >= 0) {
          res(k - burn_in, to_return[j]) = x[j];
        }
      }
    }
  }
  
  return res;
}