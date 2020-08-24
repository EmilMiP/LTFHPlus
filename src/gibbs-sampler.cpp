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
NumericVector rtmvnorm_gibbs_cpp(const NumericMatrix& P,
                                 const NumericVector& sd,
                                 const NumericVector& lower,
                                 const NumericVector& upper,
                                 const LogicalVector& fixed,
                                 NumericVector& x,
                                 int n_sim,
                                 int burn_in) {
  
 // Rcout << sum(fixed) << std::endl;
  
  std::vector<double> res;
  
  int n_total = n_sim + burn_in;
  int d = sd.size();
  
  for (int k = 0; k < n_total; k++) {
    
    // first dimension with no bounds
    double mu_0 = dot_col(P, 0, x);
    x[0] = ::Rf_rnorm(mu_0, sd[0]);
    
    if (k >= burn_in) res.push_back(mu_0);
    
    // all other dimensions
    for (int j = 1; j < d; j++) {
      if (!fixed[j]) {
        double mu_j = dot_col(P, j, x);
        double Fa = ::Rf_pnorm5(lower[j], mu_j, sd[j], 1, 0);
        double Fb = ::Rf_pnorm5(upper[j], mu_j, sd[j], 1, 0);
        double U_ab = ::Rf_runif(Fa, Fb);
        x[j] = ::Rf_qnorm5(U_ab, mu_j, sd[j], 1, 0);
      }
    }
  }
  
  return wrap(res);
}