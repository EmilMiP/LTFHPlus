/******************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/******************************************************************************/

inline void myassert_size(size_t n1, size_t n2) {
  if (n1 != n2) 
    Rcpp::stop("Tested %s == %s. %s", n1, n2, "Incompatibility between dimensions.");
}

/******************************************************************************/

double dot_col(const NumericMatrix& P, int j, const NumericVector& x) {
  
  int size = x.size();
  double cp = 0;
  
  for (int i = 0; i < size; i++) {
    cp += P(i, j) * x[i];
  }
  
  return cp;
} 

/******************************************************************************/

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
  
  int d = sd.size();
  myassert_size(P    .nrow(), d);
  myassert_size(P    .ncol(), d);
  myassert_size(lower.size(), d);
  myassert_size(upper.size(), d);
  myassert_size(fixed.size(), d);
  myassert_size(x    .size(), d);
  
  NumericMatrix res(n_sim, Rcpp::sum(to_return >= 0));
  
  for (int k = -burn_in; k < n_sim; k++) {
    
    for (int j = 0; j < d; j++) {
      
      if (!fixed[j]) {
        
        double mu_j = dot_col(P, j, x);
        double Fa = ::Rf_pnorm5(lower[j], mu_j, sd[j], 1, 0);
        double Fb = ::Rf_pnorm5(upper[j], mu_j, sd[j], 1, 0);
        double U_ab = ::Rf_runif(Fa, Fb);
        
        x[j] = ::Rf_qnorm5(U_ab, mu_j, sd[j], 1, 0);
        
        if (k >= 0 && to_return[j] >= 0)
          res(k, to_return[j]) = dot_col(P, j, x);
      }
    }
  }
  
  return res;
}

/******************************************************************************/