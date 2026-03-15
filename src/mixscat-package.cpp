// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  // Ensure symmetry
  arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma_sym);
}

// [[Rcpp::export]]
arma::mat cpp_compute_V(arma::mat X,
                        arma::vec omega,
                        arma::mat precision_matrix) {

  arma::mat XtOmegaX = X.each_col() % omega;
  XtOmegaX = XtOmegaX.t() * X;
  XtOmegaX += 1e-7 * arma::eye(X.n_cols, X.n_cols);
  arma::mat V = arma::inv(XtOmegaX + precision_matrix);

  return V;
}

// [[Rcpp::export]]
arma::mat cpp_compute_m(arma::mat V,
                        arma::mat X,
                        arma::vec z,
                        arma::vec omega,
                        arma::vec C) {
  arma::mat m = V * (X.t() * ((z - 0.5) + (omega % C)));
  return m;
}

// [[Rcpp::export]]
arma::mat sample_beta(arma::mat X,
                      arma::vec omega,
                      arma::mat inv_cov,
                      arma::vec z,
                      arma::vec C) {

  arma::mat V = cpp_compute_V(X, omega, inv_cov);
  arma::mat mu = cpp_compute_m(V, X, z, omega, C);

  return mvrnormArma(1, mu, V).t();

}

// [[Rcpp::export]]
arma::mat sample_beta2(arma::mat X,
                       arma::vec omega,
                       arma::mat precision_matrix,
                       arma::vec z,
                       arma::vec C) {

  arma::mat X_omega = X.each_col() % omega;
  arma::mat post_prec = (X_omega.t() * X) + precision_matrix;

  arma::mat L_prec = arma::chol(post_prec, "lower");

  arma::vec rhs = X.t() * ((z - 0.5) + (omega % C));
  arma::vec std_noise = arma::randn<arma::vec>(rhs.n_elem);

  return arma::solve(arma::trimatl(L_prec), rhs + std_noise);

}

// [[Rcpp::export]]
double logsumexp_cpp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}


// [[Rcpp::export]]
arma::mat fast_dummy_dense(arma::ivec x, int G) {
  int N = x.n_elem;
  arma::mat out(N, G, arma::fill::zeros);

  for(int i = 0; i < N; ++i) {
    out(i, x(i)-1) = 1.0;
  }

  return out;
}

// [[Rcpp::export]]
arma::mat sofmax_cpp(const arma::mat& x) {
  arma::mat ex = arma::exp(x);
  arma::vec row_sum = arma::sum(ex, 1);
  ex.each_col() /= row_sum;
  arma::mat P = ex;
  return P;
}

// [[Rcpp::export]]
arma::mat compute_prob_group(arma::mat& B,
                             arma::mat& beta_group,
                             arma::uvec& idx) {
  arma::mat prob = sofmax_cpp(B * beta_group);
  return prob.rows(idx);
}

// [[Rcpp::export]]
arma::mat predict_prob_cpp(int& M,
                           arma::ivec& w,
                           arma::mat& B,
                           arma::mat& beta) {

  arma::mat W = fast_dummy_dense(w, M);
  W.col(0) = 1.0;
  arma::mat X = arma::kron(W, B);
  arma::mat prob = sofmax_cpp(X * beta);
  return prob;
}

// [[Rcpp::export]]
arma::vec fast_aggregate_sum(arma::vec& log_pz, arma::ivec& id) {
  // 1. Find the range of IDs
  int min_id = id.min();
  int max_id = id.max();
  int range = max_id - min_id + 1;

  // 2. Initialize a result vector with zeros
  arma::vec sums(range, fill::zeros);

  // 3. Single-pass accumulation (The "O(n)" magic)
  for (unsigned int i = 0; i < id.n_elem; ++i) {
    // Offset by min_id so it starts at index 0
    sums(id(i) - min_id) += log_pz(i);
  }

  return sums;
}


