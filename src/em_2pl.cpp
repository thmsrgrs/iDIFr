// em_2pl.cpp
// C++ kernels for the 2PL EM algorithm.
//
// Two exported functions:
//   compute_log_lik_matrix  -- E-step log-likelihood matrix
//   nr_update_pooled        -- M-step Newton-Raphson for shared item params
//
// Loop order is chosen for column-major (Fortran) storage: the innermost
// loop always varies the row index so adjacent elements are contiguous.

#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// compute_log_lik_matrix
//
// Returns an n_persons x n_nodes matrix L where
//   L[i,k] = sum_j { obs[i,j] * [ x[i,j]*log(P) + (1-x[i,j])*log(1-P) ] }
// and P = sigmoid( a_i[i,j] * (nodes_g[g[i]-1, k] - b_i[i,j]) ).
//
// Args:
//   resp_safe  n_persons x n_items  0/1 (NAs already replaced with 0)
//   obs_mat    n_persons x n_items  logical, TRUE where observed
//   a_i        n_persons x n_items  discrimination (a_mat[g_int, ])
//   b_i        n_persons x n_items  difficulty     (b_mat[g_int, ])
//   nodes_g    n_groups  x n_nodes  scaled quadrature nodes per group
//   g_int      n_persons            group index, 1-based
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix compute_log_lik_matrix(NumericMatrix  resp_safe,
                                     LogicalMatrix  obs_mat,
                                     NumericMatrix  a_i,
                                     NumericMatrix  b_i,
                                     NumericMatrix  nodes_g,
                                     IntegerVector  g_int) {

  const int n_persons = resp_safe.nrow();
  const int n_items   = resp_safe.ncol();
  const int n_nodes   = nodes_g.ncol();

  NumericMatrix loglik(n_persons, n_nodes);   // initialised to 0

  // Loop order: item (j) -> node (k) -> person (i)
  //   - For fixed j, resp_safe(*, j), obs_mat(*, j), a_i(*, j), b_i(*, j)
  //     are contiguous columns -- cache-friendly reads.
  //   - For fixed k, loglik(*, k) is a contiguous column -- cache-friendly
  //     accumulation.

  for (int j = 0; j < n_items; j++) {
    for (int k = 0; k < n_nodes; k++) {
      for (int i = 0; i < n_persons; i++) {
        if (!obs_mat(i, j)) continue;

        int    gi    = g_int[i] - 1;                // 0-based group index
        double theta = nodes_g(gi, k);
        double eta   = a_i(i, j) * (theta - b_i(i, j));
        double P     = 1.0 / (1.0 + std::exp(-eta));

        if (P < 1e-10)       P = 1e-10;
        else if (P > 1-1e-10) P = 1.0 - 1e-10;

        double x = resp_safe(i, j);
        loglik(i, k) += x * std::log(P) + (1.0 - x) * std::log(1.0 - P);
      }
    }
  }

  return loglik;
}


// ---------------------------------------------------------------------------
// nr_update_pooled
//
// Newton-Raphson update for item discrimination (a) and difficulty (b)
// pooled across all groups.  Called once per item per EM iteration.
//
// Args:
//   x         n_obs              binary responses (0/1)
//   post      n_obs x n_nodes    posterior weights
//   nodes_ik  n_obs x n_nodes    each row = that person's scaled quad nodes
//   a, b                         current parameter values
//   fix_a, fix_b                 if TRUE the respective param is not updated
//   max_nr, tol_nr               NR stopping criteria
//
// Returns: length-2 vector c(new_a, new_b)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector nr_update_pooled(NumericVector x,
                               NumericMatrix post,
                               NumericMatrix nodes_ik,
                               double a,
                               double b,
                               bool   fix_a  = false,
                               bool   fix_b  = false,
                               int    max_nr = 20,
                               double tol_nr = 1e-5) {

  const int n_obs   = x.size();
  const int n_nodes = post.ncol();

  for (int nr = 0; nr < max_nr; nr++) {

    double g_a = 0.0, h_a = 0.0;
    double g_b = 0.0, h_b = 0.0;

    // Loop order: node (k) -> person (i)
    //   post(*, k) and nodes_ik(*, k) are contiguous columns.
    for (int k = 0; k < n_nodes; k++) {
      for (int i = 0; i < n_obs; i++) {
        double theta = nodes_ik(i, k);
        double Z     = theta - b;
        double eta   = a * Z;
        double P     = 1.0 / (1.0 + std::exp(-eta));

        if (P < 1e-10)        P = 1e-10;
        else if (P > 1-1e-10) P = 1.0 - 1e-10;

        double D     = P * (1.0 - P);
        double w     = post(i, k);
        double resid = x[i] - P;

        if (!fix_a) {
          g_a += w * resid * Z;
          h_a -= w * D * Z * Z;
        }
        if (!fix_b) {
          g_b += w * resid * (-a);
          h_b -= w * D * a * a;
        }
      }
    }

    if (!fix_a && std::abs(h_a) > 1e-10) a -= g_a / h_a;
    if (!fix_b && std::abs(h_b) > 1e-10) b -= g_b / h_b;

    if (std::abs(g_a) + std::abs(g_b) < tol_nr) break;
  }

  return NumericVector::create(a, b);
}
