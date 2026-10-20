// hpsreg computation engine
//
//  1. Hierarchy scoring kernel   (cpp_score_segment, cpp_all_pairs, cpp_index_to_pair)
//  2. Davidson model, stored pairs   (cpp_btd_eval, cpp_fit_btd, cpp_vcov_btd)
//  3. Cross-group win counts, stored pairs (cpp_cross_sums)
//  4. STREAMING engine: pairs are generated and scored inside the loops and are
//     never stored, so memory is O(N p) instead of O(N^2)
//        cpp_stream_counts, cpp_stream_fit, cpp_stream_vcov,
//        cpp_stream_cross, cpp_stream_calib
//  5. COMPACT engine: pairs are scored once into a 1-byte-per-pair code array
//     (outcome + resolving level; indices implicit from the loop order) that every
//     fit / variance pass then reads instead of re-scoring
//        cpp_compact_build, cpp_compact_fit, cpp_compact_vcov
//
// Model.  For a pair (i, j) with d = (x_i - x_j)' beta / 2 and delta = exp(gamma),
//   P(i wins) = e^d / q,  P(j wins) = e^-d / q,  P(tie) = delta / q,
//   q = e^d + e^-d + delta.
// With delta fixed at 0 (no tied pairs) this is the Bradley-Terry model.
//
// Nothing here forms the n_pairs x p matrix of covariate differences: linear
// predictors are computed per subject (N values) and differences are formed on the fly.
//
// Parallelism (streaming engine only).  Pairs are split into a fixed number of
// chunks that depends only on N. Each chunk has its own accumulators and chunks are
// combined in a fixed order, so results are bitwise identical for any number of
// threads, with or without OpenMP. OpenMP is used only if the compiler provides it.

// [[Rcpp::depends(RcppEigen)]]
#define EIGEN_NO_DEBUG 1
#define EIGEN_DONT_PARALLELIZE 1
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMat;

// [[Rcpp::export]]
int cpp_openmp_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
#endif
}

// ======================================================================
// 1. Hierarchy scoring kernel
// ======================================================================

namespace {

struct LevelSpec {
  int    type;     // 1 = numeric (binary / ordinal / continuous), 2 = time-to-event
  int    col;      // value column (numeric) or time column (tte), 0-based
  int    col2;     // tte: exact-observation indicator column; otherwise -1
  int    elig;     // eligibility column or -1
  int    na_rule;  // numeric: 0 = NA ties, 1 = NA is best, 2 = NA is worst
  double sign;     // +1 if larger value / later time is better, -1 otherwise
  double thr;      // tie margin (>= 0)
};

std::vector<LevelSpec> parse_levels(const List& levels) {
  std::vector<LevelSpec> out;
  out.reserve(levels.size());
  for (int k = 0; k < levels.size(); ++k) {
    List L = levels[k];
    LevelSpec s;
    s.type    = as<int>(L["type"]);
    s.col     = as<int>(L["col"]);
    s.col2    = as<int>(L["col2"]);
    s.elig    = as<int>(L["elig"]);
    s.na_rule = as<int>(L["na_rule"]);
    s.sign    = as<double>(L["sign"]);
    s.thr     = as<double>(L["thr"]);
    out.push_back(s);
  }
  return out;
}

// Returns +1 (a wins), -1 (b wins) or 0 (neutral / indeterminate).
// `indet` is set when a zero was caused by censoring (uninformative pair).
inline int eval_level(const LevelSpec& L, const double* V, R_xlen_t N,
                      int a, int b, bool& indet) {
  if (L.elig >= 0) {
    const double ea = V[(R_xlen_t)L.elig * N + a];
    const double eb = V[(R_xlen_t)L.elig * N + b];
    if (std::isnan(ea) || std::isnan(eb) || ea == 0.0 || eb == 0.0) return 0;
  }
  if (L.type == 1) {
    const double x = V[(R_xlen_t)L.col * N + a];
    const double y = V[(R_xlen_t)L.col * N + b];
    const bool nx = std::isnan(x), ny = std::isnan(y);
    if (nx && ny) return 0;
    if (nx || ny) {
      if (L.na_rule == 0) return 0;
      const bool na_wins = (L.na_rule == 1);        // NA treated as best
      if (nx) return na_wins ? 1 : -1;
      return na_wins ? -1 : 1;
    }
    const double diff = L.sign * (x - y);
    if (diff >  L.thr) return 1;
    if (diff < -L.thr) return -1;
    return 0;
  }
  // time-to-event (Gehan-type rule for right censoring)
  const double ta = V[(R_xlen_t)L.col  * N + a];
  const double tb = V[(R_xlen_t)L.col  * N + b];
  if (std::isnan(ta) || std::isnan(tb)) return 0;
  const bool xa = V[(R_xlen_t)L.col2 * N + a] > 0.5;   // exactly observed
  const bool xb = V[(R_xlen_t)L.col2 * N + b] > 0.5;
  int r = 0;                                           // from "later is better" view
  if (xa && xb) {
    if (ta == tb) return 0;
    const double diff = ta - tb;
    r = (diff > L.thr) ? 1 : ((diff < -L.thr) ? -1 : 0);
  } else if (xa && !xb) {            // b censored at tb
    if (tb - ta > L.thr) r = -1; else { indet = true; return 0; }
  } else if (!xa && xb) {            // a censored at ta
    if (ta - tb > L.thr) r = 1;  else { indet = true; return 0; }
  } else {
    indet = true;                    // both censored
    return 0;
  }
  return L.sign > 0 ? r : -r;
}

}  // namespace

// Score a segment of consecutive built-in levels on a set of pairs.
// id_i / id_j are 1-based subject indices; `active` (1-based, optional) selects
// the pairs to score (empty = all). Returns W in {-1,0,1}, the 1-based position of
// the resolving level inside the segment (0 = unresolved) and an indeterminacy flag.
// [[Rcpp::export]]
List cpp_score_segment(List levels, NumericMatrix V,
                       IntegerVector id_i, IntegerVector id_j,
                       IntegerVector active) {
  const std::vector<LevelSpec> L = parse_levels(levels);
  const int K = (int)L.size();
  const R_xlen_t N = V.nrow();
  const bool use_active = active.size() > 0;
  const R_xlen_t n = use_active ? (R_xlen_t)active.size() : (R_xlen_t)id_i.size();
  IntegerVector W(n), lev(n), unf(n);
  const double* Vp = V.begin();

  for (R_xlen_t m = 0; m < n; ++m) {
    const R_xlen_t k = use_active ? (R_xlen_t)(active[m] - 1) : m;
    const int a = id_i[k] - 1, b = id_j[k] - 1;
    int w = 0, lv = 0, uf = 0;
    for (int l = 0; l < K; ++l) {
      bool indet = false;
      const int r = eval_level(L[l], Vp, N, a, b, indet);
      if (r != 0) { w = r; lv = l + 1; break; }
      if (indet) uf = 1;
    }
    W[m] = w; lev[m] = lv; unf[m] = uf;
  }
  return List::create(_["W"] = W, _["level"] = lev, _["uninf"] = unf);
}

// All N(N-1)/2 pairs (i < j), 1-based, row-major order over i.
// [[Rcpp::export]]
List cpp_all_pairs(int N) {
  if (N < 2) stop("need at least 2 subjects");
  const double nd = 0.5 * (double)N * ((double)N - 1.0);
  if (nd > 4.5e15) stop("too many pairs");
  const R_xlen_t n = (R_xlen_t)nd;
  IntegerVector I(n), J(n);
  R_xlen_t k = 0;
  for (int i = 1; i < N; ++i)
    for (int j = i + 1; j <= N; ++j) { I[k] = i; J[k] = j; ++k; }
  return List::create(_["i"] = I, _["j"] = J);
}

// Map 1-based linear indices (as doubles, to allow > 2^31) of the all-pairs
// enumeration above back to (i, j).
// [[Rcpp::export]]
List cpp_index_to_pair(NumericVector k, int N) {
  const R_xlen_t n = k.size();
  IntegerVector I(n), J(n);
  const double Nd = (double)N;
  auto start = [&](long long r) -> long long {
    return r * (long long)(N - 1) - r * (r - 1) / 2;
  };
  for (R_xlen_t m = 0; m < n; ++m) {
    const long long kk = (long long)k[m] - 1;
    if (kk < 0 || kk >= (long long)(0.5 * Nd * (Nd - 1.0)))
      stop("pair index out of range");
    const double disc = (2 * Nd - 1) * (2 * Nd - 1) - 8.0 * (double)kk;
    long long r = (long long)std::floor(((2 * Nd - 1) - std::sqrt(std::max(disc, 0.0))) / 2.0);
    if (r < 0) r = 0;
    while (r > 0 && start(r) > kk) --r;
    while (start(r + 1) <= kk) ++r;
    I[m] = (int)(r + 1);
    J[m] = (int)(r + 2 + (kk - start(r)));
  }
  return List::create(_["i"] = I, _["j"] = J);
}


namespace {

// Davidson pair probabilities with a fast path and a log-sum-exp fallback.
struct PairProb { double lq, pw, pl, pt; };

inline PairProb pair_prob(double d, double gamma, double dl, bool est) {
  PairProb r;
  if (std::fabs(d) < 40.0 && (!est || std::fabs(gamma) < 40.0)) {
    const double u = std::exp(d), v = 1.0 / u;
    const double q = u + v + dl, iq = 1.0 / q;
    r.lq = std::log(q); r.pw = u * iq; r.pl = v * iq; r.pt = dl * iq;
  } else {
    const double NEG_INF = -std::numeric_limits<double>::infinity();
    const double m = std::max(std::fabs(d), est ? gamma : NEG_INF);
    double s = std::exp(d - m) + std::exp(-d - m);
    if (est) s += std::exp(gamma - m);
    r.lq = m + std::log(s);
    r.pw = std::exp(d - r.lq); r.pl = std::exp(-d - r.lq);
    r.pt = est ? std::exp(gamma - r.lq) : 0.0;
  }
  return r;
}

}  // namespace


// ======================================================================
// 2. Davidson model: accumulators shared by the stored and streaming engines
// ======================================================================

namespace {


// Accumulates log-likelihood, gradient and (optionally) information over pairs.
// NB: accumulators of neighbouring chunks live next to each other in a std::vector.
// The scalar fields are written for every pair, so they are padded to keep them on
// separate cache lines from those of other chunks (avoids false sharing between threads).
struct EvalAcc {
  char pad_head_[128];
  int p = 0; bool est = false, hess = false;
  double ll = 0.0, sg = 0.0, hgg = 0.0;
  std::vector<double> c, c2, Hbb, dx;
  char pad_tail_[128];

  void init(int N, int p_, bool est_, bool hess_) {
    p = p_; est = est_; hess = hess_;
    ll = sg = hgg = 0.0;
    c.assign(N, 0.0);
    if (hess) { c2.assign(N, 0.0); Hbb.assign((size_t)p * p, 0.0); dx.assign(p, 0.0); }
  }

  inline void add(int a, int b, int wk, const double* eta, const double* xr,
                  double gamma, double dl) {
    const double d = 0.5 * (eta[a] - eta[b]);
    const PairProb pr = pair_prob(d, gamma, dl, est);
    const double pw = pr.pw, pl = pr.pl, pt = pr.pt;
    const bool tie = (wk == 0);
    ll += d * wk + ((tie && est) ? gamma : 0.0) - pr.lq;
    const double ak = (double)wk - (pw - pl);
    c[a] += 0.5 * ak;
    c[b] -= 0.5 * ak;
    if (est) sg += (tie ? 1.0 : 0.0) - pt;
    if (hess) {
      const double f = 0.25 * (4.0 * pw * pl + pt * (pw + pl));
      const double* xa = xr + (size_t)a * p;
      const double* xb = xr + (size_t)b * p;
      for (int r = 0; r < p; ++r) dx[r] = xa[r] - xb[r];
      for (int r = 0; r < p; ++r) {
        const double cf = f * dx[r];
        double* row = &Hbb[(size_t)r * p];
        for (int s2 = r; s2 < p; ++s2) row[s2] += cf * dx[s2];
      }
      if (est) {
        hgg += pt * (pw + pl);
        const double tk = 0.5 * pt * (pw - pl);
        c2[a] += tk;
        c2[b] -= tk;
      }
    }
  }

  void merge(const EvalAcc& o) {
    ll += o.ll; sg += o.sg; hgg += o.hgg;
    for (size_t k = 0; k < c.size(); ++k) c[k] += o.c[k];
    if (hess) {
      for (size_t k = 0; k < c2.size(); ++k) c2[k] += o.c2[k];
      for (size_t k = 0; k < Hbb.size(); ++k) Hbb[k] += o.Hbb[k];
    }
  }
};

void finish_eval(const EvalAcc& A, const MatrixXd& X, bool est,
                 double& ll, VectorXd& g, MatrixXd& Hn) {
  const int N = (int)X.rows(), p = (int)X.cols();
  const int np = p + (est ? 1 : 0);
  ll = A.ll;
  const Eigen::Map<const VectorXd> c(A.c.data(), N);
  g.resize(np);
  g.head(p) = X.transpose() * c;
  if (est) g[p] = A.sg;
  if (A.hess) {
    Hn = MatrixXd::Zero(np, np);
    for (int r = 0; r < p; ++r)
      for (int s2 = r; s2 < p; ++s2) {
        Hn(r, s2) = A.Hbb[(size_t)r * p + s2];
        Hn(s2, r) = Hn(r, s2);
      }
    if (est) {
      const Eigen::Map<const VectorXd> c2(A.c2.data(), N);
      const VectorXd hbg = -(X.transpose() * c2);
      Hn.block(0, p, p, 1) = hbg;
      Hn.block(p, 0, 1, p) = hbg.transpose();
      Hn(p, p) = A.hgg;
    }
  }
}

// Accumulates per-subject influence contributions (sum of pair scores).
//
// Two entry points:
//   add()        — the original, used by cpp_vcov_btd (stored pair path).
//   add_to_row() — accumulates row a into a caller-owned local buffer and
//                  updates only row b in memory.  flush_row() writes the
//                  row-a buffer into IF once per outer i.  Used by the tiled
//                  vcov passes to keep the IF working set inside L2.
struct IFAcc {
  char pad_head_[128];
  int p = 0, np = 0; bool est = false;
  RowMat IF;
  char pad_tail_[128];

  void init(int N, int p_, bool est_) {
    p = p_; est = est_; np = p + (est ? 1 : 0);
    IF = RowMat::Zero(N, np);
  }

  inline void add(int a, int b, int wk, const double* eta, const double* xr,
                  double gamma, double dl) {
    const double d = 0.5 * (eta[a] - eta[b]);
    const PairProb pr = pair_prob(d, gamma, dl, est);
    const double ak = 0.5 * ((double)wk - (pr.pw - pr.pl));
    const double* xa = xr + (size_t)a * p;
    const double* xb = xr + (size_t)b * p;
    double* ia = IF.data() + (size_t)a * np;
    double* ib = IF.data() + (size_t)b * np;
    for (int r = 0; r < p; ++r) {
      const double v = ak * (xa[r] - xb[r]);
      ia[r] += v; ib[r] += v;
    }
    if (est) {
      const double v = (wk == 0 ? 1.0 : 0.0) - pr.pt;
      ia[p] += v; ib[p] += v;
    }
  }

  // Row-a is accumulated in `local` (length np, pre-zeroed by the caller);
  // only row b is written to memory here.
  inline void add_to_row(int a, int b, int wk, const double* eta, const double* xr,
                         double gamma, double dl, double* local) {
    const double d = 0.5 * (eta[a] - eta[b]);
    const PairProb pr = pair_prob(d, gamma, dl, est);
    const double ak = 0.5 * ((double)wk - (pr.pw - pr.pl));
    const double* xa = xr + (size_t)a * p;
    const double* xb = xr + (size_t)b * p;
    double* ib = IF.data() + (size_t)b * np;
    for (int r = 0; r < p; ++r) {
      const double v = ak * (xa[r] - xb[r]);
      local[r] += v;
      ib[r]   += v;
    }
    if (est) {
      const double v = (wk == 0 ? 1.0 : 0.0) - pr.pt;
      local[p] += v;
      ib[p]   += v;
    }
  }

  // Add the accumulated row-a buffer into IF.row(a).
  inline void flush_row(int a, const double* local) {
    double* ia = IF.data() + (size_t)a * np;
    for (int r = 0; r < np; ++r) ia[r] += local[r];
  }

  void merge(const IFAcc& o) { IF += o.IF; }
};

// Sandwich assembly: V = A^{-1} B A^{-1} / N, A = Hn / n_pairs, B from IF.
List vcov_finish(const RowMat& IF, const MatrixXd& Hn, double nn, int N) {
  const double scale = (double)N / nn;
  const MatrixXd A = Hn / nn;
  const MatrixXd B = (scale * scale) * (IF.transpose() * IF) / (double)N;
  auto inv_or_pinv = [](const MatrixXd& M, bool& singular) -> MatrixXd {
    Eigen::FullPivLU<MatrixXd> lu(M);
    if (lu.isInvertible()) return lu.inverse();
    singular = true;
    Eigen::CompleteOrthogonalDecomposition<MatrixXd> cod(M);
    return cod.pseudoInverse();
  };
  bool singular = false, s2 = false;
  const MatrixXd Ainv = inv_or_pinv(A, singular);
  const MatrixXd V = Ainv * B * Ainv.transpose() / (double)N;
  const MatrixXd Vnaive = inv_or_pinv(Hn, s2);
  return List::create(_["V_sandwich"] = V, _["V_model"] = Vnaive,
                      _["A"] = A, _["B"] = B, _["singular"] = singular || s2);
}

typedef std::function<void(const VectorXd&, double, bool, double&, VectorXd&, MatrixXd&)> EvalFn;

// Levenberg-Marquardt maximisation; `ev` evaluates (loglik, gradient, -Hessian).
List lm_driver(const EvalFn& ev, const MatrixXd& X, bool est_delta,
               const VectorXd& beta_init, double log_delta_init,
               int max_iter, double tol, double lambda0, bool verbose) {
  const int p = (int)X.cols();
  VectorXd beta = beta_init;
  double gamma = est_delta ? log_delta_init : 0.0;
  double lambda = lambda0;
  bool converged = false, stalled = false;
  double ll = 0.0, grad_norm = NA_REAL;
  int iter = 0;
  VectorXd g; MatrixXd Hn;

  // Diagnostics ----------------------------------------------------------
  std::vector<int> trials_per_iter;   // inner-loop trials per outer iter
  int n_full_evals = 0;               // gradient + Hessian
  int n_grad_evals = 0;               // gradient only (line-search trials)
  // ----------------------------------------------------------------------

  for (iter = 1; iter <= max_iter; ++iter) {
    ev(beta, gamma, true, ll, g, Hn);
    ++n_full_evals;
    grad_norm = g.cwiseAbs().maxCoeff() / (1.0 + std::fabs(ll));
    if (verbose && (iter <= 3 || iter % 10 == 0))
      Rcout << "  iter " << iter << ": loglik = " << ll
            << "  |g| = " << grad_norm << "  lambda = " << lambda << "\n";
    if (grad_norm < tol) { converged = true; break; }

    bool accepted = false;
    double ll_new = ll;
    int tr_used = 0;

    for (int tr = 0; tr < 10; ++tr) {
      ++tr_used;
      MatrixXd A = Hn;
      A.diagonal().array() += lambda;
      Eigen::LDLT<MatrixXd> ldlt(A);
      if (ldlt.info() != Eigen::Success) { lambda *= 5.0; continue; }
      const VectorXd step = ldlt.solve(g);
      if (!step.allFinite()) { lambda *= 5.0; continue; }
      const VectorXd beta_try = beta + step.head(p);
      const double gam_try = est_delta ? gamma + step[p] : 0.0;
      const VectorXd eta_try = X * beta_try;
      if (!beta_try.allFinite() || !std::isfinite(gam_try) ||
          0.5 * (eta_try.maxCoeff() - eta_try.minCoeff()) > 100.0) {
        lambda *= 5.0; continue;
      }
      double ll_try; VectorXd g_dummy; MatrixXd H_dummy;
      ev(beta_try, gam_try, false, ll_try, g_dummy, H_dummy);
      ++n_grad_evals;

      // Predicted reduction from the LM quadratic model:
      //   (H + lambda I) step = g  =>  pred = 0.5 * step' (g + lambda*step)
      const double pred = 0.5 * step.dot(g + lambda * step);
      const double rho  = (pred > 0.0)
                          ? (ll_try - ll) / pred
                          : (ll_try > ll ? 1.0 : -1.0);

      if (std::isfinite(ll_try) && ll_try > ll) {
        beta = beta_try; gamma = gam_try; ll_new = ll_try;
        // Nielsen gain-ratio shrink.
        const double shrink = 1.0 - std::pow(2.0 * rho - 1.0, 3.0);
        lambda *= std::max(shrink, 1.0 / 3.0);
        lambda  = std::max(lambda, 1e-12);
        accepted = true;
        break;
      }
      lambda *= 5.0;   // gentler than 10
    }
    trials_per_iter.push_back(tr_used);

    if (!accepted) { stalled = true; break; }
    if ((ll_new - ll) / (1.0 + std::fabs(ll)) < tol && grad_norm < 1e-6) {
      ll = ll_new; converged = true; break;
    }
  }
  if (iter > max_iter) iter = max_iter;

  ev(beta, gamma, true, ll, g, Hn);
  ++n_full_evals;
  grad_norm = g.cwiseAbs().maxCoeff() / (1.0 + std::fabs(ll));
  if (grad_norm < tol || (stalled && grad_norm < 1e-7)) converged = true;

  return List::create(
    _["beta"]          = beta,
    _["log_delta"]     = est_delta ? gamma : NA_REAL,
    _["loglik"]        = ll,
    _["converged"]     = converged,
    _["n_iter"]        = iter,
    _["grad_norm"]     = grad_norm,
    _["lambda_end"]    = lambda,
    _["neg_hessian"]   = Hn,
    _["trials_per_iter"] = wrap(trials_per_iter),
    _["n_full_evals"]  = n_full_evals,
    _["n_grad_evals"]  = n_grad_evals);
}

// ---- stored-pair evaluators -------------------------------------------------

void eval_stored(const RowMat& Xr, const MatrixXd& X,
                 const int* ii, const int* jj, const int* w, R_xlen_t n,
                 const VectorXd& beta, double gamma, bool est, bool need_hess,
                 double& ll, VectorXd& g, MatrixXd& Hn) {
  const int N = (int)Xr.rows(), p = (int)Xr.cols();
  const VectorXd eta = X * beta;
  const double dl = est ? std::exp(gamma) : 0.0;
  EvalAcc A; A.init(N, p, est, need_hess);
  for (R_xlen_t k = 0; k < n; ++k)
    A.add(ii[k] - 1, jj[k] - 1, w[k], eta.data(), Xr.data(), gamma, dl);
  finish_eval(A, X, est, ll, g, Hn);
}

}  // namespace

// Log-likelihood, gradient and negative Hessian (exposed for testing).
// [[Rcpp::export]]
List cpp_btd_eval(const Eigen::Map<Eigen::MatrixXd> X,
                  IntegerVector id_i, IntegerVector id_j, IntegerVector w,
                  Eigen::VectorXd beta, double log_delta, bool est_delta,
                  bool need_hess) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  double ll; VectorXd g; MatrixXd Hn;
  eval_stored(Xr, Xm, id_i.begin(), id_j.begin(), w.begin(), (R_xlen_t)id_i.size(),
              beta, log_delta, est_delta, need_hess, ll, g, Hn);
  return List::create(_["loglik"] = ll, _["grad"] = g, _["neg_hessian"] = Hn);
}

// [[Rcpp::export]]
List cpp_fit_btd(const Eigen::Map<Eigen::MatrixXd> X,
                 IntegerVector id_i, IntegerVector id_j, IntegerVector w,
                 bool est_delta, Eigen::VectorXd beta_init, double log_delta_init,
                 int max_iter = 200, double tol = 1e-10, double lambda0 = 1e-6,
                 bool verbose = false) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const int* ii = id_i.begin(); const int* jj = id_j.begin(); const int* ww = w.begin();
  const R_xlen_t n = (R_xlen_t)id_i.size();
  EvalFn ev = [&](const VectorXd& b, double gam, bool hess, double& ll,
                  VectorXd& g, MatrixXd& Hn) {
    eval_stored(Xr, Xm, ii, jj, ww, n, b, gam, est_delta, hess, ll, g, Hn);
  };
  return lm_driver(ev, Xm, est_delta, beta_init, log_delta_init, max_iter, tol,
                   lambda0, verbose);
}

// Sandwich variance (U-statistic projection) from stored pairs.
// IF_i = (N / n_pairs) * sum of scores of pairs containing i, which equals
// (2/(N-1)) * sum_j s_ij when all N(N-1)/2 pairs are used.
// [[Rcpp::export]]
List cpp_vcov_btd(const Eigen::Map<Eigen::MatrixXd> X,
                  IntegerVector id_i, IntegerVector id_j, IntegerVector w,
                  Eigen::VectorXd beta, double log_delta, bool est_delta) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const int N = (int)X.rows(), p = (int)X.cols();
  const int* ii = id_i.begin(); const int* jj = id_j.begin(); const int* ww = w.begin();
  const R_xlen_t n = (R_xlen_t)id_i.size();

  double ll; VectorXd g; MatrixXd Hn;
  eval_stored(Xr, Xm, ii, jj, ww, n, beta, log_delta, est_delta, true, ll, g, Hn);

  const VectorXd eta = Xm * beta;
  const double dl = est_delta ? std::exp(log_delta) : 0.0;
  IFAcc A; A.init(N, p, est_delta);
  for (R_xlen_t k = 0; k < n; ++k)
    A.add(ii[k] - 1, jj[k] - 1, ww[k], eta.data(), Xr.data(), log_delta, dl);
  return vcov_finish(A.IF, Hn, (double)n, N);
}


// ======================================================================
// 4. Cross-group win counts (for win statistics)
// ======================================================================

// For each subject, the number of wins, losses and ties against subjects of the
// *other* group. `group` is 0/1 per subject. Returns an N x 3 matrix.
// [[Rcpp::export]]
NumericMatrix cpp_cross_sums(IntegerVector id_i, IntegerVector id_j,
                             IntegerVector w, IntegerVector group) {
  const int N = group.size();
  NumericMatrix out(N, 3);
  const R_xlen_t n = id_i.size();
  for (R_xlen_t k = 0; k < n; ++k) {
    const int a = id_i[k] - 1, b = id_j[k] - 1;
    if (group[a] == group[b]) continue;
    if (w[k] == 1)       { out(a, 0) += 1; out(b, 1) += 1; }
    else if (w[k] == -1) { out(a, 1) += 1; out(b, 0) += 1; }
    else                 { out(a, 2) += 1; out(b, 2) += 1; }
  }
  return out;
}


// ======================================================================
// 4. STREAMING engine: pairs generated and scored on the fly, never stored
// ======================================================================

namespace {

// Score one pair under the whole (built-in) hierarchy.
// Returns W; `lev` = 0-based index of the resolving level (K if none); `uf` = flag.
inline int score_pair(const std::vector<LevelSpec>& L, const double* V, R_xlen_t N,
                      int a, int b, int& lev, int& uf) {
  const int K = (int)L.size();
  uf = 0; lev = K;
  for (int l = 0; l < K; ++l) {
    bool indet = false;
    const int r = eval_level(L[l], V, N, a, b, indet);
    if (r != 0) { lev = l; return r; }
    if (indet) uf = 1;
  }
  return 0;
}

// Row boundaries of the chunks: row i owns the N-1-i pairs (i, j>i); chunks hold
// (nearly) equal numbers of pairs. Depends only on (N, nchunk).
std::vector<int> make_chunks(int N, int nchunk) {
  std::vector<int> b(1, 0);
  const double total = 0.5 * (double)N * ((double)N - 1.0);
  double cum = 0.0;
  int k = 1;
  for (int i = 0; i < N - 1; ++i) {
    cum += (double)(N - 1 - i);
    while (k < nchunk && cum >= total * (double)k / (double)nchunk) {
      if (i + 1 > b.back()) b.push_back(i + 1);
      ++k;
    }
  }
  if (b.back() != N - 1) b.push_back(N - 1);
  return b;
}

// Run f(acc[chunk], i, j) over all pairs i < j; chunks may run in parallel.
template <class Acc, class F>
void stream_pairs(int N, const std::vector<int>& bd, int threads,
                  std::vector<Acc>& acc, F f) {
  const int nc = (int)bd.size() - 1;
  const int nthr = std::max(1, threads);
  (void)nthr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthr) schedule(dynamic, 1)
#endif
  for (int ch = 0; ch < nc; ++ch) {
    Acc& a = acc[ch];
    for (int i = bd[ch]; i < bd[ch + 1]; ++i)
      for (int j = i + 1; j < N; ++j) f(a, i, j);
  }
}

// Tiled enumeration of pairs i < j for IF accumulation.
//
// The pair range is blocked into T x T tiles. Within a tile, the rows of IF
// touched are exactly [i0, i1) and [j0, j1), so the working set is bounded and
// stays in L2 for a suitable T (T = 512 works well for np <= 8; peak footprint
// is about 2 * T * np doubles).
//
// Row i's contribution is accumulated into a per-chunk buffer and flushed into
// IF.row(i) once per (i, j-tile), so IF.row(i) is written ~(N/T) times per row
// instead of (N - i) times. Writes to IF.row(j) remain sequential in j.
//
// `score(i, j)` returns W in {-1, 0, 1}.  The lambda captures the engine-
// specific state (levels/V for stream, codes for compact).
template <class ScoreFn>
void stream_pairs_if_tiled(int N, const std::vector<int>& bd, int threads,
                           std::vector<IFAcc>& acc, int T,
                           const double* ep, const double* xp,
                           double log_delta, double dl,
                           ScoreFn score) {
  const int nc = (int)bd.size() - 1;
  const int nthr = std::max(1, threads);
  (void)nthr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthr) schedule(dynamic, 1)
#endif
  for (int ch = 0; ch < nc; ++ch) {
    IFAcc& a = acc[ch];
    const int np = a.np;
    std::vector<double> local(np);
    const int i_lo = bd[ch], i_hi = bd[ch + 1];
    for (int i0 = i_lo; i0 < i_hi; i0 += T) {
      const int i1 = std::min(i0 + T, i_hi);
      for (int j0 = i0; j0 < N; j0 += T) {
        const int j1 = std::min(j0 + T, N);
        for (int i = i0; i < i1; ++i) {
          std::fill(local.begin(), local.end(), 0.0);
          const int jmin = (j0 <= i) ? i + 1 : j0;
          for (int j = jmin; j < j1; ++j) {
            const int wk = score(i, j);
            a.add_to_row(i, j, wk, ep, xp, log_delta, dl, local.data());
          }
          a.flush_row(i, local.data());
        }
      }
    }
  }
}

void eval_stream(const RowMat& Xr, const MatrixXd& X,
                 const std::vector<LevelSpec>& L, const double* V,
                 int nchunk, int threads,
                 const VectorXd& beta, double gamma, bool est, bool need_hess,
                 double& ll, VectorXd& g, MatrixXd& Hn) {
  const int N = (int)Xr.rows(), p = (int)Xr.cols();
  const VectorXd eta = X * beta;
  const double dl = est ? std::exp(gamma) : 0.0;
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<EvalAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.init(N, p, est, need_hess);
  const double* ep = eta.data(); const double* xp = Xr.data();
  stream_pairs(N, bd, threads, acc, [&](EvalAcc& a, int i, int j) {
    int lev, uf;
    const int w = score_pair(L, V, N, i, j, lev, uf);
    a.add(i, j, w, ep, xp, gamma, dl);
  });
  for (size_t c = 1; c < acc.size(); ++c) acc[0].merge(acc[c]);
  finish_eval(acc[0], X, est, ll, g, Hn);
}

}  // namespace

// Pair counts per resolving level, without storing pairs.
// Returns:
//   counts    length K+1; counts[k] = pairs resolved at level k (0-based),
//             counts[K] = pairs unresolved at every level (tied in the fit)
//   uninf     pairs that hit at least one censored/informative-flagged
//             comparison at any level, whether or not later resolved
//   cens_tie  pairs unresolved at every level AND flagged by censoring,
//             i.e. ties that are attributable to censoring rather than to
//             genuine equivalence
//   n_pairs   total number of pairs
// [[Rcpp::export]]
List cpp_stream_counts(List levels, NumericMatrix V, int nchunk = 32, int threads = 1) {
  const std::vector<LevelSpec> L = parse_levels(levels);
  const int K = (int)L.size();
  const int N = V.nrow();
  const double* Vp = V.begin();
  struct CAcc {
    char pad_head_[128];
    std::vector<double> cnt;   // length K+1
    double unf = 0.0;          // any uf at any level
    double cens_tie = 0.0;     // lev == K && uf == 1
    char pad_tail_[128];
  };
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<CAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.cnt.assign(K + 1, 0.0);
  stream_pairs(N, bd, threads, acc, [&](CAcc& a, int i, int j) {
    int lev, uf;
    score_pair(L, Vp, N, i, j, lev, uf);
    a.cnt[lev] += 1.0;
    a.unf += uf;
    if (lev == K && uf) a.cens_tie += 1.0;
  });
  NumericVector cnt(K + 1);
  double unf = 0.0, cens_tie = 0.0;
  for (auto& a : acc) {
    for (int k = 0; k <= K; ++k) cnt[k] += a.cnt[k];
    unf      += a.unf;
    cens_tie += a.cens_tie;
  }
  return List::create(_["counts"]    = cnt,
                      _["uninf"]     = unf,
                      _["cens_tie"]  = cens_tie,
                      _["n_pairs"]   = 0.5 * (double)N * ((double)N - 1.0));
}

// Stored-free Levenberg-Marquardt fit.
// [[Rcpp::export]]
List cpp_stream_fit(const Eigen::Map<Eigen::MatrixXd> X, List levels, NumericMatrix V,
                    bool est_delta, Eigen::VectorXd beta_init, double log_delta_init,
                    int max_iter = 200, double tol = 1e-10, double lambda0 = 1e-6,
                    bool verbose = false, int nchunk = 32, int threads = 1) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const std::vector<LevelSpec> L = parse_levels(levels);
  const double* Vp = V.begin();
  EvalFn ev = [&](const VectorXd& b, double gam, bool hess, double& ll,
                  VectorXd& g, MatrixXd& Hn) {
    eval_stream(Xr, Xm, L, Vp, nchunk, threads, b, gam, est_delta, hess, ll, g, Hn);
  };
  return lm_driver(ev, Xm, est_delta, beta_init, log_delta_init, max_iter, tol,
                   lambda0, verbose);
}

// Stored-free sandwich variance (tiled, row-local IF accumulation).
// [[Rcpp::export]]
List cpp_stream_vcov(const Eigen::Map<Eigen::MatrixXd> X, List levels, NumericMatrix V,
                     Eigen::VectorXd beta, double log_delta, bool est_delta,
                     int nchunk = 32, int threads = 1, int tile = 512) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const int N = (int)X.rows(), p = (int)X.cols();
  const std::vector<LevelSpec> L = parse_levels(levels);
  const double* Vp = V.begin();

  // Hessian: exactly one pass.
  double ll; VectorXd g; MatrixXd Hn;
  eval_stream(Xr, Xm, L, Vp, nchunk, threads, beta, log_delta, est_delta,
              true, ll, g, Hn);

  const VectorXd eta = Xm * beta;
  const double dl = est_delta ? std::exp(log_delta) : 0.0;
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<IFAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.init(N, p, est_delta);
  const double* ep = eta.data(); const double* xp = Xr.data();

  auto score = [&](int i, int j) -> int {
    int lev, uf;
    return score_pair(L, Vp, N, i, j, lev, uf);
  };
  stream_pairs_if_tiled(N, bd, threads, acc, tile, ep, xp, log_delta, dl, score);

  for (size_t c = 1; c < acc.size(); ++c) acc[0].merge(acc[c]);
  return vcov_finish(acc[0].IF, Hn, 0.5 * (double)N * ((double)N - 1.0), N);
}

// Cross-group win counts per subject and per resolving level, for win statistics.
// group is 0/1; "treatment" is group 1. by_level columns: treatment wins,
// treatment losses, ties (counts of cross-group pairs by resolving level).
// [[Rcpp::export]]
List cpp_stream_cross(List levels, NumericMatrix V, IntegerVector group,
                      int nchunk = 32, int threads = 1) {
  const std::vector<LevelSpec> L = parse_levels(levels);
  const int K = (int)L.size();
  const int N = V.nrow();
  const double* Vp = V.begin();
  struct XAcc { std::vector<double> cs; std::vector<double> bl; };
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<XAcc> acc(bd.size() - 1);
  for (auto& a : acc) { a.cs.assign((size_t)N * 3, 0.0); a.bl.assign((size_t)(K + 1) * 3, 0.0); }
  const int* gp = group.begin();
  stream_pairs(N, bd, threads, acc, [&](XAcc& a, int i, int j) {
    if (gp[i] == gp[j]) return;
    int lev, uf;
    const int w = score_pair(L, Vp, N, i, j, lev, uf);
    // subject-level counts, wins/losses/ties from each subject's perspective
    if (w == 1)       { a.cs[i] += 1; a.cs[N + j] += 1; }
    else if (w == -1) { a.cs[N + i] += 1; a.cs[j] += 1; }
    else              { a.cs[2 * (size_t)N + i] += 1; a.cs[2 * (size_t)N + j] += 1; }
    const bool trt_wins = (gp[i] == 1) ? (w == 1) : (w == -1);
    const bool trt_loss = (gp[i] == 1) ? (w == -1) : (w == 1);
    if (trt_wins) a.bl[lev] += 1; else if (trt_loss) a.bl[(K + 1) + lev] += 1;
    else a.bl[2 * (K + 1) + lev] += 1;
  });
  NumericMatrix cs(N, 3), bl(K + 1, 3);
  for (auto& a : acc) {
    for (int s = 0; s < 3; ++s) {
      for (int i = 0; i < N; ++i) cs(i, s) += a.cs[(size_t)s * N + i];
      for (int k = 0; k <= K; ++k) bl(k, s) += a.bl[(size_t)s * (K + 1) + k];
    }
  }
  return List::create(_["cross_sums"] = cs, _["by_level"] = bl);
}

// Calibration sums by resolving level and by fine bins of the latent difference d.
// [[Rcpp::export]]
List cpp_stream_calib(List levels, NumericMatrix V, NumericVector eta, double delta,
                      int nfine = 512, int nchunk = 32, int threads = 1) {
  const std::vector<LevelSpec> L = parse_levels(levels);
  const int K = (int)L.size();
  const int N = V.nrow();
  const double* Vp = V.begin();
  double emin = eta[0], emax = eta[0];
  for (int i = 1; i < N; ++i) { emin = std::min(emin, eta[i]); emax = std::max(emax, eta[i]); }
  const double dmax = 0.5 * (emax - emin);
  const double width = dmax > 0 ? 2.0 * dmax / nfine : 1.0;
  struct KAcc { std::vector<double> lv, fb; };      // lv: K x 3, fb: nfine x 7
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<KAcc> acc(bd.size() - 1);
  for (auto& a : acc) { a.lv.assign((size_t)K * 3, 0.0); a.fb.assign((size_t)nfine * 7, 0.0); }
  const bool est = delta > 0;
  const double* ep = &eta[0];
  stream_pairs(N, bd, threads, acc, [&](KAcc& a, int i, int j) {
    int lev, uf;
    const int w = score_pair(L, Vp, N, i, j, lev, uf);
    const double d = 0.5 * (ep[i] - ep[j]);
    const PairProb pr = pair_prob(d, est ? std::log(delta) : 0.0, delta, est);
    if (lev < K) {
      a.lv[lev * 3 + 0] += 1.0;
      a.lv[lev * 3 + 1] += (w == 1) ? 1.0 : 0.0;
      a.lv[lev * 3 + 2] += pr.pw / (pr.pw + pr.pl);
    }
    int f = (int)std::floor((d + dmax) / width);
    if (f < 0) f = 0; if (f >= nfine) f = nfine - 1;
    double* b = &a.fb[(size_t)f * 7];
    b[0] += 1.0; b[1] += (w == 1); b[2] += (w == 0); b[3] += (w == -1);
    b[4] += pr.pw; b[5] += pr.pt; b[6] += pr.pl;
  });
  NumericMatrix lv(K, 3), fb(nfine, 7);
  for (auto& a : acc) {
    for (int k = 0; k < K; ++k) for (int s = 0; s < 3; ++s) lv(k, s) += a.lv[(size_t)k * 3 + s];
    for (int f = 0; f < nfine; ++f) for (int s = 0; s < 7; ++s) fb(f, s) += a.fb[(size_t)f * 7 + s];
  }
  return List::create(_["by_level"] = lv, _["fine"] = fb, _["dmax"] = dmax, _["width"] = width);
}


// ======================================================================
// 5. COMPACT engine: one byte per pair, scored once
// ======================================================================
//
// code = level * 3 + (W + 1), with level in 0..K (K = unresolved), so K <= 84.
// Pair (i, j), i < j, has index row_start(i) + (j - i - 1) in the enumeration used
// everywhere else, which is what lets us drop the index vectors.

namespace {

inline R_xlen_t row_start(int N, int i) {
  return (R_xlen_t)i * (N - 1) - (R_xlen_t)i * (i - 1) / 2;
}

// Pair index in the i-outer, j-inner enumeration (0-based i, j with i < j).
inline R_xlen_t pair_index(int N, int i, int j) {
  return row_start(N, i) + (j - i - 1);
}

// Like stream_pairs, but also passes the running pair index k.
template <class Acc, class F>
void stream_pairs_k(int N, const std::vector<int>& bd, int threads,
                    std::vector<Acc>& acc, F f) {
  const int nc = (int)bd.size() - 1;
  const int nthr = std::max(1, threads);
  (void)nthr;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthr) schedule(dynamic, 1)
#endif
  for (int ch = 0; ch < nc; ++ch) {
    Acc& a = acc[ch];
    R_xlen_t k = row_start(N, bd[ch]);
    for (int i = bd[ch]; i < bd[ch + 1]; ++i)
      for (int j = i + 1; j < N; ++j) f(a, i, j, k++);
  }
}

void eval_compact(const RowMat& Xr, const MatrixXd& X, const unsigned char* codes,
                  int nchunk, int threads,
                  const VectorXd& beta, double gamma, bool est, bool need_hess,
                  double& ll, VectorXd& g, MatrixXd& Hn) {
  const int N = (int)Xr.rows(), p = (int)Xr.cols();
  const VectorXd eta = X * beta;
  const double dl = est ? std::exp(gamma) : 0.0;
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<EvalAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.init(N, p, est, need_hess);
  const double* ep = eta.data(); const double* xp = Xr.data();
  stream_pairs_k(N, bd, threads, acc, [&](EvalAcc& a, int i, int j, R_xlen_t k) {
    a.add(i, j, (int)(codes[k] % 3) - 1, ep, xp, gamma, dl);
  });
  for (size_t c = 1; c < acc.size(); ++c) acc[0].merge(acc[c]);
  finish_eval(acc[0], X, est, ll, g, Hn);
}

}  

// Score every pair once and keep one byte per pair. Also returns the number of
// pairs resolved at each level (length K+1), the number of pairs flagged by
// censoring at any level, the number of pairs unresolved because of censoring
// only, and the number of pairs.
// [[Rcpp::export]]
List cpp_compact_build(List levels, NumericMatrix V, int nchunk = 32, int threads = 1) {
  const std::vector<LevelSpec> L = parse_levels(levels);
  const int K = (int)L.size();
  if (K > 84) stop("compact storage supports at most 84 comparisons in a hierarchy");
  const int N = V.nrow();
  const double* Vp = V.begin();
  const double nd = 0.5 * (double)N * ((double)N - 1.0);
  RawVector codes((R_xlen_t)nd);
  unsigned char* cp = (unsigned char*)RAW(codes);
  struct BAcc {
    char pad_head_[128];
    std::vector<double> cnt;   // length K+1
    double unf = 0.0;          // any uf at any level
    double cens_tie = 0.0;     // lev == K && uf == 1
    char pad_tail_[128];
  };
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<BAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.cnt.assign(K + 1, 0.0);
  stream_pairs_k(N, bd, threads, acc, [&](BAcc& a, int i, int j, R_xlen_t k) {
    int lev, uf;
    const int w = score_pair(L, Vp, N, i, j, lev, uf);
    cp[k] = (unsigned char)(lev * 3 + (w + 1));
    a.cnt[lev] += 1.0;
    a.unf += uf;
    if (lev == K && uf) a.cens_tie += 1.0;
  });
  NumericVector cnt(K + 1);
  double unf = 0.0, cens_tie = 0.0;
  for (auto& a : acc) {
    for (int k = 0; k <= K; ++k) cnt[k] += a.cnt[k];
    unf      += a.unf;
    cens_tie += a.cens_tie;
  }
  return List::create(_["codes"]     = codes,
                      _["counts"]    = cnt,
                      _["uninf"]     = unf,
                      _["cens_tie"]  = cens_tie);
}

// Levenberg-Marquardt fit reading the compact codes.
// [[Rcpp::export]]
List cpp_compact_fit(const Eigen::Map<Eigen::MatrixXd> X, RawVector codes,
                     bool est_delta, Eigen::VectorXd beta_init, double log_delta_init,
                     int max_iter = 200, double tol = 1e-10, double lambda0 = 1e-6,
                     bool verbose = false, int nchunk = 32, int threads = 1) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const unsigned char* cp = (const unsigned char*)RAW(codes);
  EvalFn ev = [&](const VectorXd& b, double gam, bool hess, double& ll,
                  VectorXd& g, MatrixXd& Hn) {
    eval_compact(Xr, Xm, cp, nchunk, threads, b, gam, est_delta, hess, ll, g, Hn);
  };
  return lm_driver(ev, Xm, est_delta, beta_init, log_delta_init, max_iter, tol,
                   lambda0, verbose);
}

// Sandwich variance reading the compact codes (tiled, row-local IF accumulation).
// [[Rcpp::export]]
List cpp_compact_vcov(const Eigen::Map<Eigen::MatrixXd> X, RawVector codes,
                      Eigen::VectorXd beta, double log_delta, bool est_delta,
                      int nchunk = 32, int threads = 1, int tile = 512) {
  const MatrixXd Xm = X;
  const RowMat Xr = X;
  const int N = (int)X.rows(), p = (int)X.cols();
  const unsigned char* cp = (const unsigned char*)RAW(codes);

  // Hessian: exactly one pass.
  double ll; VectorXd g; MatrixXd Hn;
  eval_compact(Xr, Xm, cp, nchunk, threads, beta, log_delta, est_delta,
               true, ll, g, Hn);

  const VectorXd eta = Xm * beta;
  const double dl = est_delta ? std::exp(log_delta) : 0.0;
  const std::vector<int> bd = make_chunks(N, nchunk);
  std::vector<IFAcc> acc(bd.size() - 1);
  for (auto& a : acc) a.init(N, p, est_delta);
  const double* ep = eta.data(); const double* xp = Xr.data();

  auto score = [&](int i, int j) -> int {
    return (int)(cp[pair_index(N, i, j)] % 3) - 1;
  };
  stream_pairs_if_tiled(N, bd, threads, acc, tile, ep, xp, log_delta, dl, score);

  for (size_t c = 1; c < acc.size(); ++c) acc[0].merge(acc[c]);
  return vcov_finish(acc[0].IF, Hn, 0.5 * (double)N * ((double)N - 1.0), N);
}
