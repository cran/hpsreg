## ------------------------------------------------------------------
## Davidson (Bradley-Terry with ties) regression on scored pairs
## ------------------------------------------------------------------

#' Control parameters for the Bradley-Terry-Davidson fit
#'
#' @param max_iter Maximum Levenberg-Marquardt iterations.
#' @param tol Convergence tolerance on the scaled gradient,
#'   `max|g| / (1 + |loglik|)`.
#' @param lambda0 Initial damping.
#' @param verbose Print iteration progress (C++ engine only).
#' @param threads Number of OpenMP threads for the streaming engine. `NULL`
#'   (default) uses `getOption("hpsreg.threads", 1)`. Ignored if the package was
#'   built without OpenMP (see [openmp_threads()]). Results do not depend on it.
#' @return A list.
#' @export
btd_control <- function(max_iter = 200L, tol = 1e-06, lambda0 = 1e-6,
                        verbose = FALSE, threads = NULL) {
  list(max_iter = as.integer(max_iter), tol = tol, lambda0 = lambda0,
       verbose = isTRUE(verbose), threads = threads)
}

#' Covariate-adjusted analysis of hierarchical pairwise outcomes
#'
#' Fits the Bradley-Terry-Davidson model to the win / loss / tie outcomes produced by
#' [score_pairs()]. For a pair of subjects with covariate difference
#' \eqn{x_i - x_j} and \eqn{d = (x_i - x_j)^\top \beta / 2},
#' \deqn{P(i\ \mathrm{wins}) = e^{d}/q,\quad P(j\ \mathrm{wins}) = e^{-d}/q,
#'   \quad P(\mathrm{tie}) = \delta/q,}
#' with \eqn{q = e^{d} + e^{-d} + \delta}. A coefficient is the log of the
#' multiplicative effect of a one-unit increase of the covariate on a subject's
#' "worth"; positive values mean a better position in the hierarchy. If no pair
#' is tied the model reduces to Bradley-Terry (\eqn{\delta = 0}).
#'
#' Estimation is by Levenberg-Marquardt on standardised covariates (C++ engine,
#' or a pure R reference implementation with `engine = "r"`). Standard errors are
#' sandwich (robust) estimates based on the projection of the pairwise
#' estimating equation onto subjects, which accounts for the dependence between
#' pairs that share a subject (Hoeffding / U-statistic variance). They assume the
#' complete set of N(N-1)/2 pairs; for a sampled or user-supplied subset the
#' normalisation uses the average number of pairs per subject and should be
#' treated as approximate. The model-based inverse information, which ignores
#' that dependence, is returned as `vcov_model` for reference only.
#'
#' @param formula One-sided formula (`~ x1 + x2`) for the covariates. The
#'   intercept is not identifiable and is dropped.
#' @param data Data frame with one row per subject, the same rows used in
#'   [score_pairs()].
#' @param pairs An `hb_pairs` object from [score_pairs()]. If `NULL`, supply
#'   `hierarchy` and the pairs are scored first.
#' @param hierarchy A [hierarchy()], used only if `pairs` is `NULL`.
#' @param scale Standardise covariates internally (recommended).
#' @param engine `"cpp"` (default) or `"r"` (slow reference implementation).
#' @param storage How pairs are handled when `hierarchy` is given.
#'   * `"stored"`: score and keep all pairs (about 20 bytes per pair, more at
#'     peak), available afterwards as `$pairs`.
#'   * `"compact"`: score every pair once and keep one byte per pair (outcome and
#'     resolving level; indices are implicit). Fit and variance passes read this
#'     array instead of re-scoring. About 50 MB at N = 10^4.
#'   * `"stream"`: never store anything proportional to the number of pairs; pairs
#'     are re-scored in every pass. Memory is O(N p), at the cost of extra compute.
#'   * `"auto"` (default): `"stored"` up to `getOption("hpsreg.stream_min_pairs",
#'     2.5e6)` pairs (about N = 2200), then `"compact"` up to
#'     `getOption("hpsreg.compact_max_pairs", 1e9)` pairs (about N = 44,000),
#'     then `"stream"`; for hierarchies with [lvl_custom()] levels, always `"stored"`.
#'
#'   `"compact"` and `"stream"` require `engine = "cpp"` and a `hierarchy` (not
#'   `pairs`), support at most 84 comparisons for `"compact"`, and give a fit
#'   without a `$pairs` component. Estimates and standard errors agree with the
#'   stored engine to numerical precision, and `"compact"` and `"stream"` agree
#'   bitwise.
#' @param control See [btd_control()].
#'
#' @return An object of class `btd` with `coefficients`, `log_delta`, `delta`,
#'   `vcov` (sandwich, on the original covariate scale, including `log_delta`),
#'   `vcov_model`, `loglik`, convergence information, the latent scores `eta`
#'   and the pairs used.
#'
#'   The object supports the usual generics: `summary()` (coefficient table
#'   with z-tests), `coef()`, `vcov()` (`model = TRUE` for the model-based
#'   information), `confint()` (Wald intervals), `logLik()`, `nobs()`
#'   (number of pairs; the subject count is `fit$n_subjects`), `print()`, and
#'   `predict()` (latent score for the training rows, or for a `newdata` frame). 
#'  
#'   Pair-level predicted probabilities are provided by [btd_pair_probs()]; 
#'   calibration tables by [btd_calibration()].
#'   
#' @section Two workflows:
#' `btd()` needs the scored pairs of a hierarchy. You can supply them or let
#' `btd()` compute them:
#'
#' 1. **Pre-computed pairs** (`pairs = pw`). Score once with [score_pairs()],
#'    inspect with [resolution()] or [win_stats()], then fit several covariate
#'    models on the same pairs without rescoring.
#' 2. **Inline scoring** (`hierarchy = h`). `btd()` calls [score_pairs()] for
#'    you; the pairs are kept in the fitted object (`fit$pairs`) for later use.
#'
#' In both cases `data` must be the same data frame (same rows, same order)
#' that is scored.
#'
#' @examples
#' ## From scratch: define a hierarchy from data-frame columns, score, fit.
#' dat <- simulate_hps(N = 150, seed = 11)
#'
#' h <- hierarchy(
#'   lvl_tte("death_day", better = "later", name = "survival"),
#'   lvl_tte("relapse_day", better = "later", name = "relapse"),
#'   lvl_ordinal("gvhd_grade", better = "lower", name = "GVHD grade"),
#'   lvl_binary("gvhd_resolved", favorable = 1, when = "gvhd_grade > 0",
#'              name = "GVHD resolution")
#' )
#'
#' ## Workflow 1: pre-computed pairs (reusable across models)
#' pw <- score_pairs(h, dat)
#' resolution(pw)
#' fit1 <- btd(~ donor_type + age + dri + hct_ci, data = dat, pairs = pw)
#' summary(fit1)
#' fit2 <- btd(~ donor_type + age, data = dat, pairs = pw)   # no rescoring
#'
#' ## Workflow 2: let btd() score the pairs
#' fit3 <- btd(~ donor_type + age, data = dat, hierarchy = h)
#' all.equal(coef(fit2), coef(fit3))
#'
#' ## The packaged transplant hierarchy is a shortcut for the one above
#' fit4 <- btd(~ donor_type + age, data = dat, hierarchy = hps_hierarchy())
#'
#' ## How many threads can this build use? 0 means no OpenMP, in which
#' ## case the engine runs serially regardless of the `threads` argument.
#' openmp_threads()
#'
#' \donttest{
#' ## Multiple threads: results are bitwise identical to a single-threaded
#' ## fit, so this only changes the runtime, not the estimates.
#' n <- max(1L, min(4L, openmp_threads()))
#' fit_mt <- btd(~ donor_type + age, data = dat, hierarchy = h,
#'               control = btd_control(threads = n))
#' all.equal(coef(fit_mt), coef(fit3))
#'
#' ## The same setting can be supplied through an option, which is
#' ## convenient when many calls are made in a script.
#' op <- options(hpsreg.threads = n)
#' fit_mt2 <- btd(~ donor_type + age, data = dat, hierarchy = h)
#' options(op)
#' all.equal(coef(fit_mt), coef(fit_mt2))
#' }
#' 
#' @references
#' Varadhan R. "Hierarchical Pairwise Scoring: A Regression Framework
#'   for Composite Outcomes, with Application to Hematopoietic Stem Cell
#'   Transplantation." Manuscript in preparation, 2026.
#'   
#' @export
btd <- function(formula, data, pairs = NULL, hierarchy = NULL, scale = TRUE,
                engine = c("cpp", "r"), storage = c("auto", "stored", "compact", "stream"),
                control = btd_control()) {
  cl <- match.call()
  engine <- match.arg(engine)
  storage <- match.arg(storage)
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)
  data <- as.data.frame(data)

  tt <- stats::delete.response(stats::terms(formula, data = data))
  mf <- stats::model.frame(tt, data, na.action = stats::na.fail)
  X0 <- stats::model.matrix(tt, mf)
  contr <- attr(X0, "contrasts")
  icpt <- match("(Intercept)", colnames(X0))
  X <- if (is.na(icpt)) X0 else X0[, -icpt, drop = FALSE]
  p <- ncol(X)
  if (p == 0L) stop("The model has no covariates.", call. = FALSE)
  if (nrow(X) != nrow(data)) stop("Internal error: row mismatch.", call. = FALSE)

  stream <- NULL
  if (is.null(pairs)) {
    if (is.null(hierarchy)) stop("Supply `pairs` or `hierarchy`.", call. = FALSE)
    if (!inherits(hierarchy, "hb_hierarchy"))
      stop("`hierarchy` must come from hierarchy().", call. = FALSE)
    comp <- .compile_hierarchy(hierarchy, data)
    specs <- .stream_specs(comp)
    npair <- as.numeric(nrow(X)) * (nrow(X) - 1) / 2
    if (storage %in% c("stream", "compact") && is.null(specs))
      stop("Hierarchies with lvl_custom() levels cannot use storage = \"", storage,
           "\"; use storage = \"stored\".", call. = FALSE)
    if (storage %in% c("stream", "compact") && engine == "r")
      stop("storage = \"", storage, "\" requires engine = \"cpp\".", call. = FALSE)
    mode <- storage
    if (storage == "auto") {
      mode <- if (is.null(specs) || engine != "cpp" ||
                  npair <= getOption("hpsreg.stream_min_pairs", 2.5e6)) "stored" else
              if (npair <= getOption("hpsreg.compact_max_pairs", 1e9) &&
                  length(specs) <= 84L) "compact" else "stream"
    }
    if (mode == "compact" && length(specs) > 84L)
      stop("storage = \"compact\" supports at most 84 comparisons.", call. = FALSE)
    if (mode %in% c("stream", "compact")) {
      thr <- .n_threads(control)
      if (mode == "compact") {
        cb <- cpp_compact_build(specs, comp$V, .NCHUNK, thr)
        cn <- list(counts   = cb$counts,
                   uninf    = cb$uninf,
                   cens_tie = cb$cens_tie)
        codes <- cb$codes
      } else {
        cn <- cpp_stream_counts(specs, comp$V, .NCHUNK, thr)
        codes <- NULL
      }
      stream <- list(comp = comp, specs = specs, hierarchy = hierarchy,
                     counts = cn$counts, uninf = cn$uninf,
                     cens_tie = cn$cens_tie,
                     threads = thr,
                     codes = codes, mode = mode)
    } else {
      pairs <- score_pairs(hierarchy, data)
    }
  } else {
    if (storage %in% c("stream", "compact"))
      stop("storage = \"", storage, "\" needs `hierarchy`, not pre-computed `pairs`.",
           call. = FALSE)
    if (!inherits(pairs, "hb_pairs")) stop("`pairs` must come from score_pairs().", call. = FALSE)
    if (pairs$N != nrow(X))
      stop("`pairs` was scored on ", pairs$N, " subjects but `data` has ",
           nrow(X), " rows.", call. = FALSE)
  }

  if (is.null(stream)) {
    n_pairs  <- length(pairs$W)
    n_tie    <- sum(pairs$W == 0L)
    cens_tie <- if (!is.null(pairs$uninf)) {
      sum(pairs$W == 0L & pairs$uninf == 1L)
    } else {
      NA_real_
    }
  } else {
    n_pairs  <- as.numeric(nrow(X)) * (nrow(X) - 1) / 2
    n_tie    <- stream$counts[length(stream$counts)]
    cens_tie <- stream$cens_tie
  }
  if (n_pairs - n_tie == 0) stop("Every pair is tied; nothing to estimate.", call. = FALSE)
  est_delta <- n_tie > 0
  if (!est_delta)
    message("No tied pairs: fitting the Bradley-Terry model (delta = 0).")
  
  cen <- colMeans(X)
  sc <- if (scale) apply(X, 2L, stats::sd) else rep(1, p)
  if (any(!is.finite(sc) | sc == 0))
    stop("Covariate(s) with no variation: ",
         paste(colnames(X)[!is.finite(sc) | sc == 0], collapse = ", "), call. = FALSE)
  Xs <- sweep(sweep(X, 2L, cen, "-"), 2L, sc, "/")

  if (!is.null(stream) && !is.null(stream$codes)) {
    thr <- stream$threads
    fit <- cpp_compact_fit(Xs, stream$codes, est_delta, rep(0, p), 0,
                           control$max_iter, control$tol, control$lambda0,
                           control$verbose, .NCHUNK, thr)
    vc <- cpp_compact_vcov(Xs, stream$codes, fit$beta,
                           if (est_delta) fit$log_delta else 0, est_delta, .NCHUNK, thr)
  } else if (!is.null(stream)) {
    thr <- stream$threads
    fit <- cpp_stream_fit(Xs, stream$specs, stream$comp$V, est_delta, rep(0, p), 0,
                          control$max_iter, control$tol, control$lambda0,
                          control$verbose, .NCHUNK, thr)
    vc <- cpp_stream_vcov(Xs, stream$specs, stream$comp$V, fit$beta,
                          if (est_delta) fit$log_delta else 0, est_delta, .NCHUNK, thr)
  } else if (engine == "cpp") {
    fit <- cpp_fit_btd(Xs, pairs$id_i, pairs$id_j, pairs$W, est_delta,
                       rep(0, p), 0, control$max_iter, control$tol,
                       control$lambda0, control$verbose)
    vc <- cpp_vcov_btd(Xs, pairs$id_i, pairs$id_j, pairs$W, fit$beta,
                       if (est_delta) fit$log_delta else 0, est_delta)
  } else {
    fit <- .btd_fit_r(Xs, pairs$id_i, pairs$id_j, pairs$W, est_delta, control)
    vc <- .btd_vcov_r(Xs, pairs$id_i, pairs$id_j, pairs$W, fit$beta,
                      if (est_delta) fit$log_delta else 0, est_delta)
  }
  if (!isTRUE(fit$converged))
    warning(sprintf("Fit did not converge in %d iterations (scaled |g| = %.2e).",
                    fit$n_iter, fit$grad_norm), call. = FALSE)
  if (isTRUE(vc$singular))
    warning("Information matrix is singular (collinear covariates?); ",
            "a pseudo-inverse was used for the variance.", call. = FALSE)

  scl <- c(sc, if (est_delta) 1)
  pn <- c(colnames(X), if (est_delta) "log_delta")
  back <- function(V) { V <- V / outer(scl, scl); dimnames(V) <- list(pn, pn); V }
  beta <- fit$beta / sc
  names(beta) <- colnames(X)
  eta <- drop(Xs %*% fit$beta)

  structure(list(
    coefficients = beta,
    log_delta = if (est_delta) fit$log_delta else NA_real_,
    delta = if (est_delta) exp(fit$log_delta) else 0,
    vcov = back(vc$V_sandwich),
    vcov_model = back(vc$V_model),
    loglik = fit$loglik, n_iter = fit$n_iter,
    n_full_evals    = fit$n_full_evals,
    n_grad_evals    = fit$n_grad_evals,
    trials_per_iter = fit$trials_per_iter,
    converged = isTRUE(fit$converged),
    grad_norm = fit$grad_norm, engine = engine,
    n_subjects = nrow(X), n_pairs = n_pairs, n_ties = n_tie,
    cens_tie = cens_tie,
    est_delta = est_delta,
    storage = if (is.null(stream)) "stored" else stream$mode,
    stream = if (is.null(stream)) NULL else { stream$codes <- NULL; stream },
    eta = eta, x_center = cen, x_scale = sc,
    terms = tt, xlevels = stats::.getXlevels(tt, mf), contrasts = contr,
    pairs = pairs, call = cl),
    class = "btd")
}

#' @export
coef.btd <- function(object, ...) object$coefficients

#' @export
vcov.btd <- function(object, model = FALSE, ...) {
  V <- if (model) object$vcov_model else object$vcov
  p <- length(object$coefficients)
  V[seq_len(p), seq_len(p), drop = FALSE]
}

#' @export
logLik.btd <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients) + object$est_delta,
            nobs = object$n_pairs, class = "logLik")
}

#' @export
nobs.btd <- function(object, ...) object$n_pairs

#' @export
confint.btd <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  se <- sqrt(diag(vcov(object)))
  if (missing(parm)) parm <- names(cf)
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  out <- cbind(cf[parm] - z * se[parm], cf[parm] + z * se[parm])
  dimnames(out) <- list(parm, paste0(format(100 * c(a, 1 - a), trim = TRUE,
                                            digits = 3), " %"))
  out
}

#' @rdname btd
#' @param exponentiate If `TRUE`, the `Estimate` column reports `exp(beta)`
#'   (the preference odds ratio for covariate terms, and `delta` for the
#'   `log_delta` row). `Std. Error`, `z value` and `Pr(>|z|)` remain on the
#'   log scale.
#' @param object A fitted `btd` object.
#' @param ... Unused.
#' @export
summary.btd <- function(object, exponentiate = FALSE, ...) {
  cf <- c(object$coefficients, if (object$est_delta) c(log_delta = object$log_delta))
  se <- sqrt(diag(object$vcov))[names(cf)]
  z  <- cf / se
  tab <- cbind(Estimate = cf, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  if (isTRUE(exponentiate)) tab[, "Estimate"] <- exp(tab[, "Estimate"])
  structure(list(call         = object$call,
                 coefficients = tab,
                 exponentiate = isTRUE(exponentiate),
                 delta        = object$delta,
                 loglik       = object$loglik,
                 n_subjects   = object$n_subjects,
                 n_pairs      = object$n_pairs,
                 n_ties       = object$n_ties,
                 cens_tie     = object$cens_tie,
                 converged    = object$converged,
                 n_iter       = object$n_iter,
                 engine       = object$engine,
                 storage      = object$storage),
            class = "summary.btd")
}

#' @export
print.summary.btd <- function(x, digits = 4, ...) {
  cat("Hierarchical Bradley-Terry-Davidson (BTD) model\n\nCall:\n",
      deparse(x$call), "\n\n", sep = "")
  cat(sprintf("%d subjects, %s pairs (%.1f%% tied); engine: %s, %s in %d iterations\n\n",
              x$n_subjects, format(x$n_pairs, big.mark = ","),
              100 * x$n_ties / x$n_pairs,
              paste0(x$engine, switch(x$storage,
                                      stream  = " (streaming, pairs not stored)",
                                      compact = " (compact pair codes)",
                                      stored  = " (pairs stored)",
                                      "")),
              if (x$converged) "converged" else "NOT converged", x$n_iter))

  coefs <- x$coefficients
if (isTRUE(x$exponentiate)) {
  colnames(coefs)[1L] <- "exp(Estimate)"
  rn <- rownames(coefs)
  rn[rn == "log_delta"] <- "delta"
  rownames(coefs) <- rn
}
  
  stats::printCoefmat(coefs, digits = digits,
                      signif.stars = getOption("show.signif.stars", TRUE),
                      has.Pvalue = TRUE, P.values = TRUE)

  cat("Standard errors: sandwich (U-statistic projection over subjects).\n")
  if (isTRUE(x$exponentiate))
  cat("Estimate column is on the exp scale; SE, z and p are on the log scale.\n")
  cat(sprintf("\nTie parameter delta = %.4f; log-likelihood = %.2f\n",
              x$delta, x$loglik))

  tie_total <- x$n_ties
  cens_tie  <- if (!is.null(x$cens_tie)) x$cens_tie else NA_real_
  if (tie_total > 0) {
    if (is.finite(cens_tie)) {
      pure_tie <- tie_total - cens_tie
      cat(sprintf("  Ties: %.1f%% of pairs, of which %.0f%% pure and %.0f%% censoring\n",
                  100 * tie_total / x$n_pairs,
                  100 * pure_tie  / tie_total,
                  100 * cens_tie  / tie_total))
    } else {
      cat(sprintf("  Ties: %.1f%% of pairs (censoring breakdown unavailable)\n",
                  100 * tie_total / x$n_pairs))
    }
  }
  invisible(x)
}

#' @export
print.btd <- function(x, digits = 4, ...) {
  cat("Hierarchical Bradley-Terry-Davidson (BTD) model\n\nCall:\n", deparse(x$call), "\n\n", sep = "")
  print(signif(x$coefficients, digits))
  cat(sprintf("\nTie parameter delta = %.4f (%d subjects, %s pairs)\n", x$delta,
              x$n_subjects, format(x$n_pairs, big.mark = ",")))
  invisible(x)
}

.design_new <- function(object, newdata) {
  mf <- stats::model.frame(object$terms, as.data.frame(newdata),
                           xlev = object$xlevels, na.action = stats::na.pass)
  X <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  icpt <- match("(Intercept)", colnames(X))
  if (!is.na(icpt)) X <- X[, -icpt, drop = FALSE]
  X
}

#' Predicted latent score and pair probabilities
#'
#' `predict()` returns the latent score \eqn{(x - \bar x)^\top \beta} (centred at
#' the training covariate means; only differences matter). `btd_pair_probs()`
#' returns the model probabilities that subject A wins, that the pair is tied,
#' and that subject B wins, for rows of `a` compared with the same rows of `b`.
#'
#' @param object A `btd` fit.
#' @param newdata Data frame of covariates.
#' @param a,b Data frames of equal length: first and second subject of each pair.
#' @param ... Unused.
#' @return A numeric vector, or a matrix with columns `win_a`, `tie`, `win_b`.
#' @examples
#' dat <- simulate_hps(N = 100, seed = 3)
#' fit <- btd(~ donor_type + age, dat, pairs = score_pairs(hps_hierarchy(), dat))
#' btd_pair_probs(fit, dat[1:3, ], dat[4:6, ])
#' @export
predict.btd <- function(object, newdata, ...) {
  if (missing(newdata)) return(object$eta)
  X <- .design_new(object, newdata)
  drop(sweep(X, 2L, object$x_center, "-") %*% object$coefficients)
}

#' @rdname predict.btd
#' @export
btd_pair_probs <- function(object, a, b) {
  if (nrow(a) != nrow(b)) stop("`a` and `b` must have the same number of rows.", call. = FALSE)
  d <- (predict(object, a) - predict(object, b)) / 2
  q <- exp(d) + exp(-d) + object$delta
  cbind(win_a = exp(d) / q, tie = object$delta / q, win_b = exp(-d) / q)
}

#' Calibration of the pairwise model
#'
#' Compares observed outcomes of the scored pairs with the model's predicted
#' probabilities.
#'
#' * `by = "level"`: among pairs *decided* at each hierarchy level, the observed
#'   share won by the first-listed subject versus the model's
#'   \eqn{P(i)/(P(i)+P(j))}. Systematic differences between levels suggest that
#'   a single latent scale does not describe all levels.
#' * `by = "bin"`: pairs are grouped by quantiles of the predicted latent
#'   difference; observed and predicted win / tie / loss fractions are compared.
#'   For fits made with streaming (no stored pairs) the bins are approximate
#'   quantile bins built from a fine histogram of the latent difference.
#'
#' The returned object is plotted by \code{\link{plot.btd_calibration}}, which
#' draws a dot plot for `by = "level"` and a three-panel observed-vs-predicted
#' line plot for `by = "bin"`.
#'
#' @param object A `btd` fit.
#' @param by `"level"` or `"bin"`.
#' @param n_bins Number of bins when `by = "bin"`.
#'
#' @return A `btd_calibration` object (a data frame with additional
#'   attributes) with columns
#'   * `by = "level"`: `level`, `n`, `obs_win_i`, `pred_win_i`.
#'   * `by = "bin"`: `bin`, `n`, `obs_win_i`, `pred_win_i`, `obs_tie`,
#'     `pred_tie`, `obs_win_j`, `pred_win_j`.
#'
#' @seealso \code{\link{plot.btd_calibration}} for the graphical counterpart;
#'   [btd()] to fit the model; [score_pairs()] to score the pairs.
#'
#' @examples
#' dat <- simulate_hps(N = 120, seed = 5)
#' h   <- hps_hierarchy()
#' fit <- btd(~ donor_type + dri, dat, hierarchy = h)
#'
#' ## Calibration by hierarchy level: for each level, the observed and
#' ## predicted share of pairs won by the first subject.
#' cal_level <- btd_calibration(fit, by = "level")
#' cal_level
#' plot(cal_level)
#'
#' ## Calibration by quantile bins of the latent difference d:
#' ## observed vs predicted win / tie / loss fractions.
#' cal_bin <- btd_calibration(fit, by = "bin", n_bins = 8)
#' cal_bin
#' plot(cal_bin)
#'
#' @export
btd_calibration <- function(object, by = c("level", "bin"), n_bins = 10L) {
  by <- match.arg(by)
  if (is.null(object$pairs)) return(.calibration_stream(object, by, n_bins))
  pw <- object$pairs
  d <- (object$eta[pw$id_i] - object$eta[pw$id_j]) / 2
  q <- exp(d) + exp(-d) + object$delta
  p_w <- exp(d) / q; p_l <- exp(-d) / q; p_t <- object$delta / q
  K <- length(pw$hierarchy)

  if (by == "level") {
    rows <- lapply(seq_len(K), function(k) {
      s <- pw$level == k
      data.frame(level = pw$hierarchy$level_names[k], n = sum(s),
                 obs_win_i = if (any(s)) mean(pw$W[s] == 1L) else NA_real_,
                 pred_win_i = if (any(s)) mean((p_w / (p_w + p_l))[s]) else NA_real_)
    })
    tab <- do.call(rbind, rows)
  } else {
    br <- unique(stats::quantile(d, seq(0, 1, length.out = n_bins + 1L)))
    bin <- cut(d, br, include.lowest = TRUE)
    f <- function(v) as.numeric(tapply(v, bin, mean))
    tab <- data.frame(bin = levels(bin), n = as.numeric(table(bin)),
                      obs_win_i = f(pw$W == 1L), pred_win_i = f(p_w),
                      obs_tie = f(pw$W == 0L), pred_tie = f(p_t),
                      obs_win_j = f(pw$W == -1L), pred_win_j = f(p_l))
  }

  structure(tab,
            class      = c("btd_calibration", "data.frame"),
            by         = by,
            n_bins     = n_bins,
            delta      = object$delta,
            n_subjects = object$n_subjects,
            n_pairs    = object$n_pairs)
}

## Calibration without stored pairs: sums are accumulated on the fly. Bins are
## formed by merging 512 fine equal-width bins of the latent difference so that
## each bin holds about the same number of pairs (approximate quantile bins).
.calibration_stream <- function(object, by, n_bins) {
  st <- object$stream
  h  <- st$hierarchy
  cal <- cpp_stream_calib(st$specs, st$comp$V, object$eta, object$delta, 512L,
                          .NCHUNK, st$threads)

  if (by == "level") {
    K <- length(h)
    grp <- factor(h$leaf_group, levels = seq_len(K))
    sm  <- function(col) as.numeric(tapply(cal$by_level[, col], grp, sum, default = 0))
    n <- sm(1L)
    tab <- data.frame(level = h$level_names, n = n,
                      obs_win_i  = ifelse(n > 0, sm(2L) / n, NA_real_),
                      pred_win_i = ifelse(n > 0, sm(3L) / n, NA_real_))
  } else {
    fb  <- cal$fine
    tot <- sum(fb[, 1L])
    cb  <- cumsum(fb[, 1L]) - fb[, 1L]
    grp <- pmin(n_bins, floor((cb + fb[, 1L] / 2) / tot * n_bins) + 1L)
    keep <- fb[, 1L] > 0
    grp  <- grp[keep]; fb <- fb[keep, , drop = FALSE]
    edges <- -cal$dmax + cal$width * (which(keep) - 1L)
    agg   <- rowsum(fb, grp, reorder = TRUE)
    lo    <- tapply(edges, grp, min); hi <- tapply(edges + cal$width, grp, max)
    n     <- agg[, 1L]
    tab <- data.frame(bin = sprintf("[%.3g, %.3g]", lo, hi), n = n,
                      obs_win_i = agg[, 2L] / n, pred_win_i = agg[, 5L] / n,
                      obs_tie   = agg[, 3L] / n, pred_tie   = agg[, 6L] / n,
                      obs_win_j = agg[, 4L] / n, pred_win_j = agg[, 7L] / n,
                      row.names = NULL)
  }

  structure(tab,
            class      = c("btd_calibration", "data.frame"),
            by         = by,
            n_bins     = n_bins,
            delta      = object$delta,
            n_subjects = object$n_subjects,
            n_pairs    = object$n_pairs)
}

# user arguments in `dots` override our defaults
.merge_plot_args <- function(defaults, dots) {
  nms <- names(dots)
  if (is.null(nms)) return(defaults)
  keep <- nzchar(nms)
  for (nm in nms[keep]) defaults[[nm]] <- dots[[nm]]
  defaults
}

.plot_calib_level <- function(x, ...) {
  dots <- list(...)
  k  <- nrow(x)
  yi <- rev(seq_len(k))
  lo <- min(c(x$obs_win_i, x$pred_win_i), na.rm = TRUE)
  hi <- max(c(x$obs_win_i, x$pred_win_i), na.rm = TRUE)
  rng <- c(max(0, lo - 0.05), min(1, hi + 0.05))

  op <- graphics::par(mar = c(4, 12, 3, 2))
  on.exit(par(op), add = TRUE)

  do.call(graphics::plot, .merge_plot_args(
    list(x = NA, xlim = rng, ylim = c(0.5, k + 0.5),
         yaxt = "n", xlab = "P(i wins | resolved at level)", ylab = "",
         main = "Calibration by level"),
    dots))

  graphics::axis(2, at = yi, labels = x$level, las = 1)
  graphics::abline(h = yi, col = "grey85", lty = 3)
  graphics::points(x$pred_win_i, yi, pch = 1, col = "red", cex = 1.2, lwd = 2)
  graphics::points(x$obs_win_i,  yi, pch = 16, col = "black")
  graphics::legend("bottomright",
                   legend = c("observed", "predicted"),
                   pch = c(16, 1), col = c("black", "red"), bty = "n")
}

.plot_calib_bin <- function(x, ...) {
  dots <- list(...)
  k  <- nrow(x)
  xi <- seq_len(k)

  # wider bottom margin for las = 2 bin labels
  op <- graphics::par(mfrow = c(1, 3), mar = c(7, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op), add = TRUE)

  ymax <- max(c(x$obs_win_i, x$obs_tie, x$obs_win_j,
                x$pred_win_i, x$pred_tie, x$pred_win_j), na.rm = TRUE)
  ylim <- c(0, 1.05 * ymax)

  panels <- list(
    list(obs = x$obs_win_i, pred = x$pred_win_i, main = "P(i wins)"),
    list(obs = x$obs_tie,   pred = x$pred_tie,   main = "P(tie)"),
    list(obs = x$obs_win_j, pred = x$pred_win_j, main = "P(j wins)")
  )

  for (p in panels) {
    do.call(graphics::plot, .merge_plot_args(
      list(x = xi, y = p$obs, type = "n", ylim = ylim, xaxt = "n",
           xlab = "", ylab = "Probability", main = p$main),
      dots))
    graphics::axis(1, at = xi, labels = x$bin, las = 2, cex.axis = 0.8)
    graphics::lines(xi, p$obs,  type = "b", pch = 16, col = "black")
    graphics::lines(xi, p$pred, type = "b", pch = 1,  col = "red", lty = 2)
  }
  graphics::mtext("Calibration by bin: observed (filled) vs predicted (open)",
                  side = 3, outer = TRUE, line = 0.5)
}

#' Plot a calibration table
#'
#' Two graphical summaries are produced depending on the `by` argument
#' of [btd_calibration()]:
#'
#' * `by = "level"`: a dot plot with one row per hierarchy level, showing
#'   the observed and predicted win share for pairs resolved at that level.
#'   Lines of agreement between the two are what a well-calibrated hierarchy
#'   should look like.
#' * `by = "bin"`: three line plots (win, tie, loss) of observed against
#'   predicted probabilities across quantile bins of the latent difference
#'   `d`. 
#'
#' @param x A `btd_calibration` object from [btd_calibration()].
#' @param ... Additional graphical parameters passed to the main `plot()`
#'   call of each panel. User-supplied `main`, `xlab`, `ylab` and `ylim`
#'   override the defaults.
#'
#' @details The bin plot shows observed and predicted probabilities only;
#'   no binomial error bars are drawn. Pairs share subjects, so treating
#'   them as independent Bernoulli draws would understate the sampling
#'   variability. For a variability estimate use the subject-clustered
#'   standard errors from [btd()] rather than a pair-level binomial
#'   interval.
#'   
#' @return `x`, invisibly.
#' @seealso [btd_calibration()]
#' @examples
#' dat <- simulate_hps(N = 120, seed = 5)
#' fit <- btd(~ donor_type + dri, dat, hierarchy = hps_hierarchy())
#'
#' plot(btd_calibration(fit, by = "level"))
#' plot(btd_calibration(fit, by = "bin", n_bins = 8))
#' @export
plot.btd_calibration <- function(x, ...) {
  by <- attr(x, "by")
  if (is.null(by)) by <- "bin"
  if (by == "level") .plot_calib_level(x, ...) else .plot_calib_bin(x, ...)
  invisible(x)
}

## ------------------------------------------------------------------
## Pure R reference implementation (validation; O(n_pairs * p) memory)
## ------------------------------------------------------------------

.btd_terms_r <- function(Xd, w, beta, lg, est) {
  d <- drop(Xd %*% beta) / 2
  u <- exp(d); v <- exp(-d)
  delta <- if (est) exp(lg) else 0
  q <- u + v + delta
  pw <- u / q; pl <- v / q; pt <- delta / q
  tie <- as.numeric(w == 0L)
  ll <- sum(d * w + (if (est) tie * lg else 0) - log(q))
  a <- w - (pw - pl)
  S <- 0.5 * Xd * a
  g <- colSums(S)
  h <- 4 * pw * pl + pt * (pw + pl)
  Hn <- 0.25 * crossprod(Xd, Xd * h)
  if (est) {
    st <- tie - pt
    S <- cbind(S, st)
    g <- c(g, sum(st))
    hbg <- -0.5 * crossprod(Xd, pt * (pw - pl))
    Hn <- rbind(cbind(Hn, hbg), c(hbg, sum(pt * (pw + pl))))
  }
  list(ll = ll, g = g, Hn = Hn, S = S)
}

.btd_fit_r <- function(X, id_i, id_j, w, est, control) {
  p <- ncol(X)
  Xd <- X[id_i, , drop = FALSE] - X[id_j, , drop = FALSE]
  w <- as.numeric(w)
  th <- rep(0, p + est)
  lambda <- control$lambda0
  conv <- FALSE
  for (iter in seq_len(control$max_iter)) {
    ev <- .btd_terms_r(Xd, w, th[seq_len(p)], if (est) th[p + 1L] else 0, est)
    gn <- max(abs(ev$g)) / (1 + abs(ev$ll))
    if (gn < control$tol) { conv <- TRUE; break }
    ok <- FALSE
    for (tr in 1:25) {
      step <- tryCatch(solve(ev$Hn + lambda * diag(length(th)), ev$g),
                       error = function(e) NULL)
      if (is.null(step)) { lambda <- lambda * 10; next }
      thn <- th + step
      eta_n <- drop(X %*% thn[seq_len(p)])
      if (!all(is.finite(thn)) || (max(eta_n) - min(eta_n)) / 2 > 100) {
        lambda <- lambda * 10; next
      }
      lln <- .btd_terms_r(Xd, w, thn[seq_len(p)], if (est) thn[p + 1L] else 0, est)$ll
      if (is.finite(lln) && lln > ev$ll) {
        th <- thn; lambda <- max(lambda / 3, 1e-9); ok <- TRUE; break
      }
      lambda <- lambda * 10
    }
    if (!ok) break
  }
  ev <- .btd_terms_r(Xd, w, th[seq_len(p)], if (est) th[p + 1L] else 0, est)
  gn <- max(abs(ev$g)) / (1 + abs(ev$ll))
  list(beta = th[seq_len(p)], log_delta = if (est) th[p + 1L] else NA_real_,
       loglik = ev$ll, converged = conv || gn < control$tol, n_iter = iter,
       grad_norm = gn, lambda_end = lambda)
}

.btd_vcov_r <- function(X, id_i, id_j, w, beta, lg, est) {
  N <- nrow(X)
  Xd <- X[id_i, , drop = FALSE] - X[id_j, , drop = FALSE]
  ev <- .btd_terms_r(Xd, as.numeric(w), beta, lg, est)
  n <- length(w)
  IF <- rowsum(rbind(ev$S, ev$S), group = c(id_i, id_j), reorder = TRUE)
  IF <- (N / n) * IF
  A <- ev$Hn / n
  B <- crossprod(IF) / N
  singular <- FALSE
  Ainv <- tryCatch(solve(A), error = function(e) {
    singular <<- TRUE
    solve(A + diag(1e-10, nrow(A)))
  })
  Vm <- tryCatch(solve(ev$Hn), error = function(e) {
    singular <<- TRUE
    solve(ev$Hn + diag(1e-10, nrow(ev$Hn)))
  })
  list(V_sandwich = Ainv %*% B %*% t(Ainv) / N, V_model = Vm, singular = singular)
}
