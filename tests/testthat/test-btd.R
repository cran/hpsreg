## ---- helpers ---------------------------------------------------------------
fd_grad <- function(fun, x, h = 1e-6) {
  vapply(seq_along(x), function(k) {
    e <- replace(numeric(length(x)), k, h)
    (fun(x + e) - fun(x - e)) / (2 * h)
  }, 0)
}

small_problem <- function(N = 40, p = 3, seed = 1, ties = TRUE) {
  set.seed(seed)
  X <- scale(matrix(rnorm(N * p), N, p))
  X <- unclass(X); attributes(X) <- list(dim = c(N, p))
  ids <- hpsreg:::cpp_all_pairs(N)
  beta <- c(0.6, -0.4, 0.2)[seq_len(p)]
  d <- (X[ids$i, ] %*% beta - X[ids$j, ] %*% beta) / 2
  delta <- if (ties) 0.5 else 0
  q <- exp(d) + exp(-d) + delta
  pr <- cbind(exp(d) / q, delta / q, exp(-d) / q)
  W <- apply(pr, 1, function(pp) sample(c(1L, 0L, -1L), 1, prob = pp))
  list(X = X, i = ids$i, j = ids$j, W = as.integer(W), beta = beta, delta = delta)
}

## ---- analytic derivatives -----------------------------------------------------
test_that("gradient and information match finite differences", {
  pr <- small_problem()
  th <- c(0.3, -0.2, 0.1, -0.5)
  f <- function(th) hpsreg:::cpp_btd_eval(pr$X, pr$i, pr$j, pr$W, th[1:3], th[4], TRUE, FALSE)$loglik
  ev <- hpsreg:::cpp_btd_eval(pr$X, pr$i, pr$j, pr$W, th[1:3], th[4], TRUE, TRUE)
  expect_equal(as.numeric(ev$grad), fd_grad(f, th), tolerance = 1e-5)
  gfun <- function(th) as.numeric(hpsreg:::cpp_btd_eval(pr$X, pr$i, pr$j, pr$W, th[1:3], th[4], TRUE, FALSE)$grad)
  H <- vapply(1:4, function(k) {
    e <- replace(numeric(4), k, 1e-6)
    (gfun(th + e) - gfun(th - e)) / 2e-6
  }, numeric(4))
  expect_equal(unname(ev$neg_hessian), unname(-H), tolerance = 1e-5)
})

test_that("Bradley-Terry mode (no ties) has correct derivatives", {
  pr <- small_problem(ties = FALSE, seed = 2)
  b <- c(0.2, 0.1, -0.3)
  f <- function(b) hpsreg:::cpp_btd_eval(pr$X, pr$i, pr$j, pr$W, b, 0, FALSE, FALSE)$loglik
  ev <- hpsreg:::cpp_btd_eval(pr$X, pr$i, pr$j, pr$W, b, 0, FALSE, TRUE)
  expect_equal(as.numeric(ev$grad), fd_grad(f, b), tolerance = 1e-5)
  expect_equal(dim(ev$neg_hessian), c(3L, 3L))
})

## ---- engines agree -----------------------------------------------------------
test_that("C++ and pure R engines agree on estimates and variances", {
  dat <- simulate_hps(N = 100, seed = 17)
  pw <- score_pairs(hps_hierarchy(), dat)
  fml <- ~ donor_type + age + dri + hct_ci + conditioning
  a <- btd(fml, dat, pairs = pw, engine = "cpp")
  b <- btd(fml, dat, pairs = pw, engine = "r")
  expect_true(a$converged && b$converged)
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_equal(as.numeric(a$log_delta), as.numeric(b$log_delta), tolerance = 1e-7)
  expect_equal(a$vcov, b$vcov, tolerance = 1e-6)
  expect_equal(a$vcov_model, b$vcov_model, tolerance = 1e-6)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-9)
})

test_that("estimates agree with gnm (conditional Poisson formulation)", {
  skip_if_not_installed("gnm")
  skip_on_cran()
  dat <- simulate_hps(N = 45, seed = 8)
  pw <- score_pairs(hps_hierarchy(), dat)
  fit <- btd(~ donor_type + age, dat, pairs = pw)

  X <- scale(stats::model.matrix(~ donor_type + age, dat)[, -1])
  sc <- attr(X, "scaled:scale"); p <- ncol(X); npair <- length(pw$W)
  rows <- 3L * npair
  iw <- 3L * seq_len(npair) - 2L; jw <- iw + 1L; it <- iw + 2L
  Xd <- matrix(0, rows, p + 1L)
  Xd[iw, 1:p] <- X[pw$id_i, ]; Xd[jw, 1:p] <- X[pw$id_j, ]
  Xd[it, 1:p] <- (X[pw$id_i, ] + X[pw$id_j, ]) / 2; Xd[it, p + 1L] <- 1
  colnames(Xd) <- c(paste0("b", 1:p), "ld")
  cnt <- integer(rows); cnt[iw] <- pw$W == 1L; cnt[jw] <- pw$W == -1L; cnt[it] <- pw$W == 0L
  df <- data.frame(count = cnt, pair = factor(rep(seq_len(npair), each = 3L)), Xd)
  g <- gnm::gnm(count ~ 0 + ld + b1 + b2 + b3, eliminate = pair, family = poisson,
                data = df, iterMax = 500, trace = FALSE)
  cf <- stats::coef(g)
  expect_equal(unname(cf[c("b1", "b2", "b3")] / sc), unname(coef(fit)), tolerance = 1e-4)
  expect_equal(unname(cf["ld"]), fit$log_delta, tolerance = 1e-4)
})

## ---- statistical behaviour ---------------------------------------------------
test_that("estimator recovers known parameters", {
  pr <- small_problem(N = 120, seed = 4)
  dat <- as.data.frame(pr$X); names(dat) <- c("x1", "x2", "x3")
  pw <- structure(list(id_i = pr$i, id_j = pr$j, W = pr$W, level = rep(1L, length(pr$W)),
                       leaf = NULL, uninf = integer(length(pr$W)), N = 120L, complete = TRUE,
                       hierarchy = hierarchy(lvl_binary("x1"))), class = "hb_pairs")
  fit <- btd(~ x1 + x2 + x3, dat, pairs = pw)
  se <- sqrt(diag(vcov(fit)))
  expect_lt(max(abs(coef(fit) - pr$beta) / se), 4)
  expect_lt(abs(fit$delta - 0.5), 0.1)
})

test_that("sandwich standard errors are close to Monte Carlo variability", {
  skip_on_cran()
  set.seed(77)
  N <- 60; reps <- 150
  est <- matrix(NA_real_, reps, 2); se <- est
  for (r in seq_len(reps)) {
    d <- data.frame(x = rnorm(N), z = rbinom(N, 1, 0.5))
    lat <- 0.7 * d$x - 0.5 * d$z + rnorm(N)         # hierarchy on a noisy latent score
    d$y <- lat
    pw <- score_pairs(hierarchy(lvl_continuous("y", better = "higher", threshold = 0.3)), d)
    f <- suppressMessages(btd(~ x + z, d, pairs = pw))
    est[r, ] <- coef(f); se[r, ] <- sqrt(diag(vcov(f)))
  }
  ratio <- colMeans(se) / apply(est, 2, sd)
  expect_true(all(ratio > 0.8 & ratio < 1.2))
})

test_that("summary.btd(exponentiate = TRUE) reports odds ratios", {
  dat <- simulate_hps(N = 100, seed = 21)
  fit <- btd(~ donor_type + age, dat, hierarchy = hps_hierarchy())

  s_log <- summary(fit)
  s_exp <- summary(fit, exponentiate = TRUE)

  expect_equal(s_exp$coefficients[, "Estimate"],
               exp(s_log$coefficients[, "Estimate"]),
               tolerance = 1e-12)
  expect_identical(s_exp$coefficients[, "Std. Error"],
                   s_log$coefficients[, "Std. Error"])
  expect_identical(s_exp$coefficients[, "z value"],
                   s_log$coefficients[, "z value"])
  expect_identical(s_exp$coefficients[, "Pr(>|z|)"],
                   s_log$coefficients[, "Pr(>|z|)"])
  expect_true(s_exp$exponentiate)
  expect_false(s_log$exponentiate)

  out <- capture.output(print(s_exp))
  expect_true(any(grepl("exp(Estimate)", out, fixed = TRUE)))
  expect_true(any(grepl("^delta\\s", out)))

  out_log <- capture.output(print(s_log))
  expect_false(any(grepl("exp(Estimate)", out_log, fixed = TRUE)))
})

## ---- interface ---------------------------------------------------------------
test_that("no tied pairs triggers Bradley-Terry fit with a message", {
  set.seed(9)
  d <- data.frame(x = rnorm(50)); d$y <- d$x + rnorm(50)
  pw <- score_pairs(hierarchy(lvl_continuous("y", better = "higher")), d)
  expect_message(f <- btd(~ x, d, pairs = pw), "Bradley-Terry")
  expect_false(f$est_delta)
  expect_equal(f$delta, 0)
  expect_equal(rownames(f$vcov), "x")
  expect_equal(nrow(summary(f)$coefficients), 1L)
})

test_that("btd interface: formulas, checks, methods", {
  dat <- simulate_hps(N = 90, seed = 31)
  h <- hps_hierarchy()
  expect_error(btd(~ age, dat), "pairs")
  expect_error(btd(age ~ 1, dat, pairs = score_pairs(h, dat)), "covariates")
  expect_error(btd(~ age, dat[1:50, ], pairs = score_pairs(h, dat)), "subjects")
  expect_error(btd(~ age, dat, pairs = list()), "score_pairs")
  dat$const <- 1
  expect_error(btd(~ age + const, dat, pairs = score_pairs(h, dat)), "no variation")
  # scoring inline via `hierarchy`
  f1 <- btd(~ donor_type + age, dat, hierarchy = h)
  f2 <- btd(~ donor_type + age, dat, pairs = score_pairs(h, dat))
  expect_equal(coef(f1), coef(f2))

  expect_equal(names(coef(f2)), c("donor_typeMUD", "donor_typeHID", "age"))
  expect_equal(dim(vcov(f2)), c(3L, 3L))
  expect_equal(dim(f2$vcov), c(4L, 4L))
  ci <- confint(f2)
  expect_equal(dim(ci), c(3L, 2L))
  expect_true(all(ci[, 1] < coef(f2) & coef(f2) < ci[, 2]))
  expect_s3_class(logLik(f2), "logLik")
  expect_equal(nobs(f2), 4005L)          # pairs: 90 * 89 / 2
  expect_equal(nobs(f2), f2$n_pairs)
  expect_equal(nobs(f2), attr(logLik(f2), "nobs"))
  expect_equal(f2$n_subjects, 90L)
  expect_output(print(f2), "Davidson")
  expect_output(print(summary(f2)), "sandwich")
  # scale = FALSE gives the same estimates
  f3 <- btd(~ donor_type + age, dat, pairs = score_pairs(h, dat), scale = FALSE)
  expect_equal(coef(f2), coef(f3), tolerance = 1e-6)
  # covariate rescaling transforms coefficients and SEs
  dat$age10 <- dat$age / 10
  f4 <- btd(~ donor_type + age10, dat, pairs = score_pairs(h, dat))
  expect_equal(unname(coef(f4)["age10"]), unname(10 * coef(f2)["age"]), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f4)))["age10"]),
               unname(10 * sqrt(diag(vcov(f2)))["age"]), tolerance = 1e-6)
})

test_that("predict and pair probabilities are consistent", {
  dat <- simulate_hps(N = 80, seed = 12)
  fit <- btd(~ donor_type + age + dri, dat, pairs = score_pairs(hps_hierarchy(), dat))
  s <- predict(fit, dat)
  expect_equal(s, fit$eta, tolerance = 1e-10)
  pr <- btd_pair_probs(fit, dat[1:10, ], dat[11:20, ])
  expect_equal(unname(rowSums(pr)), rep(1, 10), tolerance = 1e-12)
  d <- (s[1:10] - s[11:20]) / 2
  expect_equal(unname(pr[, "win_a"]), unname(exp(d) / (exp(d) + exp(-d) + fit$delta)))
  expect_error(btd_pair_probs(fit, dat[1:3, ], dat[1:4, ]), "same number")
  expect_equal(predict(fit, dat[1:5, ]), s[1:5], tolerance = 1e-10)
})

test_that("calibration summaries have the right shape and are sensible", {
  dat <- simulate_hps(N = 120, seed = 6)
  fit <- btd(~ donor_type + dri + age, dat, pairs = score_pairs(hps_hierarchy(), dat))
  cl <- btd_calibration(fit)
  expect_equal(nrow(cl), 4L)
  expect_true(all(cl$obs_win_i > 0.3 & cl$obs_win_i < 0.7))
  cb <- btd_calibration(fit, by = "bin", n_bins = 5)
  expect_equal(nrow(cb), 5L)
  expect_equal(sum(cb$n), length(fit$pairs$W))
  expect_equal(cb$obs_win_i + cb$obs_tie + cb$obs_win_j, rep(1, 5), tolerance = 1e-10)
  # model-implied overall tie rate equals observed (score equation for delta)
  expect_equal(sum(cb$n * cb$obs_tie) / sum(cb$n), sum(cb$n * cb$pred_tie) / sum(cb$n),
               tolerance = 1e-6)
})

## ---- calibration class and plotting ------------------------------------------
test_that("btd_calibration carries its class for every storage engine", {
  dat <- simulate_hps(N = 120, seed = 5)
  h   <- hps_hierarchy()

  for (storage in c("stored", "compact", "stream")) {
    fit <- btd(~ donor_type + dri, dat, hierarchy = h, storage = storage)
    for (by in c("level", "bin")) {
      cal <- btd_calibration(fit, by = by, n_bins = 5L)
      expect_s3_class(cal, "btd_calibration")
      expect_s3_class(cal, "data.frame")
      expect_equal(attr(cal, "by"), by)
    }
  }
})

test_that("plot.btd_calibration accepts user args and dispatches once", {
  m <- getS3method("plot", "btd_calibration")
  expect_true(is.function(m))

  dat <- simulate_hps(N = 100, seed = 6)
  fit <- btd(~ donor_type + dri, dat, hierarchy = hps_hierarchy())

  for (by in c("level", "bin")) {
    cal <- btd_calibration(fit, by = by, n_bins = 4L)
    f <- tempfile(fileext = ".pdf")
    grDevices::pdf(f); on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(plot(cal))                              # defaults work
    expect_silent(plot(cal, main = "custom"))             # user main wins
    expect_silent(plot(cal, ylim = c(0, 1)))              # user ylim wins
    expect_silent(plot(cal, main = "x", ylim = c(0, 1)))  # both
  }
})

test_that("plot.btd_calibration works on compact and stream fits too", {
  dat <- simulate_hps(N = 90, seed = 7)
  h   <- hps_hierarchy()
  for (storage in c("compact", "stream")) {
    fit <- btd(~ donor_type + dri, dat, hierarchy = h, storage = storage)
    cal <- btd_calibration(fit, by = "level")
    expect_s3_class(cal, "btd_calibration")
    f <- tempfile(fileext = ".pdf")
    grDevices::pdf(f); on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(plot(cal))
  }
})

## ---- censoring ties ----------------------------------------------------------
test_that("cens_tie agrees between stored, compact and stream paths", {
  dat <- simulate_hps(N = 150, seed = 1)
  h   <- hps_hierarchy()

  fit_stored  <- btd(~ donor_type + age, data = dat,
                     pairs = score_pairs(h, dat))
  fit_compact <- btd(~ donor_type + age, data = dat, hierarchy = h,
                     storage = "compact")
  fit_stream  <- btd(~ donor_type + age, data = dat, hierarchy = h,
                     storage = "stream")

  expect_equal(fit_stored$cens_tie, fit_compact$cens_tie)
  expect_equal(fit_stored$cens_tie, fit_stream$cens_tie)
  expect_true(fit_stored$cens_tie >= 0)
  expect_true(fit_stored$cens_tie <= fit_stored$n_ties)
})

test_that("stored-path cens_tie equals the C++ accumulator", {
  dat <- simulate_hps(N = 120, seed = 2)
  h   <- hps_hierarchy()
  fit <- btd(~ donor_type + age, data = dat, pairs = score_pairs(h, dat))

  comp <- hpsreg:::.compile_hierarchy(h, dat)
  cn   <- cpp_stream_counts(hpsreg:::.stream_specs(comp), comp$V,
                            hpsreg:::.NCHUNK, 1L)
  expect_equal(fit$cens_tie, cn$cens_tie)
})

test_that("print.summary surfaces a missing censoring breakdown", {
  x <- structure(list(
    call = quote(fit),
    coefficients = matrix(0, 1, 4,
                          dimnames = list("a", c("Estimate", "Std. Error",
                                                 "z value", "Pr(>|z|)"))),
    delta = 0.5, loglik = -1,
    n_subjects = 10L, n_pairs = 45L, n_ties = 5L,
    converged = TRUE, n_iter = 3L,
    engine = "cpp", storage = "stored"),
    class = "summary.btd")
  out <- capture.output(print(x))
  expect_true(any(grepl("censoring breakdown unavailable", out)))
})

## ---- tile and thread invariance ----------------------------------------------
test_that("tile size does not change the compact sandwich variance", {
  dat <- simulate_hps(N = 200, seed = 4)
  h   <- hps_hierarchy()
  fit <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact")

  mf <- stats::model.frame(fit$terms, dat)
  X0 <- stats::model.matrix(fit$terms, mf)
  X  <- X0[, colnames(X0) != "(Intercept)", drop = FALSE]
  Xs <- sweep(sweep(X, 2L, fit$x_center, "-"), 2L, fit$x_scale, "/")

  comp  <- hpsreg:::.compile_hierarchy(h, dat)
  codes <- hpsreg:::cpp_compact_build(hpsreg:::.stream_specs(comp), comp$V,
                                      hpsreg:::.NCHUNK, 1L)$codes
  beta_std <- fit$coefficients * fit$x_scale
  lg       <- if (fit$est_delta) fit$log_delta else 0

  v16  <- hpsreg:::cpp_compact_vcov(Xs, codes, beta_std, lg, fit$est_delta,
                                    hpsreg:::.NCHUNK, 1L, tile = 16L)
  v512 <- hpsreg:::cpp_compact_vcov(Xs, codes, beta_std, lg, fit$est_delta,
                                    hpsreg:::.NCHUNK, 1L, tile = 512L)
  expect_equal(v16$V_sandwich, v512$V_sandwich, tolerance = 1e-12)
})

test_that("thread count does not change the compact sandwich variance", {
  skip_if(openmp_threads() < 2L, "OpenMP not available")

  dat <- simulate_hps(N = 180, seed = 5)
  h   <- hps_hierarchy()
  fit <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact")

  mf <- stats::model.frame(fit$terms, dat)
  X0 <- stats::model.matrix(fit$terms, mf)
  X  <- X0[, colnames(X0) != "(Intercept)", drop = FALSE]
  Xs <- sweep(sweep(X, 2L, fit$x_center, "-"), 2L, fit$x_scale, "/")

  comp  <- hpsreg:::.compile_hierarchy(h, dat)
  codes <- hpsreg:::cpp_compact_build(hpsreg:::.stream_specs(comp), comp$V,
                                      hpsreg:::.NCHUNK, 1L)$codes
  beta_std <- fit$coefficients * fit$x_scale
  lg       <- if (fit$est_delta) fit$log_delta else 0

  v1 <- hpsreg:::cpp_compact_vcov(Xs, codes, beta_std, lg, fit$est_delta,
                                  hpsreg:::.NCHUNK, 1L)
  v4 <- hpsreg:::cpp_compact_vcov(Xs, codes, beta_std, lg, fit$est_delta,
                                  hpsreg:::.NCHUNK, 4L)
  expect_equal(v1$V_sandwich, v4$V_sandwich, tolerance = 0)
})

test_that("thread count does not change the fit", {
  skip_if(openmp_threads() < 2L, "OpenMP not available")

  dat <- simulate_hps(N = 150, seed = 8)
  h   <- hps_hierarchy()
  fit1 <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact",
              control = btd_control(threads = 1))
  fit4 <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact",
              control = btd_control(threads = 4))
  expect_identical(coef(fit1), coef(fit4))
  expect_identical(fit1$log_delta, fit4$log_delta)
  expect_identical(fit1$loglik, fit4$loglik)
})

## ---- fit diagnostics ---------------------------------------------------------
test_that("fit reports evaluation counts consistently", {
  dat <- simulate_hps(N = 100, seed = 5)
  fit <- btd(~ donor_type + age, data = dat, hierarchy = hps_hierarchy())

  expect_true(is.numeric(fit$n_full_evals) && fit$n_full_evals >= 1)
  expect_true(is.numeric(fit$n_grad_evals) && fit$n_grad_evals >= 0)
  expect_true(is.numeric(fit$trials_per_iter))
  expect_true(all(fit$trials_per_iter >= 1))
  expect_true(all(fit$trials_per_iter <= 10))       # 10 LM trials max per iter
  expect_true(length(fit$trials_per_iter) <= fit$n_iter)
  expect_true(fit$n_full_evals >= fit$n_iter)
  expect_true(fit$n_full_evals == fit$n_iter + 1L)  # one at start, one at the end
})

test_that("trials_per_iter sums to the gradient-only evaluations", {
  dat <- simulate_hps(N = 100, seed = 6)
  fit <- btd(~ donor_type + age + dri, data = dat, hierarchy = hps_hierarchy())
  expect_equal(sum(fit$trials_per_iter), fit$n_grad_evals)
})