## Tests for the 1.0.4 fixes.

h4 <- hps_hierarchy()

test_that("absurd thread counts are capped instead of crashing", {
  d <- simulate_hps(N = 60, seed = 1)
  ref <- btd(~ age, d, hierarchy = h4, storage = "compact")
  for (th in c(0, 1e5, 1e9)) {
    f <- btd(~ age, d, hierarchy = h4, storage = "compact", control = btd_control(threads = th))
    expect_identical(f$coefficients, ref$coefficients)
  }
  op <- options(hpsreg.threads = 1e7); on.exit(options(op))
  expect_equal(resolution(h4, data = d)$n, resolution(score_pairs(h4, d))$n)
  expect_error(btd(~ age, d, hierarchy = h4, storage = "compact",
                   control = btd_control(threads = NA)), "single number")
  expect_lte(hpsreg:::.n_threads(list(threads = 500)), hpsreg:::.NCHUNK)
})

test_that("by-level calibration is oriented by the fitted score", {
  d <- simulate_hps(N = 300, seed = 11)
  pw <- score_pairs(h4, d)
  tight <- btd_control(tol = 1e-12)
  fit <- btd(~ donor_type + age + dri + hct_ci, d, pairs = pw, control = tight)
  cl <- btd_calibration(fit)
  ## brute force for every level
  dd <- (fit$eta[pw$id_i] - fit$eta[pw$id_j]) / 2
  q <- exp(dd) + exp(-dd) + fit$delta
  eq <- abs(dd) <= 1e-12
  expect_true(any(eq))                       # identical profiles exist in this cohort
  hi <- ifelse(dd > 0, 1L, -1L)
  ohi <- ifelse(eq, 0.5, as.numeric(pw$W * hi == 1L))
  phi <- ifelse(eq, 0.5, ifelse(dd > 0, exp(dd), exp(-dd)) / (exp(dd) + exp(-dd)))
  for (k in seq_len(nrow(cl))) {
    s <- pw$level == k
    expect_equal(cl$obs_win_higher[k], mean(ohi[s]))
    expect_equal(cl$pred_win_higher[k], mean(phi[s]))
  }
  ## does not depend on row order (the first-listed version did)
  set.seed(3); perm <- sample.int(nrow(d))
  cl2 <- btd_calibration(btd(~ donor_type + age + dri + hct_ci, d[perm, ], hierarchy = h4,
                             control = tight))
  expect_equal(cl$obs_win_higher, cl2$obs_win_higher, tolerance = 1e-8)
  expect_equal(cl$pred_win_higher, cl2$pred_win_higher, tolerance = 1e-6)
  ## stored, compact and stream agree
  for (st in c("compact", "stream")) {
    cs <- btd_calibration(btd(~ donor_type + age + dri + hct_ci, d, hierarchy = h4, storage = st,
                              control = tight))
    expect_equal(cs$obs_win_higher, cl$obs_win_higher, tolerance = 1e-10)
    expect_equal(cs$pred_win_higher, cl$pred_win_higher, tolerance = 1e-8)
  }
  ## a model with real covariates predicts winners better than a pure-noise model
  set.seed(5); d$noise <- rnorm(nrow(d))
  cn <- btd_calibration(btd(~ noise, d, pairs = pw))
  expect_gt(cl$obs_win_higher[1], cn$obs_win_higher[1])
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(cl))
})

test_that("pre-computed pairs are rejected for data that does not match", {
  d <- simulate_hps(N = 80, seed = 2)
  pw <- score_pairs(h4, d)
  expect_silent(btd(~ donor_type + age, d, pairs = pw))
  expect_error(btd(~ donor_type + age, d[80:1, ], pairs = pw), "does not match")
  d2 <- d; d2$death_day[which(!is.na(d2$death_day))[1]] <- 1
  expect_error(btd(~ donor_type + age, d2, pairs = pw), "does not match")
  ## covariate-only data (outcome columns absent) cannot be checked and is accepted
  expect_silent(btd(~ donor_type + age, d[, c("donor_type", "age")], pairs = pw))
  ## changing a covariate does not trip the check
  d3 <- d; d3$age <- d3$age + 1
  expect_silent(btd(~ donor_type + age, d3, pairs = pw))
})

test_that("separation is reported, with the latent range", {
  set.seed(5); n <- 150
  s <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  s$y <- s$x1 + 0.5 * s$x2 + rnorm(n, 0, 0.03)
  hs <- hierarchy(lvl_continuous("y", better = "higher", threshold = 0.1))
  expect_warning(f <- btd(~ x1 + x2, s, hierarchy = hs), "separation")
  expect_gt(f$latent_range, 40)
  d <- simulate_hps(N = 150, seed = 7)
  expect_no_warning(g <- btd(~ donor_type + age, d, hierarchy = h4))
  expect_lt(g$latent_range, 20)
})

test_that("an event and a censoring at the same time follow the usual convention", {
  d <- data.frame(t = c(10, 10, 99), s = c(1, 0, 0))
  p <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "later")), d)
  expect_equal(p$W[1], -1L)                       # censored at 10 outlived the event at 10
  expect_equal(p$uninf[1], 0L)
  pu <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "later",
                                      same_time = "uninformative")), d)
  expect_equal(pu$W[1], 0L); expect_equal(pu$uninf[1], 1L)
  pe <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "earlier")), d)
  expect_equal(pe$W[1], 1L)                       # earlier is better: the event wins
  ## with a positive margin the ordering is not certain beyond the margin
  pm <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "later", threshold = 2)), d)
  expect_equal(pm$W[1], 0L); expect_equal(pm$uninf[1], 1L)
  ## the reverse orientation (censored subject listed first)
  p2 <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "later")), d[c(2, 1, 3), ])
  expect_equal(p2$W[1], 1L)
  expect_error(lvl_tte("t", status = "s", better = "later", same_time = "x"))
})

test_that("an event exactly at the horizon is still an event", {
  hz <- data.frame(t = c(180, NA, 100, 250))
  p <- score_pairs(hierarchy(lvl_tte("t", better = "later", horizon = 180)), hz)
  w <- function(i, j) p$W[p$id_i == i & p$id_j == j]
  expect_equal(w(1, 2), -1L)          # death on day 180 loses to no death
  expect_equal(w(1, 4), -1L)          # death after the horizon counts as event-free
  expect_equal(w(2, 4), 0L)           # both event-free through the horizon
  expect_equal(w(1, 3), 1L)
  hs <- data.frame(t = c(180, 180, 200, 120, 150), s = c(1, 0, 1, 0, 1))
  ps <- score_pairs(hierarchy(lvl_tte("t", status = "s", better = "later", horizon = 180)), hs)
  ws <- function(i, j) ps$W[ps$id_i == i & ps$id_j == j]
  expect_equal(ws(1, 2), -1L)         # event at 180 vs censored at 180 (event-free through 180)
  expect_equal(ws(1, 3), -1L)         # event at 180 vs event after the horizon
  expect_equal(ws(2, 3), 0L)          # both event-free through the horizon
  expect_equal(ws(3, 4), 0L)          # event-free through 180 vs censored at 120: unknown
  expect_equal(ps$uninf[ps$id_i == 3 & ps$id_j == 4], 1L)
  expect_equal(ws(3, 5), 1L)
  ## no horizon: unchanged behaviour for the transplant example
  d <- simulate_hps(N = 120, seed = 99)
  expect_identical(score_pairs(h4, d)$W, build_pairwise_fast(d)$W)
})

test_that("silent level specifications now warn", {
  d <- simulate_hps(N = 30, seed = 1)
  expect_warning(score_pairs(hierarchy(lvl_binary("donor_type", favorable = "XYZ")), d),
                 "does not occur")
  expect_warning(score_pairs(hierarchy(lvl_ordinal("gvhd_grade", better = "lower",
                                                   order = c(2, 0, 1))), d), "ignored")
})

test_that("summary notes what the delta-row test is", {
  d <- simulate_hps(N = 100, seed = 3)
  f <- btd(~ age, d, hierarchy = h4)
  expect_output(print(summary(f)), "delta = 1")
})
