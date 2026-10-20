dat_s <- simulate_hps(N = 130, seed = 41)
h_s <- hps_hierarchy()
fml_s <- ~ donor_type + age + dri + hct_ci

test_that("streaming fit equals the stored-pair fit", {
  a <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stored")
  b <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream")
  expect_equal(a$storage, "stored"); expect_equal(b$storage, "stream")
  expect_null(b$pairs)
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_equal(a$log_delta, b$log_delta, tolerance = 1e-7)
  expect_equal(a$vcov, b$vcov, tolerance = 1e-5)
  expect_equal(a$vcov_model, b$vcov_model, tolerance = 1e-5)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-9)
  expect_equal(a$n_pairs, b$n_pairs)
  expect_equal(a$n_ties, b$n_ties)
  expect_output(print(summary(b)), "streaming")
})

test_that("results are bitwise identical for any number of threads", {
  f1 <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream",
            control = btd_control(threads = 1))
  f3 <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream",
            control = btd_control(threads = 3))
  expect_identical(f1$coefficients, f3$coefficients)
  expect_identical(f1$vcov, f3$vcov)
  expect_identical(f1$loglik, f3$loglik)
  op <- options(hpsreg.threads = 2); on.exit(options(op))
  f2 <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream")
  expect_identical(f1$coefficients, f2$coefficients)
  expect_true(openmp_threads() >= 0)
})

test_that("storage = 'auto' switches on the number of pairs", {
  op <- options(hpsreg.stream_min_pairs = 2.5e6,
                hpsreg.compact_max_pairs = 1e9)
  on.exit(options(op), add = TRUE)

  expect_equal(btd(~ age, dat_s, hierarchy = h_s)$storage, "stored")

  options(hpsreg.stream_min_pairs = 100)
  expect_equal(btd(~ age, dat_s, hierarchy = h_s)$storage, "compact")

  options(hpsreg.compact_max_pairs = 100)
  expect_equal(btd(~ age, dat_s, hierarchy = h_s)$storage, "stream")

  # pre-computed pairs are always used as given
  pw <- score_pairs(h_s, dat_s)
  expect_equal(btd(~ age, dat_s, pairs = pw)$storage, "stored")
})


test_that("streaming argument checks", {
  f <- function(data, i, j) sign(data$age[i] - data$age[j])
  hc <- hierarchy(lvl_custom(f), lvl_binary("hct_ci", favorable = "0-2"))
  expect_error(btd(~ donor_type, dat_s, hierarchy = hc, storage = "stream"), "lvl_custom")
  # auto falls back to stored for custom levels even above the threshold
  op <- options(hpsreg.stream_min_pairs = 100); on.exit(options(op))
  expect_equal(btd(~ donor_type, dat_s, hierarchy = hc)$storage, "stored")
  expect_error(btd(~ age, dat_s, hierarchy = h_s, storage = "stream", engine = "r"),
               "engine")
  expect_error(btd(~ age, dat_s, pairs = score_pairs(h_s, dat_s), storage = "stream"),
               "hierarchy")
  expect_error(resolution(hc, data = dat_s), "lvl_custom")
  expect_error(resolution(h_s), "data")
})

test_that("streaming resolution and counts equal the stored ones", {
  pw <- score_pairs(h_s, dat_s)
  expect_equal(resolution(h_s, data = dat_s)$n, as.numeric(resolution(pw)$n))
  expect_equal(resolution(h_s, by = "leaf", data = dat_s)$n,
               as.numeric(resolution(pw, by = "leaf")$n))
  expect_equal(resolution(h_s, data = dat_s)$prop, resolution(pw)$prop)
  # tiny N (fewer rows than chunks)
  d5 <- dat_s[1:5, ]
  expect_equal(resolution(h_s, data = d5)$n, as.numeric(resolution(score_pairs(h_s, d5))$n))
})

test_that("streaming win statistics equal the stored ones", {
  pw <- score_pairs(h_s, dat_s)
  a <- win_stats(pw, dat_s$year)
  b <- win_stats(h_s, dat_s$year, data = dat_s)
  expect_equal(a$table, b$table, tolerance = 1e-12)
  expect_equal(a$by_level, b$by_level, tolerance = 1e-12)
  expect_equal(c(a$p_win, a$p_loss, a$p_tie), c(b$p_win, b$p_loss, b$p_tie))
  expect_error(win_stats(h_s, dat_s$year), "data")
})

test_that("streaming calibration agrees with the stored version", {
  a <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stored")
  b <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream")
  ca <- btd_calibration(a); cb <- btd_calibration(b)
  expect_equal(ca$n, cb$n)
  expect_equal(ca$obs_win_i, cb$obs_win_i, tolerance = 1e-8)
  expect_equal(ca$pred_win_i, cb$pred_win_i, tolerance = 1e-6)
  ba <- btd_calibration(a, by = "bin", n_bins = 5)
  bb <- btd_calibration(b, by = "bin", n_bins = 5)
  expect_equal(sum(bb$n), a$n_pairs)
  expect_lte(nrow(bb), 5L)
  expect_equal(bb$obs_win_i + bb$obs_tie + bb$obs_win_j, rep(1, nrow(bb)), tolerance = 1e-10)
  # overall tie rate: observed equals model-implied (score equation for delta)
  expect_equal(sum(bb$n * bb$obs_tie) / sum(bb$n), sum(bb$n * bb$pred_tie) / sum(bb$n),
               tolerance = 1e-6)
})

test_that("streaming works in Bradley-Terry mode (no ties)", {
  set.seed(9)
  d <- data.frame(x = rnorm(60)); d$y <- d$x + rnorm(60)
  h <- hierarchy(lvl_continuous("y", better = "higher"))
  expect_message(a <- btd(~ x, d, hierarchy = h, storage = "stored"), "Bradley-Terry")
  expect_message(b <- btd(~ x, d, hierarchy = h, storage = "stream"), "Bradley-Terry")
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_equal(a$vcov, b$vcov, tolerance = 1e-5)
})

test_that("streaming handles censored time-to-event levels", {
  set.seed(4)
  n <- 80
  d <- data.frame(x = rnorm(n), t0 = rexp(n), c0 = runif(n, 0.3, 2))
  d$time <- pmin(d$t0, d$c0); d$st <- as.integer(d$t0 <= d$c0)
  d$score <- sample(0:2, n, TRUE)
  h <- hierarchy(lvl_tte("time", status = "st", better = "later"),
                 lvl_ordinal("score", better = "lower"))
  a <- btd(~ x, d, hierarchy = h, storage = "stored")
  b <- btd(~ x, d, hierarchy = h, storage = "stream")
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_equal(b$stream$uninf, sum(a$pairs$uninf))
})

test_that("compact storage equals streaming bitwise and the stored fit numerically", {
  a <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stored")
  b <- btd(fml_s, dat_s, hierarchy = h_s, storage = "stream")
  c1 <- btd(fml_s, dat_s, hierarchy = h_s, storage = "compact")
  expect_equal(c1$storage, "compact")
  expect_null(c1$pairs)
  expect_null(c1$stream$codes)                    # the codes are not kept in the fit
  expect_identical(b$coefficients, c1$coefficients)
  expect_identical(b$vcov, c1$vcov)
  expect_identical(b$loglik, c1$loglik)
  expect_equal(coef(a), coef(c1), tolerance = 1e-7)
  expect_equal(a$vcov, c1$vcov, tolerance = 1e-5)
  expect_equal(a$n_pairs, c1$n_pairs); expect_equal(a$n_ties, c1$n_ties)
  expect_output(print(summary(c1)), "compact")
  # thread invariance
  c3 <- btd(fml_s, dat_s, hierarchy = h_s, storage = "compact",
            control = btd_control(threads = 3))
  expect_identical(c1$coefficients, c3$coefficients)
  expect_identical(c1$vcov, c3$vcov)
  # calibration and other summaries work from a compact fit
  expect_equal(btd_calibration(a)$n, btd_calibration(c1)$n)
})

test_that("compact storage: codes, counts and argument checks", {
  N <- 37
  d <- dat_s[seq_len(N), ]
  comp <- hpsreg:::.compile_hierarchy(h_s, d)
  specs <- hpsreg:::.stream_specs(comp)
  cb <- hpsreg:::cpp_compact_build(specs, comp$V, 32L, 1L)
  pw <- score_pairs(h_s, d)
  Kleaf <- length(h_s$leaves)
  cd <- as.integer(cb$codes); lev <- cd %/% 3L; w <- cd %% 3L - 1L
  # codes follow the all-pairs enumeration order and agree with score_pairs()
  expect_equal(w, pw$W)
  leaf <- if (is.null(pw$leaf)) pw$level else pw$leaf
  expect_equal(as.integer(lev) + 1L, leaf)
  expect_equal(cb$counts, as.numeric(tabulate(leaf, Kleaf + 1L)))
  expect_equal(length(cb$codes), N * (N - 1) / 2)
  # argument checks
  expect_error(btd(~ age, dat_s, pairs = score_pairs(h_s, dat_s), storage = "compact"),
               "hierarchy")
  expect_error(btd(~ age, dat_s, hierarchy = h_s, storage = "compact", engine = "r"),
               "engine")
  f <- function(data, i, j) sign(data$age[i] - data$age[j])
  hc <- hierarchy(lvl_custom(f), lvl_binary("hct_ci", favorable = "0-2"))
  expect_error(btd(~ donor_type, dat_s, hierarchy = hc, storage = "compact"), "lvl_custom")
  # more than 84 comparisons are rejected
  big <- do.call(hierarchy, lapply(seq_len(85), function(k)
    lvl_binary("hct_ci", favorable = "0-2", name = paste0("l", k))))
  expect_error(btd(~ age, dat_s, hierarchy = big, storage = "compact"), "84")
})
