sc <- function(h, d) score_pairs(h, d)
pw_df <- function(p) as.data.frame(p)
get_w <- function(p, i, j) p$W[p$id_i == i & p$id_j == j]

test_that("binary levels and NA rules", {
  d <- data.frame(x = c(1, 0, NA))
  p <- sc(hierarchy(lvl_binary("x", favorable = 1)), d)
  expect_equal(get_w(p, 1, 2), 1L)
  expect_equal(get_w(p, 1, 3), 0L)          # NA ties by default
  p2 <- sc(hierarchy(lvl_binary("x", favorable = 1, na = "worst")), d)
  expect_equal(get_w(p2, 1, 3), 1L)
  expect_equal(get_w(p2, 2, 3), 1L)
  p3 <- sc(hierarchy(lvl_binary("x", favorable = 1, na = "best")), d)
  expect_equal(get_w(p3, 1, 3), -1L)
  # favourable = 0 flips direction; factors / logicals work
  p4 <- sc(hierarchy(lvl_binary("x", favorable = 0)), d)
  expect_equal(get_w(p4, 1, 2), -1L)
  d5 <- data.frame(x = factor(c("yes", "no", "yes")), y = c(TRUE, FALSE, TRUE))
  expect_equal(get_w(sc(hierarchy(lvl_binary("x", favorable = "yes")), d5), 1, 2), 1L)
  expect_equal(get_w(sc(hierarchy(lvl_binary("y", favorable = TRUE)), d5), 1, 3), 0L)
})

test_that("continuous levels honour direction and tie margin", {
  d <- data.frame(v = c(10, 12, 20))
  p <- sc(hierarchy(lvl_continuous("v", better = "higher", threshold = 3)), d)
  expect_equal(p$W, c(0L, -1L, -1L))          # (1,2) within margin
  p <- sc(hierarchy(lvl_continuous("v", better = "lower")), d)
  expect_equal(p$W, c(1L, 1L, 1L))
  # margin is inclusive: |diff| == threshold is a tie
  p <- sc(hierarchy(lvl_continuous("v", better = "higher", threshold = 2)), d)
  expect_equal(get_w(p, 1, 2), 0L)
})

test_that("ordinal levels use factor order or explicit order", {
  d <- data.frame(g = factor(c("mild", "severe", "none"),
                             levels = c("none", "mild", "severe")),
                  s = c("b", "c", "a"), stringsAsFactors = FALSE)
  p <- sc(hierarchy(lvl_ordinal("g", better = "lower")), d)
  expect_equal(p$W, c(1L, -1L, -1L))
  expect_equal(get_w(p, 1, 2), 1L)           # mild beats severe
  expect_equal(get_w(p, 1, 3), -1L)          # none beats mild
  p <- sc(hierarchy(lvl_ordinal("s", better = "higher", order = c("c", "b", "a"))), d)
  expect_equal(get_w(p, 1, 2), 1L)           # in the given order "b" ranks above "c"
})

test_that("time-to-event without status: NA means event-free", {
  d <- data.frame(t = c(5, NA, 20, NA))
  p <- sc(hierarchy(lvl_tte("t", better = "later")), d)
  expect_equal(get_w(p, 1, 2), -1L)          # event vs none: none better
  expect_equal(get_w(p, 1, 3), -1L)          # later event better
  expect_equal(get_w(p, 2, 4), 0L)           # both event-free
  expect_equal(get_w(p, 3, 4), -1L)
  pe <- sc(hierarchy(lvl_tte("t", better = "earlier")), d)
  expect_equal(get_w(pe, 1, 3), 1L)
  expect_equal(get_w(pe, 1, 2), 1L)
})

test_that("time-to-event with censoring uses the Gehan rule", {
  d <- data.frame(t = c(5, 10, 8, 20), s = c(1, 1, 0, 0))
  p <- sc(hierarchy(lvl_tte("t", status = "s", better = "later")), d)
  expect_equal(get_w(p, 1, 2), -1L)          # both exact
  expect_equal(get_w(p, 1, 3), -1L)          # event at 5, censored at 8: certain
  expect_equal(get_w(p, 2, 3), 0L)           # event at 10, censored at 8: unknown
  expect_equal(get_w(p, 1, 4), -1L)
  expect_equal(get_w(p, 2, 4), -1L)
  expect_equal(get_w(p, 3, 4), 0L)           # both censored
  unf <- function(i, j) p$uninf[p$id_i == i & p$id_j == j]
  expect_equal(c(unf(2, 3), unf(3, 4), unf(1, 2)), c(1L, 1L, 0L))
  # earlier-is-better flips the decided pairs only
  pe <- sc(hierarchy(lvl_tte("t", status = "s", better = "earlier")), d)
  expect_equal(get_w(pe, 1, 3), 1L)
  expect_equal(get_w(pe, 2, 3), 0L)
  # an uninformative pair is passed down and flagged
  d$z <- c(0, 1, 0, 0)
  p2 <- sc(hierarchy(lvl_tte("t", status = "s", better = "later"),
                     lvl_binary("z", favorable = 1)), d)
  expect_equal(get_w(p2, 2, 3), 1L)
  expect_equal(p2$uninf[p2$id_i == 2 & p2$id_j == 3], 1L)
  expect_equal(p2$level[p2$id_i == 2 & p2$id_j == 3], 2L)
})

test_that("time-to-event margin and horizon", {
  d <- data.frame(t = c(50, 200, 55, NA), s = c(1, 1, 1, 0))
  p <- sc(hierarchy(lvl_tte("t", status = "s", better = "later", horizon = 100)), d)
  expect_equal(get_w(p, 1, 2), -1L)          # 200 truncated to 100 (complete)
  expect_equal(get_w(p, 1, 4), 0L)           # missing time: neutral
  p <- sc(hierarchy(lvl_tte("t", status = "s", better = "later", threshold = 10)), d)
  expect_equal(get_w(p, 1, 3), 0L)           # 50 vs 55 within margin
  expect_equal(get_w(p, 1, 2), -1L)
  expect_error(sc(hierarchy(lvl_tte("t", status = "t", better = "later")), d), "0/1")
})

test_that("eligibility guards make a level neutral for ineligible pairs", {
  d <- data.frame(grade = c(1, 1, 0, 0), res = c(1, 0, NA, NA))
  h <- hierarchy(lvl_ordinal("grade", better = "lower"),
                 lvl_binary("res", favorable = 1, when = grade > 0))
  p <- sc(h, d)
  expect_equal(get_w(p, 1, 2), 1L)
  expect_equal(p$level[p$id_i == 1 & p$id_j == 2], 2L)
  expect_equal(get_w(p, 3, 4), 0L)
  expect_equal(get_w(p, 1, 3), -1L)
  # formula / string / expression forms of `when` agree
  h2 <- hierarchy(lvl_ordinal("grade", better = "lower"),
                  lvl_binary("res", favorable = 1, when = ~ grade > 0))
  h3 <- hierarchy(lvl_ordinal("grade", better = "lower"),
                  lvl_binary("res", favorable = 1, when = "grade > 0"))
  expect_equal(sc(h2, d)$W, p$W)
  expect_equal(sc(h3, d)$W, p$W)
  # guards can use variables from the calling environment
  cut <- 0
  h4 <- hierarchy(lvl_ordinal("grade", better = "lower"),
                  lvl_binary("res", favorable = 1, when = grade > cut))
  expect_equal(sc(h4, d)$W, p$W)
})

test_that("custom levels interleave with built-in ones", {
  set.seed(3)
  d <- data.frame(a = sample(0:1, 25, TRUE), b = rnorm(25), c = sample(1:3, 25, TRUE))
  f <- function(data, i, j) sign(data$b[i] - data$b[j])
  h <- hierarchy(lvl_binary("a", favorable = 1), lvl_custom(f, "b"),
                 lvl_ordinal("c", better = "lower"))
  p <- sc(h, d)
  i <- p$id_i; j <- p$id_j
  w1 <- sign(d$a[i] - d$a[j]); w2 <- sign(d$b[i] - d$b[j]); w3 <- sign(d$c[j] - d$c[i])
  exp_w <- ifelse(w1 != 0, w1, ifelse(w2 != 0, w2, w3))
  expect_equal(p$W, as.integer(exp_w))
  expect_equal(p$level, ifelse(w1 != 0, 1L, ifelse(w2 != 0, 2L, ifelse(w3 != 0, 3L, 4L))))
  # NA from a custom level is uninformative and passed down
  g <- function(data, i, j) ifelse(data$a[i] == data$a[j], NA, 0)
  p2 <- sc(hierarchy(lvl_custom(g), lvl_ordinal("c", better = "lower")), d)
  expect_true(all(p2$uninf[d$a[p2$id_i] == d$a[p2$id_j]] == 1L))
  expect_error(sc(hierarchy(lvl_custom(function(data, i, j) rep(2, length(i)))), d),
               "\\{-1, 0, 1")
  # a custom first level is evaluated on all pairs
  p3 <- sc(hierarchy(lvl_custom(f), lvl_binary("a", favorable = 1)), d)
  wb <- sign(d$b[p3$id_i] - d$b[p3$id_j])
  expect_equal(p3$W[wb != 0], as.integer(wb[wb != 0]))
})

test_that("groups report as one level and leaves are available", {
  d <- data.frame(g = c(1, 1, 0), r = c(1, 0, NA), t = c(5, 9, 1))
  h <- hierarchy(lvl_group("G", lvl_ordinal("g", better = "lower"),
                           lvl_binary("r", favorable = 1, when = g > 0)),
                 lvl_continuous("t", better = "higher"))
  p <- sc(h, d)
  expect_equal(p$level, c(1L, 1L, 1L))
  expect_equal(p$leaf, c(2L, 1L, 1L))
  rs <- resolution(p); rl <- resolution(p, by = "leaf")
  expect_equal(rs$n, c(3L, 0L, 0L)); expect_equal(sum(rl$n), 3L)
  expect_equal(as.character(rl$level[1:3]), c("G: g", "G: r", "t"))
  expect_equal(nrow(rs), 3L)
})

test_that("pair sets: explicit, sampled, and index conversion", {
  N <- 9
  all <- t(combn(N, 2))
  k <- seq_len(nrow(all))
  r <- hpsreg:::cpp_index_to_pair(as.numeric(k), N)
  expect_equal(cbind(r$i, r$j), unname(all))
  expect_error(hpsreg:::cpp_index_to_pair(nrow(all) + 1, N), "range")
  big <- hpsreg:::cpp_index_to_pair(c(1, 4999950000), 100000L)   # last pair of 1e5 subjects
  expect_equal(c(big$i[2], big$j[2]), c(99999L, 100000L))
  set.seed(1)
  sp <- sample_pairs(30, 50)
  expect_equal(nrow(unique(sp)), 50L)
  expect_true(all(sp[, 1] < sp[, 2]))
  expect_error(sample_pairs(5, 11), "exceeds")

  d <- data.frame(v = rnorm(30))
  h <- hierarchy(lvl_continuous("v", better = "higher"))
  p <- score_pairs(h, d, pairs = sp)
  expect_false(p$complete)
  expect_equal(p$W, as.integer(sign(d$v[sp[, 1]] - d$v[sp[, 2]])))
  expect_error(score_pairs(h, d, pairs = cbind(1, 1)), "distinct")
  expect_error(score_pairs(h, d, pairs = cbind(1, 99)), "valid")
  expect_error(score_pairs(h, d, max_pairs = 10), "max_pairs")
})

test_that("print and as.data.frame methods work", {
  d <- data.frame(v = c(1, 2, 2, 3))
  p <- score_pairs(hierarchy(lvl_continuous("v", better = "higher")), d)
  expect_output(print(p), "pairs among 4")
  expect_s3_class(as.data.frame(p), "data.frame")
  expect_equal(length(p), 6L)
  expect_equal(sum(p$W == 0L), 1L)
})
