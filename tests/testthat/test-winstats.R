test_that("win statistics match a brute-force calculation", {
  set.seed(21)
  N <- 40
  d <- data.frame(g = rep(c("c", "t"), c(18, 22)), y = rnorm(N))
  d$y[d$g == "t"] <- d$y[d$g == "t"] + 0.5
  d$z <- sample(0:1, N, TRUE)
  h <- hierarchy(lvl_continuous("y", better = "higher", threshold = 0.4),
                 lvl_binary("z", favorable = 1))
  pw <- score_pairs(h, d)
  ws <- win_stats(pw, d$g, treatment = "t")

  tr <- which(d$g == "t"); co <- which(d$g == "c")
  M <- matrix(0L, length(tr), length(co))
  for (a in seq_along(tr)) for (b in seq_along(co)) {
    i <- tr[a]; j <- co[b]
    lo <- min(i, j); hi <- max(i, j)
    w <- pw$W[pw$id_i == lo & pw$id_j == hi]
    M[a, b] <- if (i < j) w else -w            # +1 = treatment wins
  }
  pw_ <- mean(M == 1L); pl_ <- mean(M == -1L)
  expect_equal(ws$p_win, pw_); expect_equal(ws$p_loss, pl_)
  expect_equal(ws$p_tie, 1 - pw_ - pl_)

  Wa <- rowMeans(M == 1L); La <- rowMeans(M == -1L)
  Wb <- colMeans(M == 1L); Lb <- colMeans(M == -1L)
  n1 <- length(tr); n0 <- length(co)
  vw <- var(Wa) / n1 + var(Wb) / n0
  vl <- var(La) / n1 + var(Lb) / n0
  cw <- cov(Wa, La) / n1 + cov(Wb, Lb) / n0
  nb <- pw_ - pl_
  tab <- ws$table
  expect_equal(tab$estimate[1], nb)
  expect_equal(tab$se_scale[1], sqrt(vw + vl - 2 * cw))
  expect_equal(tab$estimate[2], pw_ / pl_)
  expect_equal(tab$se_scale[2], sqrt(vw / pw_^2 + vl / pl_^2 - 2 * cw / (pw_ * pl_)))
  expect_equal(tab$estimate[3], (pw_ + (1 - pw_ - pl_) / 2) / (pl_ + (1 - pw_ - pl_) / 2))
  expect_true(all(tab$lower < tab$estimate & tab$estimate < tab$upper))
  expect_equal(sum(ws$by_level$treatment_wins), pw_)
  expect_equal(sum(ws$by_level$treatment_losses), pl_)
  expect_equal(sum(ws$by_level$tied), ws$p_tie)
  expect_output(print(ws), "Win statistics")
})

test_that("direction and defaults are correct", {
  set.seed(2)
  d <- data.frame(g = factor(rep(c("ctl", "trt"), each = 30)), y = c(rnorm(30), rnorm(30, 1)))
  pw <- score_pairs(hierarchy(lvl_continuous("y", better = "higher")), d)
  ws <- win_stats(pw, d$g)                      # treatment defaults to 2nd level
  expect_equal(ws$treatment, "trt")
  expect_gt(ws$table$estimate[1], 0)
  expect_gt(ws$table$estimate[2], 1)
  ws2 <- win_stats(pw, d$g, treatment = "ctl")
  expect_equal(ws2$table$estimate[1], -ws$table$estimate[1])
  expect_equal(ws2$table$estimate[2], 1 / ws$table$estimate[2])
  # lower-is-better flips the sign
  pw2 <- score_pairs(hierarchy(lvl_continuous("y", better = "lower")), d)
  expect_lt(win_stats(pw2, d$g)$table$estimate[1], 0)
})

test_that("win_stats validates its input", {
  d <- data.frame(g = c("a", "b", "c", "a", "b", "c"), y = rnorm(6))
  pw <- score_pairs(hierarchy(lvl_continuous("y", better = "higher")), d)
  expect_error(win_stats(pw, d$g), "exactly two")
  expect_error(win_stats(pw, d$g[-1]), "one value")
  expect_error(win_stats(pw, c("a", "b", NA, "a", "b", "a")), "missing")
  expect_error(win_stats(pw, rep(c("a", "b"), 3), treatment = "z"), "not a value")
  sp <- score_pairs(hierarchy(lvl_continuous("y", better = "higher")), d,
                    pairs = cbind(1:3, 4:6))
  expect_error(win_stats(sp, rep(c("a", "b"), 3)), "complete")
})
