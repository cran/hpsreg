test_that("simulate_hps reproduces the original simulator exactly", {
  beta <- hps_beta()
  a <- simulate_hps(beta, N = 150, seed = 123)
  b <- simulate_cohort(beta, N = 150, seed = 123)
  expect_identical(a, b)
})

test_that("hps_hierarchy reproduces the original pairwise scores exactly", {
  dat <- simulate_hps(N = 140, seed = 99)
  pw  <- score_pairs(hps_hierarchy(), dat)
  ref <- build_pairwise_fast(dat)
  expect_identical(pw$id_i, ref$id_i)
  expect_identical(pw$id_j, ref$id_j)
  expect_identical(pw$W, ref$W)
  expect_identical(pw$level, ref$level)
  expect_true(all(pw$uninf == 0L))
  rs <- resolution(pw)
  expect_equal(rs$n, as.numeric(tabulate(ref$level, 5L)), ignore_attr = TRUE)
})

test_that("a user-supplied pair subset gives the same scores as the full set", {
  dat <- simulate_hps(N = 60, seed = 5)
  full <- score_pairs(hps_hierarchy(), dat)
  sp <- sample_pairs(60, 300)
  sub <- score_pairs(hps_hierarchy(), dat, pairs = sp)
  key <- function(i, j) paste(i, j)
  m <- match(key(sub$id_i, sub$id_j), key(full$id_i, full$id_j))
  expect_identical(sub$W, full$W[m])
  expect_identical(sub$level, full$level[m])
})

test_that("Cox comparator runs", {
  skip_if_not_installed("survival")
  dat <- simulate_hps(N = 200, seed = 4)
  r <- hps_cox_grfs(dat, c("donor_type", "age"))
  expect_named(r, c("fit", "coef", "se"))
  expect_equal(length(r$coef), 3L)
})
