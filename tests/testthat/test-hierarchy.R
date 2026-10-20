test_that("level constructors validate input", {
  expect_error(lvl_ordinal("x"), "better")
  expect_error(lvl_continuous("x"), "better")
  expect_error(lvl_tte("t"), "better")
  expect_error(lvl_ordinal("x", better = "up"), "better")
  expect_error(lvl_tte("t", better = "lower"), "better")
  expect_error(lvl_continuous("x", better = "higher", threshold = -1), "threshold")
  expect_error(lvl_binary(1), "var")
  expect_error(lvl_custom(1), "function")
  expect_error(lvl_tte("t", better = "later", horizon = 0), "horizon")
})

test_that("`na` is not partially matched to `name`", {
  l <- lvl_binary("x", name = "foo", na = "worst")
  expect_equal(l$name, "foo")
  expect_equal(l$na, "worst")
  expect_equal(lvl_tte("t", better = "later", name = "bar")$name, "bar")
})

test_that("hierarchy() validates and flattens", {
  expect_error(hierarchy(), "at least one")
  expect_error(hierarchy(1), "levels")
  expect_error(hierarchy(lvl_binary("a"), lvl_binary("a")), "unique")
  h <- hierarchy(lvl_binary("a"),
                 lvl_group("G", lvl_ordinal("b", better = "lower"),
                           lvl_continuous("c", better = "higher")),
                 lvl_binary("d"))
  expect_equal(length(h), 3L)
  expect_equal(length(h$leaves), 4L)
  expect_equal(h$leaf_group, c(1L, 2L, 2L, 3L))
  expect_true(h$grouped)
  expect_equal(h$leaf_names, c("a", "G: b", "G: c", "d"))
  expect_output(print(h), "hierarchy")
  # a list of levels is accepted
  expect_equal(length(hierarchy(list(lvl_binary("a"), lvl_binary("b")))), 2L)
})

test_that("group guards are combined with child guards", {
  h <- hierarchy(lvl_group("G", lvl_binary("r", when = g > 0), when = ok))
  expect_length(h$leaves[[1]]$when, 2L)
})

test_that("compile errors are informative", {
  d <- data.frame(x = 1:4)
  expect_error(score_pairs(hierarchy(lvl_binary("zz")), d), "zz")
  expect_error(score_pairs(hierarchy(lvl_binary("x", when = x)), d), "logical")
  expect_error(score_pairs(hierarchy(lvl_tte("x", status = "x2", better = "later")), d), "x2")
})
