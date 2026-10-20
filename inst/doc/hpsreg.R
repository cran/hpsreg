## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4, out.width = "90%")
set.seed(1)

## -----------------------------------------------------------------------------
library(hpsreg)

## -----------------------------------------------------------------------------
dat <- simulate_hps(N = 200, seed = 11)
str(dat, list.len = 6)

## -----------------------------------------------------------------------------
h <- hierarchy(
  lvl_tte("death_day",   better = "later", name = "survival"),
  lvl_tte("relapse_day", better = "later", name = "relapse"),
  lvl_ordinal("gvhd_grade", better = "lower", name = "GVHD grade"),
  lvl_binary("gvhd_resolved", favorable = 1,
             when = "gvhd_grade > 0", name = "GVHD resolution")
)
h

## -----------------------------------------------------------------------------
pw <- score_pairs(h, dat)
pw
resolution(pw)

## -----------------------------------------------------------------------------
win_stats(pw, dat$year)

## -----------------------------------------------------------------------------
win_stats(h, dat$year, data = dat)

## -----------------------------------------------------------------------------
fit <- btd(~ donor_type + age + dri, data = dat, pairs = pw)
summary(fit)

## -----------------------------------------------------------------------------
fit2 <- btd(~ donor_type + age + dri, data = dat, hierarchy = h)
all.equal(coef(fit), coef(fit2))

## -----------------------------------------------------------------------------
cbind(sandwich = sqrt(diag(vcov(fit))),
      model    = sqrt(diag(vcov(fit, model = TRUE))))

## -----------------------------------------------------------------------------
predict(fit, dat[1:3, ])
btd_pair_probs(fit, dat[1:3, ], dat[4:6, ])

## ----fig.height = 3.5---------------------------------------------------------
btd_calibration(fit, by = "level")
plot(btd_calibration(fit, by = "level"))

## ----eval = FALSE-------------------------------------------------------------
# fit <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact")

## -----------------------------------------------------------------------------
openmp_threads()

## ----eval = FALSE-------------------------------------------------------------
# fit_mt <- btd(~ donor_type + age, data = dat, hierarchy = h,
#               control = btd_control(threads = 4))
# 
# options(hpsreg.threads = 4)          # applies to every subsequent call
# fit_mt2 <- btd(~ donor_type + age, data = dat, hierarchy = h)
# all.equal(coef(fit_mt), coef(fit_mt2))

