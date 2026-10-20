# hpsreg

Regression analysis of **hierarchically ordered clinical outcomes** via pairwise comparisons (Bradley-Terry-Davidson Model).

1. **Design** a hierarchy of outcomes (`hierarchy()`, `lvl_binary()`, `lvl_ordinal()`,
   `lvl_continuous()`, `lvl_tte()`, `lvl_custom()`, `lvl_group()`), with tie margins and
   eligibility guards (`when =`).
2. **Score** all pairs of subjects (`score_pairs()`): the first level that separates a pair
   decides it. Compiled kernel with early exit; `resolution()` shows where pairs are decided.
3. **Estimate**: covariate-adjusted Davidson (Bradley-Terry with ties) model with U-statistic
   sandwich standard errors (`btd()`); unadjusted net benefit / win ratio / win odds
   (`win_stats()`); checks (`btd_calibration()`); predictions (`predict()`, `btd_pair_probs()`).

## Quick start (from scratch)

Define the hierarchy from the columns of your data frame, score all pairs, and fit.
Levels are listed from highest to lowest priority.

```r
library(hpsreg)
dat <- simulate_hps(N = 613, seed = 1)          # stem cell transplant example data

h <- hierarchy(
  lvl_tte("death_day",   better = "later", name = "survival"),    # NA = no death by day 180
  lvl_tte("relapse_day", better = "later", name = "relapse"),
  lvl_group("GVHD",
    lvl_ordinal("gvhd_grade", better = "lower", name = "grade"),
    lvl_binary("gvhd_resolved", favorable = 1, when = "gvhd_grade > 0", name = "resolution"),
    lvl_continuous("gvhd_onset_day", better = "higher", when = "gvhd_grade > 0", name = "onset"))
)
pw  <- score_pairs(h, dat)                       # win / loss / tie for every pair
resolution(pw)                                   # where pairs are decided
fit <- btd(~ donor_type + age + dri, data = dat, pairs = pw)
summary(fit)
```

## Two workflows

**1. Pre-computed pairs.** Score once, then reuse the pairs for diagnostics, win statistics
and several covariate models. Best for exploration, because scoring is the only step that
depends on the hierarchy.

```r
pw   <- score_pairs(h, dat)
win_stats(pw, dat$year)         # unadjusted net benefit / win ratio / win odds
fit1 <- btd(~ donor_type + age, data = dat, pairs = pw)
fit2 <- btd(~ donor_type + age + dri, data = dat, pairs = pw)   # no rescoring
```

**2. Inline scoring.** Hand `btd()` the hierarchy and it scores the pairs itself. The pairs
stay available as `fit$pairs`, so `btd_calibration(fit)` still works.

```r
fit <- btd(~ donor_type + age + dri, data = dat, hierarchy = h)
```

In both workflows `data` must be the same data frame (same rows, same order) that is scored.
Changing the hierarchy (reordering levels, adding a margin) means scoring again.

Worked examples (installed under `system.file("examples", package = "hpsreg")`): `hps-example.Rmd` (render with `rmarkdown::render()`) and `hps_simulation_study.R`.

## Install

```sh
R CMD INSTALL hpsreg_0.1.0.tar.gz     # needs a C++ toolchain; Rcpp and RcppEigen
```

## Large N: compact and streaming storage

All N(N-1)/2 pairs of N = 10,000 patients are 50 million pairs. `btd()` has four storage modes
(`storage = ` argument; `"auto"` picks one):

| mode | memory | speed | notes |
|---|---|---|---|
| `"stored"` | about 20 bytes per pair (more at peak) | fastest for small N | pairs kept as `fit$pairs`; needed for `lvl_custom()` levels and pre-computed `pairs` |
| `"compact"` | 1 byte per pair | about 30% faster than streaming | pairs scored once; about 50 MB at N = 10,000 |
| `"stream"` | O(N p), about 90 MB at N = 10,000 | slowest | pairs re-scored in every pass; nothing proportional to the number of pairs is stored |

`"auto"` uses stored pairs up to about N = 2200, compact up to about N = 44,000 and streaming
beyond that. Compact and streaming need a `hierarchy` of built-in levels (at most 84 comparisons
for compact) and cannot take pre-computed pairs.

```r
fit <- btd(~ donor_type + age, data = dat, hierarchy = h)                       # auto
fit <- btd(~ donor_type + age, data = dat, hierarchy = h, storage = "compact")  # explicit
resolution(h, data = dat)                     # counts per level, nothing stored
win_stats(h, dat$arm, data = dat)             # unadjusted win statistics, nothing stored
btd_calibration(fit)                          # works on compact and streamed fits
```

Optional multithreading (OpenMP; one thread by default, identical results for any thread
count): `btd(..., control = btd_control(threads = 4))` or `options(hpsreg.threads = 4)`.
`openmp_threads()` returns 0 if your build has no OpenMP (for example the default macOS
toolchain); everything then runs single-threaded.

## Scope and limitations

* Time-to-event levels: complete follow-up, or right censoring with the Gehan rule
  (uninformative pairs pass to the next level and are flagged). IPCW (Peron) scoring is not
  implemented; the rule is biased toward the null under heavy censoring.
* Pairs are all N(N-1)/2 pairs (or a supplied subset). Sandwich variances are exact in the
  U-statistic sense only for the complete set. Stored pairs cost about 20 bytes each (more at
  peak); `max_pairs` guards `score_pairs()`. Streaming (`storage = "stream"`) avoids storing pairs.
* Compact and streaming modes need a `hierarchy` made of built-in levels (no `lvl_custom()`) and do
  not accept pre-computed or sampled pairs.
* Cross-arm-only pair sets with matching variance theory are not implemented.
* One tie parameter is shared by all levels; use `btd_calibration()` to check the single
  latent-scale assumption. Thresholds and censoring can make pairwise relations intransitive.
* The maintainer e-mail in `DESCRIPTION` is a placeholder.
