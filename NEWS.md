# hpsreg 1.0.4

Changes since the CRAN release 1.0.0 (versions 1.0.1 to 1.0.3 were not released). Fitted
coefficients and standard errors are unchanged except in the cases listed under "Behaviour
changes"; the fitting engine is about 2.5 times faster (N = 10,000, compact storage, one
thread: 46 s in 1.0.0, 18 s now, identical estimates to about 1e-16, same memory).

## Behaviour changes that can affect existing code

* **`btd_calibration(by = "level")` columns are renamed.** `obs_win_i` and `pred_win_i` are now
  `obs_win_higher` and `pred_win_higher`; code that refers to the old names must be updated.
  See the calibration item below for why the quantity changed.
* **Events and censorings at the same time** are scored differently (see below). Results change
  for data with coincident event and censoring times; `lvl_tte(same_time = "uninformative")`
  restores the 1.0.0 behaviour.
* **Events exactly at a finite `horizon`** are now events (see below). Results with
  `horizon = Inf` (the default) are unchanged.
* `fit$n_grad_evals` is now always 0 and `fit$n_full_evals` equals one initial pass plus one
  pass per trial (see "Faster fitting").

## Bug fixes

* **Thread counts are capped.** A very large thread request (for example
  `btd_control(threads = 1e5)` or `options(hpsreg.threads = 1e5)`) crashed R with a C stack
  overflow. Requests are now capped at 32 (the number of chunks the pairs are split into) and
  at the number of processors; a non-numeric request is an error.
* **By-level calibration is oriented by the fitted score.** `btd_calibration(by = "level")`
  compared the share of pairs won by the *first-listed* subject with its prediction. Row order
  is arbitrary, so that share depended on how the data were sorted and not on the model (a
  pure-noise model gave the same "observed" column). Each pair is now oriented so that the
  higher-scored subject comes first; the columns are `obs_win_higher` and `pred_win_higher`
  (renamed from `obs_win_i` and `pred_win_i`), and the plot marks 0.5 (no predictive value).
* **Pre-computed pairs are checked against `data`.** `score_pairs()` stores a fingerprint of
  the outcome data, and `btd(pairs = )` stops if `data` has different outcome values or a
  different row order. Previously a reordered data frame was accepted silently and gave wrong
  estimates.
* **(Quasi-)separation is reported again.** `btd()` warns when the fitted latent scores imply
  win probabilities indistinguishable from 1 (log-odds range above 60 for a converged fit, or
  above 40 alongside a convergence failure), and returns the range as `latent_range`.
* **Events and censorings at the same time.** With a `status` column, a pair in which one
  subject has the event and the other is censored at the same time is now scored with the usual
  convention: the censored subject outlived the event. Previously such pairs were treated as
  uninformative and passed to the next level. `lvl_tte(same_time = "uninformative")` restores
  the old behaviour. This changes results for data with coincident event and censoring times.
* **Events exactly at the horizon.** With a finite `horizon`, an event exactly at the horizon
  was tied with being event-free. It is now an event, consistent with `hps_cox_grfs()`
  (which counts events with time <= horizon). Subjects followed beyond the horizon, or censored
  exactly at it, are event-free through the horizon and beat every event at or before it,
  whatever the tie margin. Results with `horizon = Inf` (the default) are unchanged.
* `lvl_binary()` warns when the favourable value never occurs, and `lvl_ordinal()` warns that
  `order` is ignored for a numeric column.
* `summary()` notes that the z test on the `delta` row is of delta = 1, not of whether ties occur.

## Faster fitting (engine; no change to estimates)

* The pair-level contribution to the negative Hessian is accumulated per subject through the
  identity `sum_{i<j} f_ij (x_i - x_j)(x_i - x_j)' = sum_i D_i x_i x_i' - sum_i x_i m_i'`,
  where `D_i` and `m_i` are per-subject scalar and vector sums. This reduces the per-pair cost
  from O(p^2) to O(p), with one O(N p^2) reassembly at the end. Covariate columns are centered
  before accumulation and the Hessian block is symmetrised, so the identity does not lose
  precision when covariates have a large mean.
* Each Levenberg-Marquardt trial is a full pass (log-likelihood, gradient and Hessian). The
  accepted trial is cached and reused as the start of the next iteration, which removes the
  extra pass per iteration and the final pass after the loop. `fit$n_grad_evals` is therefore 0
  and `fit$n_full_evals` is `1 + sum(fit$trials_per_iter)`.

## Other changes

* Documentation: corrected the vignette's description of the era win statistics (the levels
  pull in opposite directions; there is no overall advantage), its advice on checking the tie
  parameter, and the storage-mode memory figures; noted that `AIC()` and `BIC()` from the
  pairwise likelihood are not valid model-selection criteria.

# hpsreg 1.0.0

* Initial CRAN release.