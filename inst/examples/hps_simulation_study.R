## Simulation study for the stem cell transplant example, written with the hpsreg API.
## Compares the hierarchical (BTD) analysis with the composite m-GRFS Cox model:
## bias / coverage relative to a large-sample BTD fit, and rejection rates.
##
##   source(system.file("examples", "hps_simulation_study.R", package = "hpsreg"))
##   res <- run_hps_study(n_rep = 20)         # increase n_rep for real use
##   summarise_hps_study(res)

library(hpsreg)

hps_covars <- c("donor_type", "age", "dri", "hct_ci", "sex_mismatch", "year",
                "conditioning", "cell_source")
hps_hier <- hps_hierarchy()

one_rep <- function(rep, beta = hps_beta(), N = 613, seed0 = 20241004) {
  dat <- simulate_hps(beta, N = N, seed = seed0 + rep)
  pw  <- score_pairs(hps_hier, dat)
  fit <- btd(reformulate(hps_covars), dat, pairs = pw)
  cx  <- hps_cox_grfs(dat, hps_covars)
  data.frame(rep = rep, term = names(coef(fit)),
             btd = unname(coef(fit)), btd_se = sqrt(diag(vcov(fit))),
             cox = -unname(cx$coef), cox_se = unname(cx$se))
}

run_hps_study <- function(n_rep = 20, N = 613, N_big = 3000, ...) {
  message("Population (large-N) BTD fit, N = ", N_big)
  big <- simulate_hps(N = N_big, seed = 1)
  truth <- coef(btd(reformulate(hps_covars), big, hierarchy = hps_hier))
  res <- do.call(rbind, lapply(seq_len(n_rep), one_rep, N = N, ...))
  res$truth <- truth[res$term]
  res
}

summarise_hps_study <- function(res, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  do.call(rbind, lapply(split(res, res$term), function(d) {
    data.frame(term = d$term[1], truth = d$truth[1],
               bias = mean(d$btd) - d$truth[1], emp_sd = sd(d$btd),
               mean_se = mean(d$btd_se),
               coverage = mean(abs(d$btd - d$truth) <= z * d$btd_se),
               power_btd = mean(abs(d$btd / d$btd_se) > z),
               power_cox = mean(abs(d$cox / d$cox_se) > z))
  }))
}
