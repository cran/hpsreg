## Reference implementations copied unchanged from the original study script
## (hps_sim_master.R). They are used only to verify that the package reproduces
## the original pairwise scores and simulated data exactly.

simulate_cohort <- function(beta_true_sim, N = 613, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  req <- c("intercept", "MUD", "HID", "DRI_int", "DRI_high",
           "DRI_vhigh", "sex_mismatch", "year_2011",
           "age_per_year", "hct_ci_high")
  if (!all(req %in% names(beta_true_sim)))
    stop("Missing names in beta_true_sim: ",
         paste(setdiff(req, names(beta_true_sim)), collapse = ", "))

  donor_type   <- sample(c("MRD", "MUD", "HID"), N,
                         prob = c(212, 251, 150) / N, replace = TRUE)
  age          <- round(pmin(pmax(rnorm(N, mean = 53, sd = 15), 18), 77))
  sex          <- sample(c("M", "F"), N, prob = c(0.56, 0.44), replace = TRUE)
  dri          <- sample(c("Low", "Intermediate", "High", "Very high"), N,
                         prob = c(90, 313, 171, 39) / N, replace = TRUE)
  hct_ci       <- sample(c("0-2", ">=3"), N, prob = c(0.61, 0.39), replace = TRUE)
  cell_source  <- sample(c("BM", "PBSC"), N, prob = c(118, 495) / N, replace = TRUE)
  conditioning <- sample(c("RIC", "MAC"), N, prob = c(302, 311) / N, replace = TRUE)
  year         <- sample(c("2005-2010", "2011-2015"), N,
                         prob = c(0.44, 0.56), replace = TRUE)
  sex_mismatch <- sample(c("F_to_M", "other"), N,
                         prob = c(0.18, 0.82), replace = TRUE)
  diagnosis    <- sample(c("AML", "ALL", "MDS/MPS/CML", "NHL/HD/CLL", "Other"), N,
                         prob = c(229, 82, 144, 142, 16) / N, replace = TRUE)

  eta_true <- beta_true_sim["intercept"] +
    beta_true_sim["MUD"]          * (donor_type == "MUD") +
    beta_true_sim["HID"]          * (donor_type == "HID") +
    beta_true_sim["DRI_int"]      * (dri == "Intermediate") +
    beta_true_sim["DRI_high"]     * (dri == "High") +
    beta_true_sim["DRI_vhigh"]    * (dri == "Very high") +
    beta_true_sim["sex_mismatch"] * (sex_mismatch == "F_to_M") +
    beta_true_sim["year_2011"]    * (year == "2011-2015") +
    beta_true_sim["age_per_year"] * (age - 53) +
    beta_true_sim["hct_ci_high"]  * (hct_ci == ">=3")
  eta_true <- eta_true - mean(eta_true)

  sigma_eps <- 0.8
  eta       <- eta_true + rnorm(N, 0, sigma_eps)
  eta_std   <- (eta - mean(eta)) / sd(eta)

  cal <- function(lp, target) {
    uniroot(function(a) mean(plogis(a + lp)) - target,
            interval = c(-20, 20))$root
  }

  lp_death <- -0.6 * eta_std +
    0.3 * (donor_type == "MUD") +
    0.6 * (donor_type == "HID") +
    0.3 * (hct_ci == ">=3")
  died_180 <- rbinom(N, 1, plogis(cal(lp_death, 0.12) + lp_death))
  death_day <- ifelse(died_180 == 1,
                      round(pmin(180, rgamma(N, shape = 2, scale = 40))), NA)

  lp_relapse <- -0.5 * eta_std +
    0.4 * (dri == "Intermediate") +
    0.8 * (dri == "High") +
    0.7 * (dri == "Very high")
  relapse_180 <- ifelse(died_180 == 1, 0,
                        rbinom(N, 1, plogis(cal(lp_relapse, 0.15) + lp_relapse)))
  relapse_day <- ifelse(relapse_180 == 1, round(runif(N, 30, 180)), NA)

  lp_gvhd <- -0.5 * eta_std +
     0.3 * (donor_type == "MUD") +
    -0.2 * (donor_type == "HID") +
     0.2 * (conditioning == "MAC")
  p_g1 <- plogis(cal(lp_gvhd, 0.50) + lp_gvhd)
  p_g2 <- plogis(cal(lp_gvhd, 0.20) + lp_gvhd)

  gvhd_grade <- rep(NA_integer_, N)
  alive      <- which(died_180 == 0)
  u          <- runif(length(alive))
  gvhd_grade[alive] <- ifelse(u < p_g2[alive], 2L,
                       ifelse(u < p_g1[alive], 1L, 0L))
  gvhd_resolved <- ifelse(is.na(gvhd_grade) | gvhd_grade == 0, NA,
                   ifelse(gvhd_grade == 1, rbinom(N, 1, 0.80),
                                          rbinom(N, 1, 0.50)))
  gvhd_onset_day <- ifelse(!is.na(gvhd_grade) & gvhd_grade > 0,
                           round(runif(N, 14, 120)), NA)

  lp_infection <- -0.4 * eta_std +
    0.2 * (hct_ci == ">=3") +
    0.2 * (conditioning == "MAC")
  p_i1 <- plogis(cal(lp_infection, 0.45) + lp_infection)
  p_i2 <- plogis(cal(lp_infection, 0.20) + lp_infection)

  infection_grade <- rep(NA_integer_, N)
  u <- runif(length(alive))
  infection_grade[alive] <- ifelse(u < p_i2[alive], 2L,
                            ifelse(u < p_i1[alive], 1L, 0L))
  infection_resolved <- ifelse(is.na(infection_grade) | infection_grade == 0, NA,
                        ifelse(infection_grade == 1, rbinom(N, 1, 0.85),
                                                      rbinom(N, 1, 0.60)))
  infection_onset_day <- ifelse(!is.na(infection_grade) & infection_grade > 0,
                                round(runif(N, 7, 180)), NA)

  data.frame(
    id                  = 1:N,
    donor_type          = factor(donor_type,   levels = c("MRD", "MUD", "HID")),
    age                 = age,
    sex                 = factor(sex,          levels = c("F", "M")),
    dri                 = factor(dri,          levels = c("Low", "Intermediate",
                                                         "High", "Very high")),
    hct_ci              = factor(hct_ci,       levels = c("0-2", ">=3")),
    cell_source         = factor(cell_source,  levels = c("BM", "PBSC")),
    conditioning        = factor(conditioning, levels = c("MAC", "RIC")),
    year                = factor(year,         levels = c("2005-2010", "2011-2015")),
    sex_mismatch        = factor(sex_mismatch, levels = c("other", "F_to_M")),
    diagnosis           = factor(diagnosis,    levels = c("AML", "ALL",
                                                          "MDS/MPS/CML",
                                                          "NHL/HD/CLL", "Other")),
    eta_true            = eta_true,
    eta                 = eta,
    alpha               = exp(eta),
    died_180            = died_180,
    death_day           = death_day,
    relapse_180         = relapse_180,
    relapse_day         = relapse_day,
    gvhd_grade          = gvhd_grade,
    gvhd_resolved       = gvhd_resolved,
    gvhd_onset_day      = gvhd_onset_day,
    infection_grade     = infection_grade,
    infection_resolved  = infection_resolved,
    infection_onset_day = infection_onset_day
  )
}


compare_survival_vec <- function(d1, d2) {
  a1 <- is.na(d1); a2 <- is.na(d2)
  w  <- integer(length(d1))
  w[ a1 & !a2] <-  1L; w[!a1 & a2] <- -1L
  both <- !a1 & !a2
  if (any(both)) {
    x <- d1[both]; y <- d2[both]
    w[both] <- ifelse(x > y, 1L, ifelse(x < y, -1L, 0L))
  }
  w
}

compare_relapse_vec <- function(r1, r2) {
  a1 <- is.na(r1); a2 <- is.na(r2)
  w  <- integer(length(r1))
  w[ a1 & !a2] <-  1L; w[!a1 & a2] <- -1L
  both <- !a1 & !a2
  if (any(both)) {
    x <- r1[both]; y <- r2[both]
    w[both] <- ifelse(x > y, 1L, ifelse(x < y, -1L, 0L))
  }
  w
}

compare_gvhd_vec <- function(g1, r1, o1, g2, r2, o2) {
  n <- length(g1); w <- integer(n)
  ok <- !is.na(g1) & !is.na(g2)
  if (any(ok)) {
    a <- g1[ok]; b <- g2[ok]
    w[ok] <- ifelse(a < b, 1L, ifelse(a > b, -1L, 0L))
  }
  gt <- which(ok & g1 == g2 & g1 > 0L)
  if (length(gt) > 0L) {
    r1_g <- r1[gt]; r2_g <- r2[gt]
    w_g  <- integer(length(gt))
    both_r <- !is.na(r1_g) & !is.na(r2_g)
    w_g[both_r & r1_g == 1L & r2_g == 0L] <-  1L
    w_g[both_r & r1_g == 0L & r2_g == 1L] <- -1L
    still <- which(both_r & w_g == 0L)
    if (length(still) > 0L) {
      o1_s <- o1[gt][still]; o2_s <- o2[gt][still]
      both_o <- !is.na(o1_s) & !is.na(o2_s)
      if (any(both_o)) {
        idx_o <- still[both_o]
        x <- o1[gt][idx_o]; y <- o2[gt][idx_o]
        w_g[idx_o] <- ifelse(x > y, 1L, ifelse(x < y, -1L, 0L))
      }
    }
    w[gt] <- w_g
  }
  w
}

compare_infection_vec <- function(i1, r1, i2, r2) {
  n <- length(i1); w <- integer(n)
  ok <- !is.na(i1) & !is.na(i2)
  if (any(ok)) {
    a <- i1[ok]; b <- i2[ok]
    w[ok] <- ifelse(a < b, 1L, ifelse(a > b, -1L, 0L))
  }
  gt <- which(ok & i1 == i2 & i1 > 0L)
  if (length(gt) > 0L) {
    r1_g <- r1[gt]; r2_g <- r2[gt]
    w_g  <- integer(length(gt))
    both_r <- !is.na(r1_g) & !is.na(r2_g)
    w_g[both_r & r1_g == 1L & r2_g == 0L] <-  1L
    w_g[both_r & r1_g == 0L & r2_g == 1L] <- -1L
    w[gt] <- w_g
  }
  w
}

build_pairwise_fast <- function(dat, pair_idx = NULL) {
  N <- nrow(dat)
  if (is.null(pair_idx)) {
    pair_idx <- t(combn(N, 2L))
    storage.mode(pair_idx) <- "integer"
    colnames(pair_idx) <- c("i", "j")
  }
  ii <- pair_idx[, 1]; jj <- pair_idx[, 2]

  w1 <- compare_survival_vec(dat$death_day[ii], dat$death_day[jj])
  w2 <- compare_relapse_vec(dat$relapse_day[ii], dat$relapse_day[jj])
  w3 <- compare_gvhd_vec(dat$gvhd_grade[ii], dat$gvhd_resolved[ii],
                         dat$gvhd_onset_day[ii],
                         dat$gvhd_grade[jj], dat$gvhd_resolved[jj],
                         dat$gvhd_onset_day[jj])
  w4 <- compare_infection_vec(dat$infection_grade[ii], dat$infection_resolved[ii],
                              dat$infection_grade[jj], dat$infection_resolved[jj])

  W <- as.integer(ifelse(w1 != 0L, w1,
                ifelse(w2 != 0L, w2,
                ifelse(w3 != 0L, w3,
                ifelse(w4 != 0L, w4, 0L)))))
  level <- as.integer(ifelse(w1 != 0L, 1L,
                     ifelse(w2 != 0L, 2L,
                     ifelse(w3 != 0L, 3L,
                     ifelse(w4 != 0L, 4L, 5L)))))
  data.frame(id_i = ii, id_j = jj, W = W, level = level)
}

