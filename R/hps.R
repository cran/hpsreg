## ------------------------------------------------------------------
## Worked example: outcomes after allogeneic stem cell transplant
## (survival > relapse > GVHD > infection hierarchy, 180-day horizon)
## ------------------------------------------------------------------

#' Hierarchy for the stem cell transplant example
#'
#' Four levels, evaluated in this order: survival (later death better,
#' non-death best), relapse (later relapse better, non-relapse best), GVHD
#' (lower grade; then resolution; then later onset) and infection (lower
#' grade; then resolution). Resolution and onset are only compared between
#' patients who had the event.
#'
#' @return A [hierarchy()].
#' @examples
#' hps_hierarchy()
#' @export
hps_hierarchy <- function() {
  hierarchy(
    lvl_tte("death_day", better = "later", name = "survival"),
    lvl_tte("relapse_day", better = "later", name = "relapse"),
    lvl_group("GVHD",
      lvl_ordinal("gvhd_grade", better = "lower", name = "grade"),
      lvl_binary("gvhd_resolved", favorable = 1, when = "gvhd_grade > 0",
                 name = "resolution"),
      lvl_continuous("gvhd_onset_day", better = "higher", when = "gvhd_grade > 0",
                     name = "onset")),
    lvl_group("infection",
      lvl_ordinal("infection_grade", better = "lower", name = "grade"),
      lvl_binary("infection_resolved", favorable = 1, when = "infection_grade > 0",
                 name = "resolution"))
  )
}

#' Data-generating coefficients for the transplant example
#'
#' Starting values for [simulate_hps()]: the negative log hazard ratios of
#' published m-GRFS effects (uncalibrated), with fixed age and comorbidity
#' effects. Use these as a convenient default; effects should be re-calibrated
#' to a target if exact hazard ratios matter.
#'
#' @return A named numeric vector.
#' @examples
#' hps_beta()
#' @export
hps_beta <- function() {
  c(intercept = 0,
    MUD = -log(1.39), HID = -log(1.10),
    DRI_int = -log(1.45), DRI_high = -log(2.18), DRI_vhigh = -log(1.91),
    sex_mismatch = -log(1.52), year_2011 = -log(0.73),
    age_per_year = -0.010, hct_ci_high = -0.15)
}

#' Simulate a transplant cohort
#'
#' Simulates covariates and 180-day outcomes (death, relapse, graft-versus-host
#' disease, infection) for the worked example. Outcomes share a latent
#' frailty-like score driven by the covariates in `beta`.
#'
#' @param beta Named coefficient vector, see [hps_beta()].
#' @param N Number of patients.
#' @param seed Optional random seed.
#' @return A data frame with covariates (`donor_type`, `age`, `sex`, `dri`,
#'   `hct_ci`, `cell_source`, `conditioning`, `year`, `sex_mismatch`,
#'   `diagnosis`), the latent score, and outcome columns (`death_day`,
#'   `relapse_day`, `gvhd_grade`, `gvhd_resolved`, `gvhd_onset_day`,
#'   `infection_grade`, `infection_resolved`, `infection_onset_day`).
#'   `NA` in a `*_day` column means the event did not occur by day 180.
#' @examples
#' dat <- simulate_hps(N = 200, seed = 1)
#' head(dat[, c("donor_type", "age", "death_day", "gvhd_grade")])
#' @export
simulate_hps <- function(beta = hps_beta(), N = 613, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  req <- c("intercept", "MUD", "HID", "DRI_int", "DRI_high",
           "DRI_vhigh", "sex_mismatch", "year_2011",
           "age_per_year", "hct_ci_high")
  if (!all(req %in% names(beta)))
    stop("Missing names in `beta`: ", paste(setdiff(req, names(beta)), collapse = ", "),
         call. = FALSE)
  b <- function(nm) unname(beta[nm])

  donor_type   <- sample(c("MRD", "MUD", "HID"), N,
                         prob = c(212, 251, 150) / N, replace = TRUE)
  age          <- round(pmin(pmax(stats::rnorm(N, mean = 53, sd = 15), 18), 77))
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

  eta_true <- b("intercept") +
    b("MUD")          * (donor_type == "MUD") +
    b("HID")          * (donor_type == "HID") +
    b("DRI_int")      * (dri == "Intermediate") +
    b("DRI_high")     * (dri == "High") +
    b("DRI_vhigh")    * (dri == "Very high") +
    b("sex_mismatch") * (sex_mismatch == "F_to_M") +
    b("year_2011")    * (year == "2011-2015") +
    b("age_per_year") * (age - 53) +
    b("hct_ci_high")  * (hct_ci == ">=3")
  eta_true <- eta_true - mean(eta_true)

  sigma_eps <- 0.8
  eta       <- eta_true + stats::rnorm(N, 0, sigma_eps)
  eta_std   <- (eta - mean(eta)) / stats::sd(eta)

  cal <- function(lp, target) {
    stats::uniroot(function(a) mean(stats::plogis(a + lp)) - target,
                   interval = c(-20, 20))$root
  }

  lp_death <- -0.6 * eta_std +
    0.3 * (donor_type == "MUD") +
    0.6 * (donor_type == "HID") +
    0.3 * (hct_ci == ">=3")
  died_180 <- stats::rbinom(N, 1, stats::plogis(cal(lp_death, 0.12) + lp_death))
  death_day <- ifelse(died_180 == 1,
                      round(pmin(180, stats::rgamma(N, shape = 2, scale = 40))), NA)

  lp_relapse <- -0.5 * eta_std +
    0.4 * (dri == "Intermediate") +
    0.8 * (dri == "High") +
    0.7 * (dri == "Very high")
  relapse_180 <- ifelse(died_180 == 1, 0,
                        stats::rbinom(N, 1, stats::plogis(cal(lp_relapse, 0.15) + lp_relapse)))
  relapse_day <- ifelse(relapse_180 == 1, round(stats::runif(N, 30, 180)), NA)

  lp_gvhd <- -0.5 * eta_std +
     0.3 * (donor_type == "MUD") +
    -0.2 * (donor_type == "HID") +
     0.2 * (conditioning == "MAC")
  p_g1 <- stats::plogis(cal(lp_gvhd, 0.50) + lp_gvhd)
  p_g2 <- stats::plogis(cal(lp_gvhd, 0.20) + lp_gvhd)

  gvhd_grade <- rep(NA_integer_, N)
  alive      <- which(died_180 == 0)
  u          <- stats::runif(length(alive))
  gvhd_grade[alive] <- ifelse(u < p_g2[alive], 2L,
                       ifelse(u < p_g1[alive], 1L, 0L))
  gvhd_resolved <- ifelse(is.na(gvhd_grade) | gvhd_grade == 0, NA,
                   ifelse(gvhd_grade == 1, stats::rbinom(N, 1, 0.80),
                                           stats::rbinom(N, 1, 0.50)))
  gvhd_onset_day <- ifelse(!is.na(gvhd_grade) & gvhd_grade > 0,
                           round(stats::runif(N, 14, 120)), NA)

  lp_infection <- -0.4 * eta_std +
    0.2 * (hct_ci == ">=3") +
    0.2 * (conditioning == "MAC")
  p_i1 <- stats::plogis(cal(lp_infection, 0.45) + lp_infection)
  p_i2 <- stats::plogis(cal(lp_infection, 0.20) + lp_infection)

  infection_grade <- rep(NA_integer_, N)
  u <- stats::runif(length(alive))
  infection_grade[alive] <- ifelse(u < p_i2[alive], 2L,
                            ifelse(u < p_i1[alive], 1L, 0L))
  infection_resolved <- ifelse(is.na(infection_grade) | infection_grade == 0, NA,
                        ifelse(infection_grade == 1, stats::rbinom(N, 1, 0.85),
                                                      stats::rbinom(N, 1, 0.60)))
  infection_onset_day <- ifelse(!is.na(infection_grade) & infection_grade > 0,
                                round(stats::runif(N, 7, 180)), NA)

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
    eta_true            = as.numeric(eta_true),
    eta                 = as.numeric(eta),
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

#' Composite m-GRFS Cox model (comparator for the example)
#'
#' Fits a Cox model for modified GVHD- and relapse-free survival: the first of
#' death, relapse, grade >= 2 GVHD or grade >= 2 infection within the horizon.
#' This is the conventional composite analysis that the hierarchical analysis is
#' compared with in the example. Requires the \pkg{survival} package.
#'
#' @param dat Data from [simulate_hps()].
#' @param covars Character vector of covariate names.
#' @param horizon Follow-up horizon in days.
#' @return A list with `fit`, `coef` and `se`.
#' @examples
#' dat <- simulate_hps(N = 200, seed = 4)
#' hps_cox_grfs(dat, c("donor_type", "age"))$coef
#' @export
hps_cox_grfs <- function(dat, covars, horizon = 180) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required.", call. = FALSE)
  t_death     <- ifelse(is.na(dat$death_day), horizon, dat$death_day)
  t_relapse   <- ifelse(is.na(dat$relapse_day), horizon, dat$relapse_day)
  t_gvhd      <- ifelse(!is.na(dat$gvhd_grade) & dat$gvhd_grade >= 2,
                        ifelse(is.na(dat$gvhd_onset_day), horizon,
                               dat$gvhd_onset_day), horizon)
  t_infection <- ifelse(!is.na(dat$infection_grade) & dat$infection_grade >= 2,
                        ifelse(is.na(dat$infection_onset_day), horizon,
                               dat$infection_onset_day), horizon)
  t_grfs <- pmin(t_death, t_relapse, t_gvhd, t_infection, horizon)
  grfs_event <- as.integer(
    (!is.na(dat$death_day)    & dat$death_day    <= horizon) |
    (!is.na(dat$relapse_day)  & dat$relapse_day  <= horizon) |
    (!is.na(dat$gvhd_grade)   & dat$gvhd_grade >= 2 &
       !is.na(dat$gvhd_onset_day) & dat$gvhd_onset_day <= horizon) |
    (!is.na(dat$infection_grade) & dat$infection_grade >= 2 &
       !is.na(dat$infection_onset_day) & dat$infection_onset_day <= horizon))
  df <- cbind(dat[, covars, drop = FALSE], t_grfs = t_grfs, grfs_event = grfs_event)
  fml <- stats::as.formula(paste("survival::Surv(t_grfs, grfs_event) ~",
                                 paste(covars, collapse = " + ")))
  fit <- survival::coxph(fml, data = df)
  list(fit = fit, coef = stats::coef(fit), se = sqrt(diag(stats::vcov(fit))))
}
