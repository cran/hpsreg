#' hpsreg: Regression for Hierarchically Ordered Clinical Outcomes
#' (Bradley-Terry-Davidson Model)
#'
#' A three-stage workflow:
#' \enumerate{
#'   \item **Design the hierarchy** with [hierarchy()] and the `lvl_*()`
#'     constructors ([lvl_binary()], [lvl_ordinal()], [lvl_continuous()],
#'     [lvl_tte()], [lvl_custom()]): binary, ordinal, continuous and
#'     time-to-event outcomes, tie margins, eligibility guards, groups.
#'   \item **Score pairs** with [score_pairs()] (compiled kernel, early exit at
#'     the first level that separates a pair), summarised by [resolution()].
#'   \item **Estimate**: covariate-adjusted Davidson model with sandwich standard
#'     errors, [btd()]; unadjusted [win_stats()]; checks with [btd_calibration()].
#' }
#' [simulate_hps()] and [hps_hierarchy()] provide a worked stem cell transplant
#' example (worked example: `system.file("examples", "hps-example.Rmd")`).
#'
#' @useDynLib hpsreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef vcov logLik nobs confint predict
#' @importFrom graphics par axis points lines segments abline legend mtext
#' @keywords internal
"_PACKAGE"

