## ------------------------------------------------------------------
## Unadjusted win statistics for a two-group comparison
## ------------------------------------------------------------------

#' Win statistics for two groups
#'
#' Summarises the cross-group pairs of a scored hierarchy: the proportion of
#' pairs won by each group, and the derived net benefit, win ratio and win odds
#' of `treatment` versus the reference group, with U-statistic (Hoeffding)
#' variances that account for each subject appearing in many pairs.
#'
#' With \eqn{p_w} and \eqn{p_l} the proportions of cross-group pairs won by the
#' treatment and reference groups and \eqn{p_t = 1 - p_w - p_l}:
#' net benefit \eqn{p_w - p_l}; win ratio \eqn{p_w/p_l};
#' win odds \eqn{(p_w + p_t/2)/(p_l + p_t/2)}. Confidence intervals for the two
#' ratios are computed on the log scale (delta method).
#'
#' Requires the complete set of pairs. Pass a hierarchy together with `data`
#' instead of scored pairs to compute everything on the fly without storing
#' pairs (memory O(N)); this is the way to analyse very large N.
#'
#' @param pairs An `hb_pairs` object from [score_pairs()] (complete pair set),
#'   or a [hierarchy()] (then `data` is required).
#' @param group A vector with exactly two distinct values, one per subject.
#' @param treatment Value of `group` defining the treatment group; defaults to
#'   the second level (for a factor) or the larger value.
#' @param conf.level Confidence level.
#' @param data The data frame, when `pairs` is a hierarchy.
#' @return An object of class `hb_win` with the summary `table`, the
#'   `by_level` contribution of each hierarchy level to wins and losses, and the
#'   counts.
#' @examples
#' dat <- simulate_hps(N = 150, seed = 7)
#' pw <- score_pairs(hps_hierarchy(), dat)
#' win_stats(pw, dat$year)
#' win_stats(hps_hierarchy(), dat$year, data = dat)   # streaming, same result
#' @export
win_stats <- function(pairs, group, treatment = NULL, conf.level = 0.95, data = NULL) {
  streaming <- inherits(pairs, "hb_hierarchy")
  if (!streaming && !inherits(pairs, "hb_pairs"))
    stop("`pairs` must come from score_pairs() or be a hierarchy().", call. = FALSE)
  if (streaming) {
    if (is.null(data)) stop("`data` is required when `pairs` is a hierarchy.", call. = FALSE)
    h <- pairs; N <- nrow(data)
  } else {
    if (!pairs$complete)
      stop("win_stats() requires the complete set of pairs.", call. = FALSE)
    h <- pairs$hierarchy; N <- pairs$N
  }
  if (length(group) != N) stop("`group` must have one value per subject.", call. = FALSE)
  if (anyNA(group)) stop("`group` must not contain missing values.", call. = FALSE)
  lv <- if (is.factor(group)) levels(droplevels(group)) else sort(unique(group))
  if (length(lv) != 2L) stop("`group` must have exactly two distinct values.", call. = FALSE)
  trt <- if (is.null(treatment)) lv[2L] else treatment
  if (!as.character(trt) %in% as.character(lv))
    stop("`treatment` is not a value of `group`.", call. = FALSE)
  ref <- setdiff(as.character(lv), as.character(trt))
  g <- as.integer(as.character(group) == as.character(trt))
  n1 <- sum(g == 1L); n0 <- sum(g == 0L)
  if (n1 < 2L || n0 < 2L) stop("Each group needs at least 2 subjects.", call. = FALSE)
  ncross <- as.numeric(n1) * n0
  K <- length(h)
  lab <- c(h$level_names, "(tied)")

  if (streaming) {
    comp <- .compile_hierarchy(h, as.data.frame(data))
    specs <- .stream_specs(comp)
    if (is.null(specs))
      stop("Hierarchies with lvl_custom() levels cannot be streamed; use score_pairs().",
           call. = FALSE)
    cr <- cpp_stream_cross(specs, comp$V, g, .NCHUNK, .n_threads())
    cs <- cr$cross_sums
    Kleaf <- length(h$leaves)
    grp <- factor(h$leaf_group, levels = seq_len(K))
    agg <- function(col) c(as.numeric(tapply(cr$by_level[seq_len(Kleaf), col], grp, sum,
                                             default = 0)), 0)
    tw <- agg(1L) / ncross; tl <- agg(2L) / ncross
    tied <- c(rep(0, K), cr$by_level[Kleaf + 1L, 3L] / ncross)
  } else {
    cs <- cpp_cross_sums(pairs$id_i, pairs$id_j, pairs$W, g)   # wins, losses, ties
    cross <- g[pairs$id_i] != g[pairs$id_j]
    trt_wins <- ifelse(g[pairs$id_i] == 1L, pairs$W == 1L, pairs$W == -1L)
    trt_loss <- ifelse(g[pairs$id_i] == 1L, pairs$W == -1L, pairs$W == 1L)
    tw <- tabulate(pairs$level[cross & trt_wins], K + 1L) / ncross
    tl <- tabulate(pairs$level[cross & trt_loss], K + 1L) / ncross
    tied <- c(rep(0, K), sum(cross & pairs$W == 0L) / ncross)
  }

  trt_idx <- which(g == 1L); ref_idx <- which(g == 0L)
  # per-subject proportions of the n1*n0 cross pairs won / lost by the treatment group
  Wa <- cs[trt_idx, 1L] / n0; La <- cs[trt_idx, 2L] / n0
  Wb <- cs[ref_idx, 2L] / n1; Lb <- cs[ref_idx, 1L] / n1
  pw <- mean(Wa); pl <- mean(La); pt <- 1 - pw - pl
  vw <- stats::var(Wa) / n1 + stats::var(Wb) / n0
  vl <- stats::var(La) / n1 + stats::var(Lb) / n0
  cwl <- stats::cov(Wa, La) / n1 + stats::cov(Wb, Lb) / n0

  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  nb <- pw - pl
  vnb <- vw + vl - 2 * cwl
  lwr <- log(pw / pl)
  vlwr <- vw / pw^2 + vl / pl^2 - 2 * cwl / (pw * pl)
  lwo <- log((1 + nb) / (1 - nb))
  vlwo <- vnb * (2 / (1 - nb^2))^2

  est <- c(net_benefit = nb, win_ratio = pw / pl, win_odds = exp(lwo))
  se  <- c(sqrt(vnb), sqrt(vlwr), sqrt(vlwo))     # SE on the (log) scale of the CI
  low <- c(nb - z * se[1L], exp(lwr - z * se[2L]), exp(lwo - z * se[3L]))
  upp <- c(nb + z * se[1L], exp(lwr + z * se[2L]), exp(lwo + z * se[3L]))
  zst <- c(nb / se[1L], lwr / se[2L], lwo / se[3L])
  tab <- data.frame(statistic = names(est), estimate = unname(est),
                    lower = low, upper = upp, se_scale = se,
                    z = zst, p = 2 * stats::pnorm(-abs(zst)), row.names = NULL)

  by_level <- data.frame(level = factor(lab, levels = lab),
                         treatment_wins = tw, treatment_losses = tl, tied = tied,
                         row.names = NULL)

  structure(list(table = tab, p_win = pw, p_loss = pl, p_tie = pt,
                 n_treatment = n1, n_reference = n0,
                 treatment = as.character(trt), reference = ref,
                 by_level = by_level, conf.level = conf.level),
            class = "hb_win")
}

#' @export
print.hb_win <- function(x, digits = 4, ...) {
  cat(sprintf("Win statistics: '%s' (n = %d) vs '%s' (n = %d)\n", x$treatment,
              x$n_treatment, x$reference, x$n_reference))
  cat(sprintf("Cross-group pairs: win %.3f | loss %.3f | tie %.3f\n\n",
              x$p_win, x$p_loss, x$p_tie))
  tb <- x$table
  num <- c("estimate", "lower", "upper", "z", "p")
  tb[num] <- lapply(tb[num], signif, digits)
  print(tb[c("statistic", num)], row.names = FALSE)
  cat(sprintf("\n%g%% CI; ratios use log-scale intervals. Contribution by level:\n",
              100 * x$conf.level))
  bl <- x$by_level
  bl[-1L] <- lapply(bl[-1L], round, digits)
  print(bl, row.names = FALSE)
  invisible(x)
}
