## ------------------------------------------------------------------
## Pairwise scoring
## ------------------------------------------------------------------

#' Score all pairs of subjects under a hierarchy
#'
#' Applies a [hierarchy()] to every pair of subjects (or to a supplied set of
#' pairs). Built-in levels are evaluated by a compiled kernel that stops at the
#' first level that separates a pair; `lvl_custom()` levels are evaluated in R
#' on the pairs still unresolved.
#'
#' @param hierarchy A [hierarchy()].
#' @param data A data frame with one row per subject.
#' @param pairs `NULL` (all N(N-1)/2 pairs, the default) or a two-column integer
#'   matrix of row indices `(i, j)`, e.g. from [sample_pairs()]. See the note on
#'   variance in [btd()] for non-complete pair sets.
#' @param max_pairs Safety limit on the number of pairs held in memory.
#'
#' @return An object of class `hb_pairs` with components `id_i`, `id_j`
#'   (row indices), `W` (1 if `i` wins, -1 if `j` wins, 0 if tied on every
#'   level), `level` (index of the resolving top-level level, `K + 1` if
#'   unresolved), `leaf` (finer index when the hierarchy contains groups, else
#'   `NULL`), `uninf` (1 if an earlier level was uninformative because of
#'   censoring), `N`, `complete` and the `hierarchy`.
#' @examples
#' dat <- simulate_hps(N = 100, seed = 1)
#' pw <- score_pairs(hps_hierarchy(), dat)
#' pw
#' resolution(pw)
#' @export
score_pairs <- function(hierarchy, data, pairs = NULL, max_pairs = 5e7) {
  data <- as.data.frame(data)
  N <- nrow(data)
  if (N < 3L) stop("Need at least 3 subjects.", call. = FALSE)

  if (is.null(pairs)) {
    n <- N * (N - 1) / 2
    if (n > max_pairs)
      stop(sprintf(paste0("%.3g pairs exceed `max_pairs` (%.3g). Use ",
                          "sample_pairs() or raise `max_pairs`."), n, max_pairs),
           call. = FALSE)
    ids <- cpp_all_pairs(N)
    complete <- TRUE
  } else {
    pairs <- as.matrix(pairs)
    if (ncol(pairs) != 2L) stop("`pairs` must have two columns.", call. = FALSE)
    ids <- list(i = as.integer(pairs[, 1L]), j = as.integer(pairs[, 2L]))
    if (anyNA(ids$i) || anyNA(ids$j) || any(ids$i < 1L | ids$i > N) ||
        any(ids$j < 1L | ids$j > N) || any(ids$i == ids$j))
      stop("`pairs` must contain valid, distinct row indices.", call. = FALSE)
    complete <- FALSE
  }
  n <- length(ids$i)

  comp <- .compile_hierarchy(hierarchy, data)
  Kleaf <- comp$K
  W <- integer(n)
  leaf <- rep.int(Kleaf + 1L, n)
  uninf <- integer(n)
  active <- integer(0)
  all_active <- TRUE

  for (seg in comp$segments) {
    if (!all_active && length(active) == 0L) break
    if (seg$kind == "cpp") {
      r <- cpp_score_segment(seg$specs, comp$V, ids$i, ids$j,
                             if (all_active) integer(0) else active)
      res <- r$level > 0L
      if (all_active) {
        W <- r$W
        uninf <- r$uninf
        leaf[res] <- seg$leaf_idx[r$level[res]]
        active <- which(!res)
      } else {
        kk <- active[res]
        W[kk] <- r$W[res]
        leaf[kk] <- seg$leaf_idx[r$level[res]]
        uninf[active] <- pmax.int(uninf[active], r$uninf)
        active <- active[!res]
      }
    } else {
      idx <- if (all_active) seq_len(n) else active
      out <- seg$fun(data, ids$i[idx], ids$j[idx])
      if (length(out) != length(idx) || !all(out[!is.na(out)] %in% c(-1, 0, 1)))
        stop("lvl_custom() functions must return one value in {-1, 0, 1, NA} per pair.",
             call. = FALSE)
      if (seg$elig >= 0L) {
        e <- comp$V[, seg$elig + 1L]
        ok <- e[ids$i[idx]] == 1 & e[ids$j[idx]] == 1
        out[!ok] <- 0L
      }
      unf <- as.integer(is.na(out))
      out[is.na(out)] <- 0L
      res <- out != 0L
      if (all_active) {
        W <- as.integer(out); uninf <- unf
        leaf[res] <- seg$leaf_idx
        active <- which(!res)
      } else {
        kk <- idx[res]
        W[kk] <- as.integer(out[res])
        leaf[kk] <- seg$leaf_idx
        uninf[idx] <- pmax.int(uninf[idx], unf)
        active <- idx[!res]
      }
    }
    all_active <- FALSE
  }

  level <- if (hierarchy$grouped) hierarchy$leaf_group[pmin(leaf, Kleaf)] else leaf
  if (hierarchy$grouped) level[leaf == Kleaf + 1L] <- length(hierarchy) + 1L

  structure(list(id_i = ids$i, id_j = ids$j, W = W, level = level,
                 leaf = if (hierarchy$grouped) leaf else NULL,
                 uninf = uninf, N = N, complete = complete,
                 hierarchy = hierarchy),
            class = "hb_pairs")
}

#' Sample a random set of pairs
#'
#' Draws `n` distinct pairs uniformly from the N(N-1)/2 possible pairs without
#' enumerating them. Useful when N is too large for all pairs.
#'
#' @param N Number of subjects.
#' @param n Number of pairs to draw.
#' @return A two-column integer matrix `(i, j)` with `i < j`.
#' @examples
#' set.seed(1); head(sample_pairs(1000, 5))
#' @export
sample_pairs <- function(N, n) {
  total <- N * (N - 1) / 2
  if (n > total) stop("`n` exceeds the number of available pairs.", call. = FALSE)
  k <- sample.int(total, n)
  r <- cpp_index_to_pair(as.numeric(k), as.integer(N))
  cbind(i = r$i, j = r$j)
}

#' @export
print.hb_pairs <- function(x, ...) {
  K <- length(x$hierarchy)
  cat("<hb_pairs> ", format(length(x$W), big.mark = ","), " pairs among ",
      x$N, " subjects", if (!x$complete) " (subset of all pairs)", "\n", sep = "")
  tb <- tabulate(x$level, nbins = K + 1L)
  cat(sprintf("  wins/losses: %d resolved, %d tied on every level (%.1f%%)\n",
              sum(tb[seq_len(K)]), tb[K + 1L], 100 * tb[K + 1L] / length(x$W)))
  if (any(x$uninf == 1L))
    cat(sprintf("  %.1f%% of pairs passed an uninformative (censored) comparison\n",
                100 * mean(x$uninf == 1L)))
  invisible(x)
}

#' @export
length.hb_pairs <- function(x) length(x$W)

#' @export
as.data.frame.hb_pairs <- function(x, ...) {
  out <- data.frame(id_i = x$id_i, id_j = x$id_j, W = x$W, level = x$level)
  if (!is.null(x$leaf)) out$leaf <- x$leaf
  out$uninf <- x$uninf
  out
}

#' Where are pairs decided?
#'
#' Number and share of pairs resolved at each level of the hierarchy, and the
#' share left tied. Pass either scored pairs, or a hierarchy together with the
#' data: in that case the pairs are generated and counted on the fly and are
#' never stored, so this works for very large N.
#'
#' @param pairs An `hb_pairs` object from [score_pairs()], or a [hierarchy()].
#' @param by `"level"` (top-level levels, default) or `"leaf"` (each comparison
#'   inside groups).
#' @param data The data frame, when `pairs` is a hierarchy.
#' @return A data frame with columns `level`, `n` (pairs decided there) and
#'   `prop` (share of all pairs). The last row, `(tied)`, counts pairs that no
#'   level separates.
#' @examples
#' dat <- simulate_hps(N = 80, seed = 2)
#' resolution(score_pairs(hps_hierarchy(), dat))
#' resolution(hps_hierarchy(), data = dat)      # streaming, nothing stored
#' @export
resolution <- function(pairs, by = c("level", "leaf"), data = NULL) {
  by <- match.arg(by)
  if (inherits(pairs, "hb_hierarchy")) {
    h <- pairs
    if (is.null(data)) stop("`data` is required when `pairs` is a hierarchy.", call. = FALSE)
    comp <- .compile_hierarchy(h, as.data.frame(data))
    specs <- .stream_specs(comp)
    if (is.null(specs))
      stop("Hierarchies with lvl_custom() levels cannot be streamed; use score_pairs().",
           call. = FALSE)
    cn <- cpp_stream_counts(specs, comp$V, .NCHUNK, .n_threads())
    leaf_counts <- cn$counts
    Kleaf <- length(h$leaves)
    if (by == "leaf" && h$grouped) {
      lab <- c(h$leaf_names, "(tied)"); n <- leaf_counts
    } else {
      grp <- factor(h$leaf_group, levels = seq_along(h$level_names))
      n <- c(as.numeric(tapply(leaf_counts[seq_len(Kleaf)], grp, sum, default = 0)),
             leaf_counts[Kleaf + 1L])
      lab <- c(h$level_names, "(tied)")
    }
    return(data.frame(level = factor(lab, levels = lab), n = n,
                      prop = n / cn$n_pairs, row.names = NULL))
  }
  h <- pairs$hierarchy
  if (by == "leaf" && is.null(pairs$leaf)) by <- "level"
  if (by == "level") {
    lab <- c(h$level_names, "(tied)"); lv <- pairs$level
  } else {
    lab <- c(h$leaf_names, "(tied)"); lv <- pairs$leaf
  }
  K <- length(lab) - 1L
  n <- tabulate(lv, nbins = K + 1L)
  data.frame(level = factor(lab, levels = lab), n = n,
             prop = n / length(lv), row.names = NULL)
}
