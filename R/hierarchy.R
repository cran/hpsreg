## ------------------------------------------------------------------
## Hierarchy specification module
##
## A hierarchy is an ordered list of levels. For a pair of subjects the
## levels are examined top to bottom; the first level that separates them
## decides the pair (win / loss), otherwise the pair is a tie.
## ------------------------------------------------------------------

#' Outcome levels for a hierarchy
#'
#' Constructors for the building blocks of a [hierarchy()]. Each level compares
#' two subjects on one outcome and returns "first wins", "second wins" or
#' "neutral" (the pair is passed to the next level).
#'
#' @section Eligibility guards:
#' `when` is an expression evaluated in the data (e.g. `when = gvhd_grade > 0`).
#' Inside package code, where an unquoted column name would trigger an
#' "undefined global variable" note in `R CMD check`, pass it as a string
#' (`when = "gvhd_grade > 0"`); both forms are equivalent.
#' If either subject of a pair is ineligible (guard `FALSE` or `NA`) the level
#' is neutral for that pair. This expresses conditional comparisons such as
#' "compare resolution of GVHD only for patients who had GVHD".
#'
#' @section Missing values:
#' For `lvl_binary()`, `lvl_ordinal()` and `lvl_continuous()` the argument `na`
#' controls pairs in which exactly one subject has a missing value: `"tie"`
#' (default) makes the level neutral, `"best"` / `"worst"` treat a missing value
#' as the best / worst possible outcome. Pairs in which both are missing are
#' always neutral.
#'
#' @section Time-to-event levels:
#' With `status = NULL` the time variable is read as "time of the event, `NA` if
#' the event did not occur" and follow-up is assumed complete (no censoring);
#' non-events are placed at `horizon` (or infinity). With a `status` column
#' (1 = event, 0 = censored) pairs are compared with the Gehan rule: a pair is
#' decided only if the ordering is certain under right censoring, otherwise it
#' is *uninformative* and passed to the next level (and flagged, see
#' [score_pairs()]). Follow-up beyond `horizon` is truncated and treated as
#' complete at `horizon`. This rule is simple and biased toward the null under
#' heavy censoring; inverse-probability-of-censoring scoring is not yet
#' implemented.
#'
#' @param var Name of the outcome column in the data.
#' @param favorable Value of `var` that is the favourable outcome (binary
#'   levels). Logical, numeric, character and factor columns are supported.
#' @param better Direction that is favourable: `"higher"` or `"lower"` (ordinal,
#'   continuous); `"later"` or `"earlier"` (time-to-event). Must be supplied.
#' @param order Optional character vector giving the levels of an ordinal
#'   variable from lowest to highest. Defaults to the factor levels (or the
#'   sorted unique values).
#' @param threshold Tie margin (>= 0): differences up to this size are treated
#'   as ties (clinically meaningful difference).
#' @param time,status Names of the time and event-status columns
#'   (`lvl_tte()`); see Details.
#' @param horizon Administrative horizon for time-to-event levels.
#' @param when Optional eligibility guard, an expression, formula or string.
#' @param name Label used in printed output and resolution tables.
#' @param na How to treat a missing value in exactly one subject (see below).
#' @param fun For `lvl_custom()`: a function `function(data, i, j)` returning,
#'   for the pairs `(i[k], j[k])`, a vector in `{-1, 0, 1}` (1 = subject `i`
#'   wins) or `NA` for an uninformative pair.
#'
#' @return An object of class `hb_level`.
#' @examples
#' lvl_tte("death_day", better = "later", name = "survival")
#' lvl_ordinal("gvhd_grade", better = "lower")
#' lvl_binary("resolved", favorable = 1, when = gvhd_grade > 0)
#' lvl_continuous("pain_score", better = "lower", threshold = 1)
#' @name hb_levels
NULL

.capture_when <- function(expr, env) {
  if (is.null(expr)) return(NULL)
  if (is.character(expr)) expr <- parse(text = expr, keep.source = FALSE)[[1L]]
  if (is.call(expr) && identical(expr[[1L]], as.name("~")))
    expr <- expr[[length(expr)]]
  list(list(expr = expr, env = env))
}

.check_scalar_chr <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop("`", what, "` must be a single non-empty character string.", call. = FALSE)
}

.check_better <- function(better, choices) {
  if (missing(better))
    stop("`better` must be given explicitly (one of ",
         paste0("\"", choices, "\"", collapse = ", "),
         "): state which direction is favourable.", call. = FALSE)
  if (!is.character(better) || length(better) != 1L || !better %in% choices)
    stop("`better` must be one of ", paste0("\"", choices, "\"", collapse = ", "),
         ".", call. = FALSE)
  better
}

## NB: formals carry a leading dot so that extra arguments such as `na = ` can
## never be partially matched to them.
.new_level <- function(.type, .name, .var, .when, ...) {
  structure(c(list(type = .type, name = .name, var = .var, when = .when), list(...)),
            class = "hb_level")
}

.check_na <- function(na) match.arg(na, c("tie", "best", "worst"))

.check_threshold <- function(threshold) {
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0)
    stop("`threshold` must be a single number >= 0.", call. = FALSE)
  threshold
}

#' @rdname hb_levels
#' @export
lvl_binary <- function(var, favorable = 1, when = NULL, name = var, na = "tie") {
  .check_scalar_chr(var, "var")
  w <- .capture_when(substitute(when), parent.frame())
  force(name)
  .new_level("binary", name, var, w, favorable = favorable, na = .check_na(na),
             threshold = 0)
}

#' @rdname hb_levels
#' @export
lvl_ordinal <- function(var, better, order = NULL, threshold = 0,
                        when = NULL, name = var, na = "tie") {
  .check_scalar_chr(var, "var")
  better <- .check_better(better, c("higher", "lower"))
  w <- .capture_when(substitute(when), parent.frame())
  force(name)
  .new_level("ordinal", name, var, w, better = better, order = order,
             na = .check_na(na), threshold = .check_threshold(threshold))
}

#' @rdname hb_levels
#' @export
lvl_continuous <- function(var, better, threshold = 0, when = NULL,
                           name = var, na = "tie") {
  .check_scalar_chr(var, "var")
  better <- .check_better(better, c("higher", "lower"))
  w <- .capture_when(substitute(when), parent.frame())
  force(name)
  .new_level("continuous", name, var, w, better = better,
             na = .check_na(na), threshold = .check_threshold(threshold))
}

#' @rdname hb_levels
#' @export
lvl_tte <- function(time, status = NULL, better, horizon = Inf, threshold = 0,
                    when = NULL, name = time) {
  .check_scalar_chr(time, "time")
  if (!is.null(status)) .check_scalar_chr(status, "status")
  better <- .check_better(better, c("later", "earlier"))
  if (!is.numeric(horizon) || length(horizon) != 1L || is.na(horizon) || horizon <= 0)
    stop("`horizon` must be a single positive number (or Inf).", call. = FALSE)
  w <- .capture_when(substitute(when), parent.frame())
  force(name)
  .new_level("tte", name, time, w, status = status, better = better,
             horizon = horizon, threshold = .check_threshold(threshold),
             na = "tie")
}

#' @rdname hb_levels
#' @export
lvl_custom <- function(fun, name = "custom", when = NULL) {
  if (!is.function(fun)) stop("`fun` must be a function(data, i, j).", call. = FALSE)
  w <- .capture_when(substitute(when), parent.frame())
  .new_level("custom", name, NA_character_, w, fun = fun)
}

#' Group several levels into one reporting level
#'
#' A group is evaluated as consecutive levels (in the order given) but is
#' reported as a single level, e.g. GVHD grade, then resolution, then onset
#' time are three comparisons forming one "GVHD" level of the hierarchy.
#'
#' @param name Label of the group.
#' @param ... Levels created by the `lvl_*()` constructors.
#' @param when Optional eligibility guard applied to every level in the group.
#' @return An object of class `hb_group`.
#' @export
lvl_group <- function(name, ..., when = NULL) {
  .check_scalar_chr(name, "name")
  kids <- list(...)
  if (!length(kids)) stop("A group needs at least one level.", call. = FALSE)
  if (!all(vapply(kids, inherits, NA, what = "hb_level")))
    stop("Groups may only contain levels created by lvl_*().", call. = FALSE)
  w <- .capture_when(substitute(when), parent.frame())
  structure(list(name = name, levels = kids, when = w), class = "hb_group")
}

#' Define a hierarchy of outcomes
#'
#' Levels are listed from highest to lowest priority. When two subjects are
#' compared, the first level on which they differ decides who wins; if no level
#' separates them the pair is a tie.
#'
#' @param ... Levels (`lvl_*()`) and/or groups ([lvl_group()]).
#' @return An object of class `hb_hierarchy`.
#' @examples
#' h <- hierarchy(
#'   lvl_tte("death_day", better = "later", name = "survival"),
#'   lvl_group("GVHD",
#'     lvl_ordinal("gvhd_grade", better = "lower"),
#'     lvl_binary("gvhd_resolved", favorable = 1, when = gvhd_grade > 0)
#'   )
#' )
#' h
#' @export
hierarchy <- function(...) {
  items <- list(...)
  if (length(items) == 1L && is.list(items[[1L]]) &&
      !inherits(items[[1L]], c("hb_level", "hb_group")))
    items <- items[[1L]]
  if (!length(items)) stop("A hierarchy needs at least one level.", call. = FALSE)
  ok <- vapply(items, inherits, NA, what = c("hb_level", "hb_group"))
  # inherits(what = vector) is TRUE if any class matches
  if (!all(ok)) stop("All elements must be levels (lvl_*) or groups (lvl_group).",
                     call. = FALSE)

  top_names <- vapply(items, function(x) x$name, "")
  if (anyNA(top_names) || any(!nzchar(top_names)))
    stop("Every level needs a non-empty name.", call. = FALSE)
  if (anyDuplicated(top_names))
    stop("Level names must be unique: ",
         paste(unique(top_names[duplicated(top_names)]), collapse = ", "),
         call. = FALSE)

  leaves <- list(); leaf_group <- integer(); leaf_names <- character()
  for (k in seq_along(items)) {
    it <- items[[k]]
    if (inherits(it, "hb_group")) {
      for (child in it$levels) {
        lf <- child
        lf$when <- c(it$when, lf$when)
        leaves[[length(leaves) + 1L]] <- lf
        leaf_group <- c(leaf_group, k)
        leaf_names <- c(leaf_names, paste0(it$name, ": ", child$name))
      }
    } else {
      leaves[[length(leaves) + 1L]] <- it
      leaf_group <- c(leaf_group, k)
      leaf_names <- c(leaf_names, it$name)
    }
  }
  structure(list(items = items, leaves = leaves,
                 level_names = top_names,
                 leaf_names = make.unique(leaf_names),
                 leaf_group = leaf_group,
                 grouped = length(leaves) != length(items)),
            class = "hb_hierarchy")
}

#' @export
length.hb_hierarchy <- function(x) length(x$level_names)

.describe_leaf <- function(lf) {
  dir <- if (!is.null(lf$better)) paste0(lf$better, " better") else
    if (lf$type == "binary") paste0("favourable = ", format(lf$favorable)) else ""
  extra <- character()
  if (!is.null(lf$threshold) && lf$threshold > 0)
    extra <- c(extra, paste0("margin ", format(lf$threshold)))
  if (lf$type == "tte") {
    extra <- c(extra, if (is.null(lf$status)) "complete follow-up" else
      paste0("status = ", lf$status, " (Gehan)"))
    if (is.finite(lf$horizon)) extra <- c(extra, paste0("horizon ", format(lf$horizon)))
  }
  if (!is.null(lf$na) && lf$na != "tie") extra <- c(extra, paste0("NA = ", lf$na))
  if (length(lf$when))
    extra <- c(extra, paste0("when ", paste(vapply(lf$when, function(g)
      paste(deparse(g$expr), collapse = ""), ""), collapse = " & ")))
  paste0(lf$type, if (!is.na(lf$var)) paste0(" [", lf$var, "]") else "",
         if (nzchar(dir)) paste0(", ", dir) else "",
         if (length(extra)) paste0(" (", paste(extra, collapse = "; "), ")") else "")
}

#' @export
print.hb_hierarchy <- function(x, ...) {
  cat("<hierarchy> ", length(x), " level", if (length(x) != 1L) "s", " ",
      "(first level that separates a pair decides it)\n", sep = "")
  for (k in seq_along(x$items)) {
    it <- x$items[[k]]
    if (inherits(it, "hb_group")) {
      cat(sprintf("%2d. %s\n", k, it$name))
      g <- if (length(it$when)) paste0("    when ", paste(vapply(it$when, function(g)
        paste(deparse(g$expr), collapse = ""), ""), collapse = " & "), "\n") else ""
      cat(g, sep = "")
      for (ch in it$levels)
        cat(sprintf("      - %s: %s\n", ch$name, .describe_leaf(ch)))
    } else {
      cat(sprintf("%2d. %s: %s\n", k, it$name, .describe_leaf(it)))
    }
  }
  invisible(x)
}

#' @export
print.hb_level <- function(x, ...) {
  cat("<level> ", x$name, ": ", .describe_leaf(x), "\n", sep = "")
  invisible(x)
}

#' @export
print.hb_group <- function(x, ...) {
  cat("<level group> ", x$name, "\n", sep = "")
  for (ch in x$levels) cat("  - ", ch$name, ": ", .describe_leaf(ch), "\n", sep = "")
  invisible(x)
}

## ------------------------------------------------------------------
## Compilation: turn a hierarchy + data into a numeric matrix and
## kernel specifications for the C++ engine.
## ------------------------------------------------------------------

.get_col <- function(data, var, what) {
  if (!var %in% names(data))
    stop("Column '", var, "' (", what, ") not found in `data`.", call. = FALSE)
  data[[var]]
}

.numeric_values <- function(lf, data) {
  x <- .get_col(data, lf$var, lf$name)
  switch(lf$type,
    binary = {
      as.numeric(x == lf$favorable)
    },
    ordinal = {
      if (is.factor(x) || is.character(x)) {
        lv <- if (!is.null(lf$order)) lf$order else
          if (is.factor(x)) levels(x) else sort(unique(x))
        f <- factor(x, levels = lv)
        if (any(is.na(f) & !is.na(x)))
          stop("Level '", lf$name, "': values of '", lf$var,
               "' not covered by `order`.", call. = FALSE)
        as.numeric(f)
      } else as.numeric(x)
    },
    continuous = as.numeric(x),
    stop("internal: unexpected level type"))
}

.tte_values <- function(lf, data) {
  N <- nrow(data)
  tt  <- as.numeric(.get_col(data, lf$var, lf$name))
  tau <- lf$horizon
  if (is.null(lf$status)) {
    e <- rep(1, N)
    t <- tt
    t[is.na(t)] <- tau
    t[t > tau]  <- tau
  } else {
    st <- .get_col(data, lf$status, paste0(lf$name, " status"))
    if (is.logical(st)) st <- as.numeric(st)
    if (!all(st[!is.na(st)] %in% c(0, 1)))
      stop("Status column '", lf$status, "' must be coded 0/1 (or logical).",
           call. = FALSE)
    e <- as.numeric(st)
    t <- tt
    beyond <- !is.na(t) & t >= tau
    t[beyond] <- tau
    e[beyond] <- 1
    t[is.na(e)] <- NA_real_
  }
  list(t = t, e = e)
}

.eval_guards <- function(guards, data) {
  N <- nrow(data)
  ok <- rep(TRUE, N)
  for (g in guards) {
    v <- eval(g$expr, data, g$env)
    if (!is.logical(v) || length(v) != N)
      stop("Guard `", paste(deparse(g$expr), collapse = ""),
           "` must evaluate to a logical vector of length nrow(data).", call. = FALSE)
    ok <- ok & !is.na(v) & v
  }
  as.numeric(ok)
}

.compile_hierarchy <- function(h, data) {
  if (!inherits(h, "hb_hierarchy")) stop("`hierarchy` must come from hierarchy().", call. = FALSE)
  data <- as.data.frame(data)
  N <- nrow(data)
  cols <- list(); add_col <- function(v) {
    cols[[length(cols) + 1L]] <<- as.numeric(v)
    length(cols) - 1L          # 0-based index
  }

  K <- length(h$leaves)
  specs <- vector("list", K)
  for (k in seq_len(K)) {
    lf <- h$leaves[[k]]
    elig <- if (length(lf$when)) add_col(.eval_guards(lf$when, data)) else -1L
    if (lf$type == "custom") {
      specs[[k]] <- list(kind = "custom", fun = lf$fun, elig = elig)
      next
    }
    if (lf$type == "tte") {
      tv <- .tte_values(lf, data)
      specs[[k]] <- list(kind = "cpp", type = 2L, col = add_col(tv$t),
                         col2 = add_col(tv$e), elig = elig, na_rule = 0L,
                         sign = if (lf$better == "later") 1 else -1,
                         thr = lf$threshold)
    } else {
      v <- .numeric_values(lf, data)
      better <- if (lf$type == "binary") "higher" else lf$better
      specs[[k]] <- list(kind = "cpp", type = 1L, col = add_col(v),
                         col2 = -1L, elig = elig,
                         na_rule = match(lf$na, c("tie", "best", "worst")) - 1L,
                         sign = if (better == "higher") 1 else -1,
                         thr = lf$threshold)
    }
  }
  V <- if (length(cols)) do.call(cbind, cols) else matrix(0, N, 1L)
  storage.mode(V) <- "double"

  ## consecutive built-in levels form one segment evaluated in a single C++ pass
  segments <- list()
  k <- 1L
  while (k <= K) {
    if (specs[[k]]$kind == "custom") {
      segments[[length(segments) + 1L]] <- list(kind = "custom", leaf_idx = k,
                                                fun = specs[[k]]$fun,
                                                elig = specs[[k]]$elig)
      k <- k + 1L
    } else {
      idx <- k
      while (idx[length(idx)] < K && specs[[idx[length(idx)] + 1L]]$kind == "cpp")
        idx <- c(idx, idx[length(idx)] + 1L)
      segments[[length(segments) + 1L]] <- list(
        kind = "cpp", leaf_idx = idx,
        specs = lapply(specs[idx], function(s) s[c("type", "col", "col2", "elig",
                                                   "na_rule", "sign", "thr")]))
      k <- idx[length(idx)] + 1L
    }
  }
  list(V = V, segments = segments, N = N, K = K, elig_cols = V)
}

## Streaming needs every level to be a built-in kernel level (no R callbacks).
.stream_specs <- function(comp) {
  if (!all(vapply(comp$segments, function(s) identical(s$kind, "cpp"), NA))) return(NULL)
  unlist(lapply(comp$segments, `[[`, "specs"), recursive = FALSE)
}

.NCHUNK <- 32L

.n_threads <- function(control = NULL) {
  t <- control$threads
  if (is.null(t)) t <- getOption("hpsreg.threads", 1L)
  as.integer(max(1L, t))
}

#' OpenMP availability
#'
#' Number of threads OpenMP can use in this build of the package, or `0` if the
#' package was compiled without OpenMP. The streaming engine uses one thread by
#' default; set `options(hpsreg.threads = n)` (or `btd_control(threads = n)`) to
#' use more. Results are bitwise identical for any number of threads.
#'
#' @return An integer.
#' @examples
#' openmp_threads()
#' @export
openmp_threads <- function() cpp_openmp_threads()
