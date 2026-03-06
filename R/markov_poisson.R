.build_markov_oe <- function(paths, t_grid, transitions = NULL) {
  .validate_grid(t_grid, "t_grid")

  segments <- .flatten_paths(paths)
  if (nrow(segments) == 0) {
    stop("No path segments found.")
  }

  trans_df <- .normalize_transitions(transitions, segments)
  n_intervals <- length(t_grid) - 1

  interval_exposure <- matrix(0, nrow = nrow(segments), ncol = n_intervals)
  for (m in seq_len(n_intervals)) {
    interval_exposure[, m] <- pmax(
      0,
      pmin(segments$exit, t_grid[m + 1]) - pmax(segments$entry, t_grid[m])
    )
  }

  exposure_long <- data.frame(
    from = rep(segments$from, times = n_intervals),
    interval = rep(seq_len(n_intervals), each = nrow(segments)),
    exposure = as.vector(interval_exposure),
    stringsAsFactors = FALSE
  )
  exposure <- stats::aggregate(exposure ~ from + interval, data = exposure_long, FUN = sum)

  events <- segments[segments$from != segments$to, , drop = FALSE]
  if (nrow(events) > 0) {
    events$transition <- .format_transition(events$from, events$to)
    events$interval <- cut(events$exit, breaks = t_grid, labels = FALSE, right = TRUE)
    events <- events[!is.na(events$interval), c("transition", "interval"), drop = FALSE]
  }

  if (nrow(events) == 0) {
    counts <- data.frame(
      transition = character(0),
      interval = integer(0),
      count = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    counts <- stats::aggregate(rep(1, nrow(events)) ~ transition + interval, data = events, FUN = sum)
    names(counts)[3] <- "count"
  }

  design <- expand.grid(
    transition = trans_df$transition,
    interval = seq_len(n_intervals),
    stringsAsFactors = FALSE
  )
  design <- merge(design, trans_df, by = "transition", all.x = TRUE, sort = FALSE)

  oe <- merge(design, exposure, by = c("from", "interval"), all.x = TRUE, sort = FALSE)
  oe <- merge(oe, counts, by = c("transition", "interval"), all.x = TRUE, sort = FALSE)
  oe$exposure[is.na(oe$exposure)] <- 0
  oe$count[is.na(oe$count)] <- 0

  oe <- oe[order(oe$transition, oe$interval), ]
  rownames(oe) <- NULL
  oe
}

#' Fit Markov Poisson regressions on a time grid
#'
#' Fits one Poisson regression per transition using interval-specific effects and
#' offset `log(exposure)`. Intensities are piecewise constant over `t_grid`.
#'
#' @param paths List of jump paths. Each element must contain `times` and
#'   `states` vectors.
#' @param t_grid Strictly increasing time grid.
#' @param transitions Optional transitions to fit. Use `"i->j"`, shorthand
#'   `"ij"`, or a two-column matrix/data.frame with `(from, to)`.
#'
#' @return An object of class `jump_pois_markov_fit` containing fitted `glm`
#'   objects, aggregated occurrence-exposure data, and fitted interval rates.
#' @export
fit_markov_poisson <- function(paths, t_grid, transitions = NULL) {
  oe <- .build_markov_oe(paths, t_grid, transitions)
  trans_df <- unique(oe[, c("transition", "from", "to"), drop = FALSE])
  trans_df <- trans_df[order(trans_df$from, trans_df$to), , drop = FALSE]

  models <- stats::setNames(vector("list", nrow(trans_df)), trans_df$transition)
  model_info <- stats::setNames(vector("list", nrow(trans_df)), trans_df$transition)
  rate_list <- vector("list", nrow(trans_df))

  for (i in seq_len(nrow(trans_df))) {
    tr <- trans_df$transition[i]
    sub_all <- oe[oe$transition == tr, , drop = FALSE]
    sub_fit <- sub_all[sub_all$exposure > 0, , drop = FALSE]

    rate_out <- sub_all[, c("transition", "from", "to", "interval")]
    rate_out$rate <- NA_real_

    if (nrow(sub_fit) == 0) {
      model_info[[tr]] <- list(has_interval_factor = FALSE, interval_levels = character(0))
      rate_list[[i]] <- rate_out
      next
    }

    has_factor <- length(unique(sub_fit$interval)) > 1
    if (has_factor) {
      sub_fit$interval_factor <- factor(sub_fit$interval)
      model <- stats::glm(
        count ~ interval_factor + offset(log(exposure)),
        family = stats::poisson(),
        data = sub_fit
      )
      interval_levels <- levels(sub_fit$interval_factor)
      nd <- data.frame(
        interval_factor = factor(sub_all$interval, levels = interval_levels),
        exposure = rep(1, nrow(sub_all))
      )
      rate_out$rate <- as.numeric(stats::predict(model, newdata = nd, type = "response"))
    } else {
      model <- stats::glm(count ~ 1 + offset(log(exposure)), family = stats::poisson(), data = sub_fit)
      interval_levels <- as.character(unique(sub_fit$interval))
      rate_out$rate <- rep(as.numeric(stats::predict(model, newdata = data.frame(exposure = 1), type = "response"))[1], nrow(sub_all))
    }

    models[[tr]] <- model
    model_info[[tr]] <- list(has_interval_factor = has_factor, interval_levels = interval_levels)
    rate_list[[i]] <- rate_out
  }

  rates <- do.call(rbind, rate_list)
  rownames(rates) <- NULL

  out <- list(
    models = models,
    model_info = model_info,
    oe = oe,
    rates = rates,
    transitions = trans_df,
    t_grid = t_grid
  )
  class(out) <- "jump_pois_markov_fit"
  out
}

#' Predict Markov transition intensity
#'
#' Predicts fitted piecewise-constant intensity for one transition.
#'
#' @param object Fitted object from [fit_markov_poisson()].
#' @param time Numeric vector of times where prediction is needed.
#' @param transition Transition label (`"i->j"`). If omitted and only one
#'   transition was fitted, that transition is used.
#' @param type `"rate"` for intensity or `"count"` for expected count.
#' @param exposure Exposure used when `type = "count"`. Scalar or vector.
#'
#' @return Numeric vector of predicted values.
#' @export
predict_markov_poisson <- function(object, time, transition = NULL,
                                   type = c("rate", "count"), exposure = 1) {
  if (!inherits(object, "jump_pois_markov_fit")) {
    stop("`object` must come from `fit_markov_poisson()`.")
  }

  type <- match.arg(type)
  if (is.null(transition)) {
    if (length(object$models) != 1) {
      stop("Please provide `transition` when multiple transitions are fitted.")
    }
    transition <- names(object$models)[1]
  }
  if (!transition %in% names(object$models)) {
    stop("`transition` not found in fitted object.")
  }

  time <- as.numeric(time)
  interval <- cut(time, breaks = object$t_grid, labels = FALSE, right = TRUE)
  pred <- rep(NA_real_, length(time))

  model <- object$models[[transition]]
  info <- object$model_info[[transition]]
  if (is.null(model)) {
    return(pred)
  }

  valid <- !is.na(interval)
  if (any(valid)) {
    if (isTRUE(info$has_interval_factor)) {
      nd <- data.frame(
        interval_factor = factor(interval[valid], levels = info$interval_levels),
        exposure = rep(1, sum(valid))
      )
      pred[valid] <- as.numeric(stats::predict(model, newdata = nd, type = "response"))
    } else {
      pred[valid] <- as.numeric(stats::predict(
        model,
        newdata = data.frame(exposure = rep(1, sum(valid))),
        type = "response"
      ))
    }
  }

  if (type == "count") {
    if (length(exposure) == 1) {
      exposure <- rep(exposure, length(pred))
    }
    if (length(exposure) != length(pred)) {
      stop("`exposure` must be scalar or same length as `time`.")
    }
    pred <- pred * exposure
  }

  pred
}
