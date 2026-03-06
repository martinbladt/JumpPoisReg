.segment_exposure_2d <- function(entry, exit, t_grid, u_grid) {
  if (!is.finite(entry) || !is.finite(exit) || exit <= entry) {
    return(numeric(0))
  }

  t_cuts <- t_grid[t_grid > entry & t_grid < exit]
  u_cuts <- entry + u_grid[u_grid > 0]
  u_cuts <- u_cuts[u_cuts > entry & u_cuts < exit]

  cuts <- sort(unique(c(entry, exit, t_cuts, u_cuts)))
  if (length(cuts) < 2) {
    return(numeric(0))
  }

  expo <- numeric(0)
  for (k in seq_len(length(cuts) - 1)) {
    left <- cuts[k]
    right <- cuts[k + 1]
    if (right <= left) {
      next
    }

    mid_t <- (left + right) / 2
    mid_u <- mid_t - entry

    t_box <- cut(mid_t, breaks = t_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)
    u_box <- cut(mid_u, breaks = u_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)

    if (is.na(t_box) || is.na(u_box)) {
      next
    }

    box <- paste0(t_box, "_", u_box)
    if (box %in% names(expo)) {
      expo[box] <- expo[box] + (right - left)
    } else {
      expo[box] <- right - left
    }
  }

  expo
}

.build_semi_markov_oe <- function(paths, t_grid, u_grid, transitions = NULL) {
  .validate_grid(t_grid, "t_grid")
  .validate_grid(u_grid, "u_grid")

  segments <- .flatten_paths(paths)
  if (nrow(segments) == 0) {
    stop("No path segments found.")
  }

  trans_df <- .normalize_transitions(transitions, segments)

  exposure_rows <- vector("list", nrow(segments))
  for (i in seq_len(nrow(segments))) {
    ex <- .segment_exposure_2d(segments$entry[i], segments$exit[i], t_grid, u_grid)
    if (length(ex) == 0) {
      next
    }

    exposure_rows[[i]] <- data.frame(
      from = segments$from[i],
      box = names(ex),
      exposure = as.numeric(ex),
      stringsAsFactors = FALSE
    )
  }

  exposure_rows <- Filter(Negate(is.null), exposure_rows)
  if (length(exposure_rows) == 0) {
    stop("No exposure found inside the supplied (t_grid, u_grid).")
  }

  exposure <- do.call(rbind, exposure_rows)
  exposure <- stats::aggregate(exposure ~ from + box, data = exposure, FUN = sum)

  events <- segments[segments$from != segments$to, , drop = FALSE]
  if (nrow(events) > 0) {
    events$transition <- .format_transition(events$from, events$to)
    events$t_box <- cut(events$exit, breaks = t_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)
    events$u_box <- cut(events$duration, breaks = u_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)
    events <- events[!is.na(events$t_box) & !is.na(events$u_box), , drop = FALSE]
    events$box <- paste0(events$t_box, "_", events$u_box)
    events <- events[, c("transition", "box"), drop = FALSE]
  }

  if (nrow(events) == 0) {
    counts <- data.frame(
      transition = character(0),
      box = character(0),
      count = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    counts <- stats::aggregate(rep(1, nrow(events)) ~ transition + box, data = events, FUN = sum)
    names(counts)[3] <- "count"
  }

  boxes <- sort(unique(exposure$box))
  design <- expand.grid(transition = trans_df$transition, box = boxes, stringsAsFactors = FALSE)
  design <- merge(design, trans_df, by = "transition", all.x = TRUE, sort = FALSE)

  oe <- merge(design, exposure, by = c("from", "box"), all.x = TRUE, sort = FALSE)
  oe <- merge(oe, counts, by = c("transition", "box"), all.x = TRUE, sort = FALSE)
  oe$exposure[is.na(oe$exposure)] <- 0
  oe$count[is.na(oe$count)] <- 0

  box_index <- .box_to_indices(oe$box)
  oe$t_box <- box_index$t_box
  oe$u_box <- box_index$u_box

  oe <- oe[order(oe$transition, oe$t_box, oe$u_box), ]
  rownames(oe) <- NULL
  oe
}

#' Fit semi-Markov Poisson regressions on a time-duration grid
#'
#' Fits one Poisson regression per transition with grid-box-specific effects and
#' offset `log(exposure)`. Intensities are piecewise constant on the
#' `(time, duration)` boxes induced by `t_grid` and `u_grid`.
#'
#' @param paths List of jump paths. Each element must contain `times` and
#'   `states` vectors.
#' @param t_grid Strictly increasing time grid.
#' @param u_grid Strictly increasing duration grid.
#' @param transitions Optional transitions to fit. Use `"i->j"`, shorthand
#'   `"ij"`, or a two-column matrix/data.frame with `(from, to)`.
#'
#' @return An object of class `jump_pois_semi_markov_fit` with fitted `glm`
#'   objects, occurrence-exposure table, and fitted box rates.
#' @export
fit_semi_markov_poisson <- function(paths, t_grid, u_grid, transitions = NULL) {
  oe <- .build_semi_markov_oe(paths, t_grid, u_grid, transitions)
  trans_df <- unique(oe[, c("transition", "from", "to"), drop = FALSE])
  trans_df <- trans_df[order(trans_df$from, trans_df$to), , drop = FALSE]

  models <- stats::setNames(vector("list", nrow(trans_df)), trans_df$transition)
  model_info <- stats::setNames(vector("list", nrow(trans_df)), trans_df$transition)
  rate_list <- vector("list", nrow(trans_df))

  for (i in seq_len(nrow(trans_df))) {
    tr <- trans_df$transition[i]
    sub_all <- oe[oe$transition == tr, , drop = FALSE]
    sub_fit <- sub_all[sub_all$exposure > 0, , drop = FALSE]

    rate_out <- sub_all[, c("transition", "from", "to", "box", "t_box", "u_box")]
    rate_out$rate <- NA_real_

    if (nrow(sub_fit) == 0) {
      model_info[[tr]] <- list(has_box_factor = FALSE, box_levels = character(0))
      rate_list[[i]] <- rate_out
      next
    }

    has_factor <- length(unique(sub_fit$box)) > 1
    if (has_factor) {
      sub_fit$box_factor <- factor(sub_fit$box)
      model <- stats::glm(
        count ~ box_factor + offset(log(exposure)),
        family = stats::poisson(),
        data = sub_fit
      )
      box_levels <- levels(sub_fit$box_factor)
      nd <- data.frame(
        box_factor = factor(sub_all$box, levels = box_levels),
        exposure = rep(1, nrow(sub_all))
      )
      rate_out$rate <- as.numeric(stats::predict(model, newdata = nd, type = "response"))
    } else {
      model <- stats::glm(count ~ 1 + offset(log(exposure)), family = stats::poisson(), data = sub_fit)
      box_levels <- as.character(unique(sub_fit$box))
      rate_out$rate <- rep(as.numeric(stats::predict(model, newdata = data.frame(exposure = 1), type = "response"))[1], nrow(sub_all))
    }

    models[[tr]] <- model
    model_info[[tr]] <- list(has_box_factor = has_factor, box_levels = box_levels)
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
    t_grid = t_grid,
    u_grid = u_grid
  )
  class(out) <- "jump_pois_semi_markov_fit"
  out
}

#' Predict semi-Markov transition intensity
#'
#' Predicts fitted piecewise-constant intensity for one transition at
#' `(time, duration)` points.
#'
#' @param object Fitted object from [fit_semi_markov_poisson()].
#' @param time Numeric vector of times.
#' @param duration Numeric vector of durations since the last jump.
#' @param transition Transition label (`"i->j"`). If omitted and only one
#'   transition was fitted, that transition is used.
#' @param type `"rate"` for intensity or `"count"` for expected count.
#' @param exposure Exposure used when `type = "count"`. Scalar or vector.
#'
#' @return Numeric vector of predicted values.
#' @export
predict_semi_markov_poisson <- function(object, time, duration, transition = NULL,
                                        type = c("rate", "count"), exposure = 1) {
  if (!inherits(object, "jump_pois_semi_markov_fit")) {
    stop("`object` must come from `fit_semi_markov_poisson()`.")
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
  duration <- as.numeric(duration)
  if (length(time) != length(duration)) {
    stop("`time` and `duration` must have equal length.")
  }

  t_box <- cut(time, breaks = object$t_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)
  u_box <- cut(duration, breaks = object$u_grid, labels = FALSE, right = FALSE, include.lowest = TRUE)

  pred <- rep(NA_real_, length(time))

  model <- object$models[[transition]]
  info <- object$model_info[[transition]]
  if (is.null(model)) {
    return(pred)
  }

  valid <- !is.na(t_box) & !is.na(u_box) & (duration <= time)
  if (any(valid)) {
    box <- paste0(t_box[valid], "_", u_box[valid])
    if (isTRUE(info$has_box_factor)) {
      nd <- data.frame(
        box_factor = factor(box, levels = info$box_levels),
        exposure = rep(1, length(box))
      )
      pred[valid] <- as.numeric(stats::predict(model, newdata = nd, type = "response"))
    } else {
      pred[valid] <- as.numeric(stats::predict(
        model,
        newdata = data.frame(exposure = rep(1, length(box))),
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
