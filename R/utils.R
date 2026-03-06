.validate_grid <- function(grid, name) {
  if (!is.numeric(grid) || length(grid) < 2 || anyNA(grid)) {
    stop(sprintf("`%s` must be a numeric vector with at least two entries.", name))
  }
  if (any(diff(grid) <= 0)) {
    stop(sprintf("`%s` must be strictly increasing.", name))
  }
  invisible(TRUE)
}

.validate_paths <- function(paths) {
  if (!is.list(paths) || length(paths) == 0) {
    stop("`paths` must be a non-empty list of trajectories.")
  }

  for (idx in seq_along(paths)) {
    path <- paths[[idx]]
    if (!is.list(path) || is.null(path$times) || is.null(path$states)) {
      stop(sprintf("Path %d must be a list containing `times` and `states`.", idx))
    }

    times <- as.numeric(path$times)
    states <- as.integer(path$states)

    if (length(times) != length(states) || length(times) < 2) {
      stop(sprintf("Path %d must have equal-length `times`/`states` with length >= 2.", idx))
    }
    if (anyNA(times) || anyNA(states)) {
      stop(sprintf("Path %d contains missing values.", idx))
    }
    if (any(diff(times) < 0)) {
      stop(sprintf("Path %d has decreasing jump times.", idx))
    }
  }

  invisible(TRUE)
}

.flatten_paths <- function(paths) {
  .validate_paths(paths)

  out <- vector("list", length(paths))
  for (idx in seq_along(paths)) {
    times <- as.numeric(paths[[idx]]$times)
    states <- as.integer(paths[[idx]]$states)

    out[[idx]] <- data.frame(
      id = idx,
      entry = utils::head(times, -1),
      exit = utils::tail(times, -1),
      from = utils::head(states, -1),
      to = utils::tail(states, -1),
      stringsAsFactors = FALSE
    )
  }

  segments <- do.call(rbind, out)
  rownames(segments) <- NULL
  segments$duration <- segments$exit - segments$entry

  if (any(segments$duration < 0)) {
    stop("At least one path has negative segment duration.")
  }

  segments
}

.format_transition <- function(from, to) {
  paste0(from, "->", to)
}

.parse_transition_string <- function(x) {
  x <- gsub("\\s+", "", x)

  if (grepl("^[0-9]+->[0-9]+$", x)) {
    parts <- strsplit(x, "->", fixed = TRUE)[[1]]
    return(c(from = as.integer(parts[1]), to = as.integer(parts[2])))
  }

  # Backward-compatible shorthand like "12" for 1->2.
  if (grepl("^[0-9]{2}$", x)) {
    return(c(from = as.integer(substr(x, 1, 1)), to = as.integer(substr(x, 2, 2))))
  }

  stop(
    "Transition labels must be in format 'i->j' (or shorthand 'ij' for single-digit states)."
  )
}

.normalize_transitions <- function(transitions, segments) {
  if (is.null(transitions)) {
    jumps <- unique(segments[segments$from != segments$to, c("from", "to"), drop = FALSE])
    if (nrow(jumps) == 0) {
      stop("No observed jumps found in `paths`.")
    }
    trans_df <- jumps
  } else if (is.matrix(transitions) || is.data.frame(transitions)) {
    if (ncol(transitions) < 2) {
      stop("`transitions` matrix/data.frame must have two columns: from, to.")
    }
    trans_df <- data.frame(
      from = as.integer(transitions[, 1]),
      to = as.integer(transitions[, 2]),
      stringsAsFactors = FALSE
    )
  } else if (is.character(transitions)) {
    parsed <- lapply(transitions, .parse_transition_string)
    trans_df <- data.frame(
      from = as.integer(vapply(parsed, `[[`, numeric(1), "from")),
      to = as.integer(vapply(parsed, `[[`, numeric(1), "to")),
      stringsAsFactors = FALSE
    )
  } else {
    stop("`transitions` must be NULL, character, matrix, or data.frame.")
  }

  if (anyNA(trans_df$from) || anyNA(trans_df$to)) {
    stop("`transitions` contains invalid state labels.")
  }
  if (any(trans_df$from == trans_df$to)) {
    stop("Transitions must satisfy `from != to`.")
  }

  trans_df <- unique(trans_df)
  trans_df <- trans_df[order(trans_df$from, trans_df$to), , drop = FALSE]
  trans_df$transition <- .format_transition(trans_df$from, trans_df$to)
  trans_df[, c("transition", "from", "to")]
}

.box_to_indices <- function(box) {
  split_box <- strsplit(box, "_", fixed = TRUE)
  data.frame(
    t_box = as.integer(vapply(split_box, `[[`, character(1), 1)),
    u_box = as.integer(vapply(split_box, `[[`, character(1), 2)),
    stringsAsFactors = FALSE
  )
}
