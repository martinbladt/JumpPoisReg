library(AalenJohansen)

plot_dir <- if (dir.exists("Plots")) {
  "Plots"
} else if (dir.exists("../Plots")) {
  "../Plots"
} else if (dir.exists("../../Plots")) {
  "../../Plots"
} else if (dir.exists("../../../Plots")) {
  "../../../Plots"
} else {
  stop("Could not locate the 'Plots' directory.")
}
plot_file <- function(name) file.path(plot_dir, name)

prepare_oe_markov <- function(paths, grid, transitions) {
  n_intervals <- length(grid) - 1

  len <- vapply(paths, function(x) length(x$times), integer(1))
  frame <- data.frame(
    id = rep(seq_along(paths), times = len),
    time = unlist(lapply(paths, `[[`, "times"), use.names = FALSE),
    state = unlist(lapply(paths, `[[`, "states"), use.names = FALSE)
  )
  frame <- frame[order(frame$id, frame$time), ]

  frame$lag_time <- ave(frame$time, frame$id, FUN = function(x) c(0, x[-length(x)]))
  frame$lag_state <- ave(frame$state, frame$id, FUN = function(x) c(0, x[-length(x)]))
  frame$entry <- frame$lag_time
  frame$from <- frame$lag_state
  frame$to <- frame$state
  frame$trans <- paste0(frame$from, frame$to)

  proc <- frame[frame$time > 0, ]

  for (m in seq_len(n_intervals)) {
    interval_start <- grid[m]
    interval_end <- grid[m + 1]
    overlap_start <- pmax(proc$entry, interval_start)
    overlap_end <- pmin(proc$time, interval_end)
    proc[[paste0("I_", m)]] <- pmax(0, overlap_end - overlap_start)
  }

  jump_rows <- proc[proc$from != proc$to, ]
  jump_rows$int_id <- cut(jump_rows$time, breaks = grid, labels = FALSE, right = TRUE)
  counts <- aggregate(id ~ trans + int_id, data = jump_rows, length)
  names(counts)[names(counts) == "id"] <- "O"

  exposure <- aggregate(proc[paste0("I_", seq_len(n_intervals))],
                        by = list(from = proc$from),
                        FUN = sum)
  exp_long <- reshape(exposure,
                      direction = "long",
                      varying = paste0("I_", seq_len(n_intervals)),
                      v.names = "E",
                      timevar = "int_id",
                      times = seq_len(n_intervals))
  exp_long <- exp_long[, c("from", "int_id", "E")]
  rownames(exp_long) <- NULL

  design <- expand.grid(
    trans = transitions,
    int_id = seq_len(n_intervals),
    stringsAsFactors = FALSE
  )
  design$from <- as.integer(substr(design$trans, 1, 1))

  oe <- merge(design, exp_long, by = c("from", "int_id"), all.x = TRUE)
  oe <- merge(oe, counts, by = c("trans", "int_id"), all.x = TRUE)
  oe$O[is.na(oe$O)] <- 0
  oe$E[is.na(oe$E)] <- 0
  oe
}

fit_poisson_markov <- function(oe, transitions, n_intervals) {
  fits <- data.frame()

  for (tr in transitions) {
    sub <- oe[oe$trans == tr & oe$E > 0, ]
    sub$int_id <- factor(sub$int_id)
    model <- glm(O ~ int_id + offset(log(E)), family = poisson(), data = sub)

    fitted_rates <- predict(model, type = "response") / sub$E
    ids <- as.integer(as.character(sub$int_id))
    fits <- rbind(fits, data.frame(trans = tr, int_id = ids, rate = fitted_rates))
  }

  complete <- expand.grid(trans = transitions, int_id = seq_len(n_intervals), stringsAsFactors = FALSE)
  complete <- merge(complete, fits, by = c("trans", "int_id"), all.x = TRUE)
  complete
}

fill_missing_step <- function(x) {
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  first_ok <- which(!is.na(x))[1]
  x[seq_len(first_ok)] <- x[first_ok]
  for (i in seq_len(length(x))[-1]) {
    if (is.na(x[i])) x[i] <- x[i - 1]
  }
  x
}

mu_12 <- function(t) 0.09 + 0.0018 * t + 0.045 * sin(t / 2)
mu_13 <- function(t) 0.01 + 0.0002 * t + 0.005 * sin(t / 2)
mu_23 <- function(t) 0.06 + 0.0020 * t + 0.050 * sin(t / 2)

true_rate <- function(tr, t) {
  if (tr == "12") {
    mu_12(t)
  } else if (tr == "13") {
    mu_13(t)
  } else {
    mu_23(t)
  }
}

plot_markov_transition <- function(tr, fit_rates, grid) {
  d <- fit_rates[fit_rates$trans == tr, ]
  d <- d[order(d$int_id), ]
  d$rate <- fill_missing_step(d$rate)

  x_step <- c(grid[-length(grid)], tail(grid, 1))
  y_step <- c(d$rate, tail(d$rate, 1))
  t_seq <- seq(0, 40, length.out = 400)
  y_true <- true_rate(tr, t_seq)

  pdf(plot_file(paste0("poisson_markov_transition", tr, ".pdf")),
      width = 6.0, height = 6.0)

  plot(x_step, y_step, type = "s", col = "black", lwd = 2.0,
       xlab = "", ylab = "",
       ylim = c(0, 1.05 * max(c(y_step, y_true), na.rm = TRUE)))
  lines(t_seq, y_true, lty = 2, lwd = 2.0, col = "grey30")
  abline(v = grid, lty = 3, col = "grey80", lwd = 0.8)

  dev.off()
}

set.seed(1)

jump_rate_markov <- function(i, t, u) {
  if (i == 1) {
    0.1 + 0.002 * t + 0.05 * sin(t / 2)
  } else if (i == 2) {
    0.06 + 0.002 * t + 0.05 * sin(t / 2)
  } else {
    0
  }
}

mark_dist_markov <- function(i, s, v) {
  if (i == 1) {
    c(0, 0.9, 0.1)
  } else if (i == 2) {
    c(0, 0, 1)
  } else {
    0
  }
}

n <- 100000
censoring <- runif(n, 10, 40)
markov_sim <- vector("list", n)

for (i in seq_len(n)) {
  markov_sim[[i]] <- sim_path(
    1,
    rates = jump_rate_markov,
    dists = mark_dist_markov,
    tn = censoring[i],
    bs = c(
      0.1 + 0.002 * censoring[i] + 0.05,
      0.06 + 0.002 * censoring[i] + 0.05,
      0
    )
  )
}

grid <- seq(0, 40, length.out = 40)
transitions <- c("12", "13", "23")

oe <- prepare_oe_markov(markov_sim, grid, transitions)
fit_rates <- fit_poisson_markov(oe, transitions, n_intervals = length(grid) - 1)

for (tr in transitions) {
  plot_markov_transition(tr, fit_rates, grid)
}

cat("Generated Markov PDFs in", normalizePath(plot_dir), "for transitions",
    paste(transitions, collapse = ", "), "\n")
