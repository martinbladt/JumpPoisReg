library(JumpPoisReg)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
example_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

plot_dir <- file.path(example_dir, "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
plot_file <- function(name) file.path(plot_dir, name)

fill_missing_step <- function(x) {
  if (all(is.na(x))) {
    return(rep(0, length(x)))
  }
  first_ok <- which(!is.na(x))[1]
  x[seq_len(first_ok)] <- x[first_ok]
  for (i in seq_len(length(x))[-1]) {
    if (is.na(x[i])) {
      x[i] <- x[i - 1]
    }
  }
  x
}

mu_12 <- function(t) 0.09 + 0.0018 * t + 0.045 * sin(t / 2)
mu_13 <- function(t) 0.01 + 0.0002 * t + 0.005 * sin(t / 2)
mu_23 <- function(t) 0.06 + 0.0020 * t + 0.050 * sin(t / 2)

true_rate <- function(tr, t) {
  if (tr == "1->2") {
    mu_12(t)
  } else if (tr == "1->3") {
    mu_13(t)
  } else {
    mu_23(t)
  }
}

plot_markov_transition <- function(tr, fit, grid) {
  d <- fit$rates[fit$rates$transition == tr, c("interval", "rate")]
  d <- d[order(d$interval), ]
  d$rate <- fill_missing_step(d$rate)

  x_step <- c(grid[-length(grid)], tail(grid, 1))
  y_step <- c(d$rate, tail(d$rate, 1))
  t_seq <- seq(0, 40, length.out = 400)
  y_true <- true_rate(tr, t_seq)

  file_tag <- gsub("->", "", tr, fixed = TRUE)

  pdf(plot_file(paste0("poisson_markov_transition", file_tag, ".pdf")), width = 6, height = 6)
  plot(
    x_step,
    y_step,
    type = "s",
    col = "black",
    lwd = 2,
    xlab = "",
    ylab = "",
    ylim = c(0, 1.05 * max(c(y_step, y_true), na.rm = TRUE))
  )
  lines(t_seq, y_true, lty = 2, lwd = 2, col = "grey30")
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

# Increase n to 100000 for chapter-grade Monte Carlo precision.
n <- 20000
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
transitions <- c("1->2", "1->3", "2->3")
fit <- fit_markov_poisson(markov_sim, t_grid = grid, transitions = transitions)

for (tr in transitions) {
  plot_markov_transition(tr, fit, grid)
}

cat(
  "Generated Markov PDFs in", normalizePath(plot_dir),
  "for transitions", paste(transitions, collapse = ", "), "\n"
)
