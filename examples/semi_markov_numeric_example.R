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

mu_12 <- function(t, u) 0.09 + 0.0018 * t
mu_13 <- function(t, u) 0.01 + 0.0002 * t
mu_23 <- function(t, u) 0.09 + 0.001 * t * (1 + 0.1 * u) + 0.2 / (1 + exp(0.5 * (u - 4)))

build_surface <- function(t_vals, u_vals, mu_fun) {
  z_true <- outer(t_vals, u_vals, Vectorize(mu_fun))
  for (i in seq_along(t_vals)) {
    for (j in seq_along(u_vals)) {
      if (u_vals[j] > t_vals[i]) {
        z_true[i, j] <- NA
      }
    }
  }
  z_true
}

build_prediction_cloud <- function(fit, transition, t_vals, u_vals) {
  pred_grid <- expand.grid(time = t_vals, duration = u_vals)
  pred_grid <- pred_grid[pred_grid$duration <= pred_grid$time, ]
  pred_grid$fit <- predict_semi_markov_poisson(
    fit,
    time = pred_grid$time,
    duration = pred_grid$duration,
    transition = transition
  )
  pred_grid[!is.na(pred_grid$fit), ]
}

plot_surface_with_points <- function(filename, mu_fun, fits, t_vals, u_vals, theta, phi) {
  z_true <- build_surface(t_vals, u_vals, mu_fun)

  pdf(plot_file(filename), width = 6, height = 6)
  par(mar = c(3.6, 3.6, 1.0, 1.4))

  pmat <- persp(
    x = t_vals,
    y = u_vals,
    z = z_true,
    theta = theta,
    phi = phi,
    expand = 1,
    col = "lightblue",
    axes = TRUE,
    ticktype = "detailed",
    d = 100,
    xlab = "Time",
    ylab = "Duration",
    zlab = "Intensity",
    main = ""
  )

  pts <- trans3d(fits$time, fits$duration, fits$fit, pmat = pmat)
  points(pts, pch = 16, col = rgb(0.75, 0, 0, 0.45), cex = 0.5)
  dev.off()
}

build_diagonal_slice <- function(fits, target_diff, tol = 0.5) {
  slice <- fits[abs((fits$time - fits$duration) - target_diff) < tol, ]
  slice <- slice[order(slice$time), ]

  if (nrow(slice) < 2) {
    return(NULL)
  }

  smooth_t <- seq(min(slice$time), max(slice$time), length.out = 500)
  smooth_u <- smooth_t - target_diff
  list(slice = slice, smooth_t = smooth_t, smooth_u = smooth_u)
}

plot_single_slice <- function(filename, fits, mu_fun, target_diff) {
  obj <- build_diagonal_slice(fits, target_diff)
  if (is.null(obj)) {
    return(invisible(NULL))
  }

  smooth_true <- mu_fun(obj$smooth_t, obj$smooth_u)
  y_lim <- range(c(obj$slice$fit, smooth_true), na.rm = TRUE)

  pdf(plot_file(filename), width = 6, height = 6)
  plot(obj$smooth_t, smooth_true, type = "l", col = "grey40", lwd = 1.5,
       xlab = "", ylab = "", ylim = y_lim)
  lines(obj$slice$time, obj$slice$fit, col = "black", lwd = 1.5, lty = 2, type = "s")
  dev.off()
}

plot_multi_slice <- function(filename, fits, mu_fun, target_diffs, cols,
                             y_lim = NULL, fig_width = 6, fig_height = 6) {
  objs <- lapply(target_diffs, function(td) build_diagonal_slice(fits, td))
  keep <- vapply(objs, function(x) !is.null(x), logical(1))
  objs <- objs[keep]
  diffs <- target_diffs[keep]
  cols <- cols[seq_along(objs)]

  if (length(objs) == 0) {
    return(invisible(NULL))
  }

  x_min <- min(vapply(objs, function(x) min(x$smooth_t), numeric(1)))
  x_max <- max(vapply(objs, function(x) max(x$smooth_t), numeric(1)))

  y_all <- c()
  for (i in seq_along(objs)) {
    y_all <- c(y_all, objs[[i]]$slice$fit, mu_fun(objs[[i]]$smooth_t, objs[[i]]$smooth_u))
  }
  auto_y_lim <- range(y_all, na.rm = TRUE)
  if (is.null(y_lim)) {
    y_lim <- auto_y_lim
  }

  pdf(plot_file(filename), width = fig_width, height = fig_height)
  plot(NA, NA, type = "n", xlim = c(x_min, x_max), ylim = y_lim, xlab = "", ylab = "")

  for (i in seq_along(objs)) {
    obj <- objs[[i]]
    y_true <- mu_fun(obj$smooth_t, obj$smooth_u)
    lines(obj$smooth_t, y_true, col = cols[i], lwd = 1.5, lty = 1)
    lines(obj$slice$time, obj$slice$fit, col = cols[i], lwd = 1.5, lty = 2, type = "s")
  }

  legend("topleft", legend = paste0("t-u = ", diffs), col = cols, lty = 1, bty = "n", cex = 0.9)
  dev.off()
}

set.seed(3)

jump_rate <- function(i, t, u) {
  if (i == 1) {
    0.1 + 0.002 * t
  } else if (i == 2) {
    0.09 + 0.001 * t * (1 + 0.1 * u) +
      0.2 / (1 + exp(0.5 * (u - 4)))
  } else {
    0
  }
}

mark_dist <- function(i, s, v) {
  if (i == 1) {
    c(0, 0.9, 0.1)
  } else if (i == 2) {
    c(0, 0, 1)
  } else {
    0
  }
}

# Increase n to 100000 for chapter-grade Monte Carlo precision.
n <- 100000
censoring <- runif(n, 10, 40)
semi_markov_sim <- vector("list", n)

for (i in seq_len(n)) {
  semi_markov_sim[[i]] <- sim_path(
    1,
    rates = jump_rate,
    dists = mark_dist,
    tn = censoring[i],
    bs = c(
      0.1 + 0.002 * censoring[i],
      0.29 + 0.001 * censoring[i],
      0
    )
  )
}

t_grid <- seq(0, 40, by = 2)
u_grid <- seq(0, 40, by = 2)
transitions <- c("1->2", "1->3", "2->3")

fit <- fit_semi_markov_poisson(
  semi_markov_sim,
  t_grid = t_grid,
  u_grid = u_grid,
  transitions = transitions
)

t_vals <- seq(min(t_grid), max(t_grid), length.out = 80)
u_vals <- seq(min(u_grid), max(u_grid), length.out = 80)

pred_list <- list(
  "1->2" = build_prediction_cloud(fit, "1->2", t_vals, u_vals),
  "1->3" = build_prediction_cloud(fit, "1->3", t_vals, u_vals),
  "2->3" = build_prediction_cloud(fit, "2->3", t_vals, u_vals)
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition12.pdf",
  mu_fun = mu_12,
  fits = pred_list[["1->2"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition13.pdf",
  mu_fun = mu_13,
  fits = pred_list[["1->3"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition23_a.pdf",
  mu_fun = mu_23,
  fits = pred_list[["2->3"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition23_b.pdf",
  mu_fun = mu_23,
  fits = pred_list[["2->3"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -65,
  phi = 20
)

slice_cols <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")

plot_single_slice(
  filename = "poisson_semi_markov_slice12.pdf",
  fits = pred_list[["1->2"]],
  mu_fun = mu_12,
  target_diff = 0
)

plot_single_slice(
  filename = "poisson_semi_markov_slice13.pdf",
  fits = pred_list[["1->3"]],
  mu_fun = mu_13,
  target_diff = 0
)

plot_multi_slice(
  filename = "poisson_semi_markov_slice23.pdf",
  fits = pred_list[["2->3"]],
  mu_fun = mu_23,
  target_diffs = c(1, 5, 10, 20),
  cols = slice_cols,
  y_lim = c(0, 0.5),
  fig_width = 12,
  fig_height = 6
)

cat(
  "Generated semi-Markov PDFs in", normalizePath(plot_dir),
  "for transitions 1->2, 1->3, and 2->3.\n"
)
