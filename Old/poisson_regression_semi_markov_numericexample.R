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

prepare_OE <- function(semi_markov_sim, t_grid, u_grid) {
  lens <- vapply(semi_markov_sim, function(x) length(x$times), integer(1))
  path_frame <- data.frame(
    Id = rep(seq_along(lens), lens),
    Times = unlist(lapply(semi_markov_sim, `[[`, "times"), use.names = FALSE),
    States = unlist(lapply(semi_markov_sim, `[[`, "states"), use.names = FALSE)
  )

  path_frame$Dur <- with(
    path_frame,
    Times - ave(Times, Id, FUN = function(x) c(0, head(x, -1)))
  )

  path_frame$Trans_id <- with(
    path_frame,
    paste0(ave(States, Id, FUN = function(x) c(1, head(x, -1))), States)
  )

  path_frame <- path_frame[path_frame$Times != 0, ]
  path_frame$From_state <- substr(path_frame$Trans_id, 1, 1)
  path_frame$Entry_time <- path_frame$Times - path_frame$Dur

  b1 <- 1
  b2 <- 1
  t <- 0
  grid2d <- character(0)

  repeat {
    grid2d <- c(grid2d, paste0(b1, "_", b2))
    if (b1 == length(t_grid) - 1 && b2 == length(u_grid) - 1) {
      break
    }

    if (u_grid[b2 + 1] - t < t_grid[b1 + 1] - t) {
      t <- u_grid[b2 + 1]
      b2 <- b2 + 1
    } else {
      t <- t_grid[b1 + 1]
      b1 <- b1 + 1
      bx <- 1
      while (bx < b2) {
        grid2d <- c(grid2d, paste0(b1, "_", bx))
        bx <- bx + 1
      }
    }
  }

  path_frame$t_box <- cut(path_frame$Times, t_grid, labels = FALSE, right = FALSE)
  path_frame$u_box <- cut(path_frame$Dur, u_grid, labels = FALSE, right = FALSE)
  path_frame$Box <- paste0(path_frame$t_box, "_", path_frame$u_box)

  Inc_Expo <- function(entry, exit, tgrid, ugrid) {
    v <- c(
      tgrid[tgrid > entry & tgrid <= exit],
      entry + ugrid[entry + ugrid <= exit],
      exit
    )

    v <- sort(unique(v - entry))
    v <- diff(v)

    pathbox <- character(length(v))
    b1 <- cut(entry, tgrid, labels = FALSE, right = FALSE)
    b2 <- 1
    x <- entry
    y <- 0

    for (i in seq_along(v)) {
      pathbox[i] <- paste0(b1, "_", b2)
      if (i == length(v)) {
        break
      }
      if (tgrid[b1 + 1] - x == ugrid[b2 + 1] - y) {
        x <- tgrid[b1 + 1]
        y <- ugrid[b2 + 1]
        b1 <- b1 + 1
        b2 <- b2 + 1
      } else if (tgrid[b1 + 1] - x < ugrid[b2 + 1] - y) {
        y <- y + tgrid[b1 + 1] - x
        x <- tgrid[b1 + 1]
        b1 <- b1 + 1
      } else {
        x <- x + ugrid[b2 + 1] - y
        y <- ugrid[b2 + 1]
        b2 <- b2 + 1
      }
    }

    hit <- grid2d %in% pathbox
    hit[hit] <- v
    hit
  }

  df_expo <- matrix(0, nrow(path_frame), length(grid2d))
  colnames(df_expo) <- paste0("Box", grid2d)

  for (i in seq_len(nrow(path_frame))) {
    df_expo[i, ] <- Inc_Expo(path_frame$Entry_time[i], path_frame$Times[i], t_grid, u_grid)
  }

  path_frame <- cbind(path_frame, as.data.frame(df_expo))

  tmp <- path_frame[path_frame$From_state != path_frame$States, ]
  Occurences <- as.data.frame(
    table(tmp$Trans_id, tmp$From_state, tmp$Box),
    stringsAsFactors = FALSE
  )

  names(Occurences) <- c("Trans_id", "From_state", "Box", "O_jk_m")
  Occurences <- Occurences[Occurences$O_jk_m > 0, ]

  box_cols <- grep("^Box", names(path_frame), value = TRUE)
  X <- path_frame[, box_cols][, -1]
  X[] <- lapply(X, as.numeric)
  Exposure_wide <- rowsum(X, path_frame$From_state)

  Exposure <- data.frame(
    From_state = rownames(Exposure_wide),
    stack(Exposure_wide),
    row.names = NULL
  )

  names(Exposure)[2:3] <- c("E_j_m", "Box")
  Exposure$Box <- sub("^Box", "", Exposure$Box)

  OE <- merge(Occurences, Exposure, by = c("From_state", "Box"), all.x = TRUE)
  list(OE = OE, grid2d = grid2d)
}

fit_predict_Poisson <- function(OE, grid2d, t_grid, u_grid, t_vals, u_vals) {
  models <- list()

  for (jk in unique(OE$Trans_id)) {
    models[[jk]] <- glm(
      O_jk_m ~ Box + offset(log(E_j_m)),
      family = poisson,
      data = OE[OE$Trans_id == jk, ]
    )
  }

  pred_grid <- expand.grid(time = t_vals, duration = u_vals)
  pred_grid <- pred_grid[pred_grid$duration <= pred_grid$time, ]

  pred_grid$t_box <- cut(pred_grid$time, t_grid, labels = FALSE, right = FALSE)
  pred_grid$u_box <- cut(pred_grid$duration, u_grid, labels = FALSE, right = FALSE)
  pred_grid$Box <- paste0(pred_grid$t_box, "_", pred_grid$u_box)
  pred_grid <- pred_grid[pred_grid$Box %in% grid2d, ]

  add_fit <- function(jk) {
    model <- models[[jk]]
    valid_boxes <- unique(OE$Box[OE$Trans_id == jk])
    out <- pred_grid[pred_grid$Box %in% valid_boxes, ]
    out$E_j_m <- 1
    out$fit <- as.numeric(predict(model, newdata = out, type = "response"))
    out
  }

  pred_list <- lapply(names(models), add_fit)
  names(pred_list) <- names(models)
  pred_list
}

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

plot_surface_with_points <- function(filename, mu_fun, fits, t_vals, u_vals, theta, phi) {
  z_true <- build_surface(t_vals, u_vals, mu_fun)

  pdf(plot_file(filename), width = 6.0, height = 6.0)
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

  pdf(plot_file(filename), width = 6.0, height = 6.0)
  plot(obj$smooth_t, smooth_true, type = "l", col = "grey40", lwd = 1.5,
       xlab = "", ylab = "", ylim = y_lim)
  lines(obj$slice$time, obj$slice$fit, col = "black", lwd = 1.5, lty = 2, type = "s")
  dev.off()
}

plot_multi_slice <- function(filename, fits, mu_fun, target_diffs, cols,
                             y_lim = NULL, fig_width = 6.0, fig_height = 6.0) {
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

  legend("topleft",
         legend = paste0("t-u = ", diffs),
         col = cols,
         lty = 1,
         bty = "n",
         cex = 0.9)
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

res <- prepare_OE(semi_markov_sim, t_grid, u_grid)
OE <- res$OE
grid2d <- res$grid2d

t_vals <- seq(min(t_grid), max(t_grid), length.out = 80)
u_vals <- seq(min(u_grid), max(u_grid), length.out = 80)
pred_list <- fit_predict_Poisson(OE, grid2d, t_grid, u_grid, t_vals, u_vals)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition12.pdf",
  mu_fun = mu_12,
  fits = pred_list[["12"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition13.pdf",
  mu_fun = mu_13,
  fits = pred_list[["13"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition23_a.pdf",
  mu_fun = mu_23,
  fits = pred_list[["23"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -10,
  phi = 15
)

plot_surface_with_points(
  filename = "poisson_semi_markov_3d_transition23_b.pdf",
  mu_fun = mu_23,
  fits = pred_list[["23"]],
  t_vals = t_vals,
  u_vals = u_vals,
  theta = -65,
  phi = 20
)

slice_cols <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a")

plot_single_slice(
  filename = "poisson_semi_markov_slice12.pdf",
  fits = pred_list[["12"]],
  mu_fun = mu_12,
  target_diff = 0
)

plot_single_slice(
  filename = "poisson_semi_markov_slice13.pdf",
  fits = pred_list[["13"]],
  mu_fun = mu_13,
  target_diff = 0
)

plot_multi_slice(
  filename = "poisson_semi_markov_slice23.pdf",
  fits = pred_list[["23"]],
  mu_fun = mu_23,
  target_diffs = c(1, 5, 10, 20),
  cols = slice_cols,
  y_lim = c(0, 0.5),
  fig_width = 12.0,
  fig_height = 6.0
)

cat(
  "Generated semi-Markov 3D PDFs in", normalizePath(plot_dir),
  "for transitions 12, 13, and 23 (two perspectives), plus diagonal slice plots.\n"
)
