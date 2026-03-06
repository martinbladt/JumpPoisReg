sim_jump <- function(i, t, u = 0, tn, rate, dist, b = NA) {
  if (is.na(b)) {
    b <- -stats::optimize(function(x) -rate(t + x, u + x), c(0, tn - t))$objective
  }

  r <- stats::rexp(1, rate = b)
  s <- t + r
  v <- u + r
  y <- stats::runif(1)

  while (y > rate(s, v) / b) {
    r <- stats::rexp(1, rate = b)
    s <- s + r
    v <- v + r
    y <- stats::runif(1)
  }

  pr <- dist(s, v)
  j <- sample(length(pr), 1, prob = pr)
  list(time = s, mark = j)
}

#' Simulate a jump-process path
#'
#' Simulates the path of a time-inhomogeneous Markov or semi-Markov jump
#' process up to a maximum time.
#'
#' @param i Initial state (integer).
#' @param rates Function `rates(state, time, duration)` returning total jump rate.
#' @param dists Function `dists(state, time, duration)` returning next-state probabilities.
#' @param t Initial time. Default is `0`.
#' @param u Initial duration since the last jump. Default is `0`.
#' @param tn Maximum simulation time. Default is `Inf`.
#' @param abs Logical vector indicating absorbing states. If omitted, the last
#'   state is treated as absorbing.
#' @param bs Optional upper bounds for rates by state. If omitted, bounds are
#'   approximated by optimization.
#'
#' @return A list with elements `times` and `states`.
#' @export
#'
#' @examples
#' jump_rate <- function(i, t, u) if (i == 1) 2 + t else if (i == 2) 1 + t else 0
#' mark_dist <- function(i, s, v) if (i == 1) c(0, 0.7, 0.3) else if (i == 2) c(0, 0, 1) else 0
#' sim_path(1, rates = jump_rate, dists = mark_dist, tn = 2)
sim_path <- function(i, rates, dists, t = 0, u = 0, tn = Inf, abs = numeric(0), bs = NA) {
  times <- t
  marks <- i

  if (length(abs) == 0) {
    abs <- c(rep(FALSE, length(dists(i, t, u)) - 1), TRUE)
  }

  while (!abs[utils::tail(marks, 1)]) {
    z <- sim_jump(
      utils::tail(marks, 1),
      utils::tail(times, 1),
      u,
      tn,
      function(s, v) rates(utils::tail(marks, 1), s, v),
      function(s, v) dists(utils::tail(marks, 1), s, v),
      bs[utils::tail(marks, 1)]
    )

    if (z$time > tn) {
      break
    }

    times <- c(times, z$time)
    marks <- c(marks, z$mark)
    u <- 0
  }

  if (!abs[utils::tail(marks, 1)]) {
    times <- c(times, tn)
    marks <- c(marks, utils::tail(marks, 1))
  }

  list(times = times, states = marks)
}
