# JumpPoisReg

`JumpPoisReg` is a minimal package for Poisson regression of transition intensities from jump-path data.

## Main functions

- `fit_markov_poisson()`: piecewise-constant Markov intensity fit on a time grid.
- `predict_markov_poisson()`: predicts fitted Markov intensities for new times.
- `fit_semi_markov_poisson()`: piecewise-constant semi-Markov intensity fit on a time-duration grid.
- `predict_semi_markov_poisson()`: predicts fitted semi-Markov intensities for `(time, duration)` pairs.
- `sim_path()`: simulates jump paths (ported from the `AalenJohansen` package template).

## Example scripts

See `examples/markov_numeric_example.R` and `examples/semi_markov_numeric_example.R`.
