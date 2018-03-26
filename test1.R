m=4500
sigmaPrior= 0.01
betaPrior = rnorm(m, 1, sigmaPrior)
hist(betaPrior, freq = FALSE)
lines(density(betaPrior))
