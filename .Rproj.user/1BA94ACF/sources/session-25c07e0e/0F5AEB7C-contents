## ----echo = FALSE, results = 'hide', echo = FALSE, warning = FALSE, message = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/INLAmethod-'
)
library(lattice) 

## ----eval = FALSE--------------------------------------------------------
## # Set CRAN mirror and INLA repository
## options(repos = c(getOption("repos"),
##   INLA = "https://inla.r-inla-download.org/R/testing"))
## # Install INLA and dependencies
## install.packages("INLA", dependencies = TRUE)

## ------------------------------------------------------------------------
library(INLA)
data(SPDEtoy)

## ----label = "SPDEtoy", results = 'hide', fig.cap =  '(ref:SPDEtoy)'-----

SPDEtoy.sp <- SPDEtoy
coordinates(SPDEtoy.sp) <- ~ s1 + s2

bubble(SPDEtoy.sp, "y", key.entries = c(5, 7.5, 10, 12.5, 15), 
       maxsize = 2, xlab = "s1", ylab = "s2")

## ------------------------------------------------------------------------
m0 <- inla(y ~ s1 + s2, data = SPDEtoy)

## ---- R.options=list(digits=2)-------------------------------------------
summary(m0)

## ----label = "SPDEtoym0",  echo = FALSE, results = 'hide', fig.cap = '(ref:SPDEtoym0)'----

par(mfrow = c(2, 2), mar=c(3, 3, 1, 1), mgp=c(2,1,0))

# Fixed effects
  plot(m0$marginals.fixed[[1]], type = "l", xlab = expression(alpha),
    ylab = "density")
  plot(m0$marginals.fixed[[2]], type = "l", xlab = expression(beta[1]),
    ylab = "density")
  plot(m0$marginals.fixed[[3]], type = "l", xlab = expression(beta[2]),
    ylab = "density")
# Precision
plot(m0$marginals.hyperpar[[1]], type = "l", xlab = expression(tau), 
  ylab = "density")


## ------------------------------------------------------------------------
f.rw1 <- y ~ f(s1, model = "rw1", scale.model = TRUE) +
  f(s2, model = "rw1", scale.model = TRUE)

## ---- R.options=list(digits=2)-------------------------------------------
m1 <- inla(f.rw1, data = SPDEtoy)

summary(m1)

## ----label = "SPDEtoym1",  echo = FALSE, results = 'hide', fig.cap =  '(ref:SPDEtoym1)'----

par(mfrow = c(3, 2), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))

# Plot RW1 effect
# X: model$summary.random components to plot.
plot.rw1 <- function(X, ...) {
  plot(X[1:2], type = "l", 
       ylim=range(X[, c(4,6)]),  ...)
  lines(X$ID, X[,4], lty = 2)
  lines(X$ID, X[,6], lty = 2)
}

# Intercept
  plot(m1$marginals.fixed[[1]], type = "l", xlab = expression(alpha),
    ylab = "density")
# Precision
plot(m1$marginals.hyperpar[[1]], type = "l", xlab = expression(tau),
  ylab = "density")

# Non-linear effects
plot.rw1(m1$summary.random$s1, xlab = expression(s[1][",i"]),
  ylab = expression(u[1][",i"]))
plot.rw1(m1$summary.random$s2, xlab = expression(s[2][",i"]),
  ylab = expression(u[2][",i"]))

# Precisions on non-linear effects
plot(m1$marginals.hyperpar[[2]], type = "l", xlab = expression(tau[1]),
  ylab = "density", xlim = c(0, 75))
plot(m1$marginals.hyperpar[[3]], type = "l", xlab = expression(tau[2]),
  ylab = "density", xlim = c(0, 75))

## ---- R.options=list(digits=3)-------------------------------------------
m0$summary.fixed

## ----eval = FALSE--------------------------------------------------------
## plot(m0$marginals.fixed[[1]], type = "l",
##   xlab = expression(alpha), ylab = "density")

## ------------------------------------------------------------------------
SPDEtoy.pred <- rbind(SPDEtoy, c(NA, 0.5, 0.5))

## ------------------------------------------------------------------------
m0.pred <- inla(y ~ s1 + s2, data = SPDEtoy.pred,
  control.predictor = list(compute = TRUE))

## ----eval = FALSE--------------------------------------------------------
## m0.pred$marginals.fitted.values[[201]]

## ----label = "inlamissing", echo = FALSE, results = 'hide', out.width="69%", fig.cap =  "Predictive distribution of the response at location (0.5, 0.5)."----

plot(m0.pred$marginals.fitted.values[[201]], type = "l", 
     xlab = expression(y[201]), ylab = "density")

## ------------------------------------------------------------------------
m0.opts <- inla(y ~ s1 + s2, data = SPDEtoy,
  control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE)
)

## ---- R.options=list(digits=2)-------------------------------------------
summary(m0.opts)

## ----eval = FALSE--------------------------------------------------------
## m0.opts$cpo$cpo
## m0.opts$cpo$pit

## ----label = "cpopit", echo = FALSE, results = 'hide', fig.width = 10, fig.height = 5, fig.cap =  "Histograms of CPO and PIT values for the model with fixed effects."----

par(mfrow = c(1, 2), mar=c(3,3,2,1), mgp=c(2,1,0))

hist(m0.opts$cpo$cpo, main = "CPO", xlab = "cpo")
hist(m0.opts$cpo$pit, main = "PIT", xlab = "pit")

## ------------------------------------------------------------------------
m1.ccd <- inla(f.rw1, data = SPDEtoy,
  control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla = list(int.strategy = "ccd"))
m1.eb <- inla(f.rw1, data = SPDEtoy,
  control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
  control.inla = list(int.strategy = "eb"))

## ------------------------------------------------------------------------
# CCD strategy
m1.ccd$cpu.used

# EB strategy
m1.eb$cpu.used

## ------------------------------------------------------------------------
# Prior on the fixed effects
prior.fixed <- list(mean.intercept = 0, prec.intercept = 1,
  mean = 0, prec = 1)

## ------------------------------------------------------------------------
# Prior on the likelihood precision (log-scale)
prior.prec <- list(initial = 0, prior = "normal", param = c(0, 1),
  fixed = FALSE)

## ------------------------------------------------------------------------
# Prior on the precision of the RW1
prior.rw1 <- list(initial = 0, fixed = TRUE)

## ------------------------------------------------------------------------
f.hyper <- y ~ 1 +
  f(s1, model = "rw1", hyper = list(prec = prior.rw1),
    scale.model = TRUE) +
  f(s2, model = "rw1", hyper = list(prec = prior.rw1),
    scale.model = TRUE)

m1.hyper <- inla(f.hyper, data = SPDEtoy,
  control.fixed = prior.fixed,
  control.family = list(hyper = list(prec = prior.prec)))

## ------------------------------------------------------------------------
summary(m1.hyper)

## ------------------------------------------------------------------------
# Compute posterior marginal of variance
post.var <- inla.tmarginal(function(x) exp(-x), 
  m0$internal.marginals.hyperpar[[1]])

## ------------------------------------------------------------------------
# Compute summary statistics
inla.zmarginal(post.var)

## ------------------------------------------------------------------------
inla.hpdmarginal(0.95, post.var)

