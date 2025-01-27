source("R/header.R")

summary(cement)
cement <- rbind(cement, data.frame(x1 = 24, x2 = 24, x3 = 24, x4 = 24, y  = NA))

formula_1 <- y ~ 0 + 1 + x1 + x2 + x3 + x4

model_1_1 <- inla(formula = formula_1, 
                  family  = "gaussian",
                  data = cement,
                  control.predictor = list(compute = TRUE), # Should the marginals for the linear predictor be returned?
                  control.compute   = list(return.marginals.predictor = TRUE), # Should the marginals for the linear predictor be returned? 
                  control.fixed     = list(mean.intercept = 0, prec.intercept = 0, mean = 0, prec = 0.001)) # Priors

saveRDS(object = model_1_1, file = "models/model_1_1.RDS")

summary(model_1_1)

model_1_1$summary.fixed
model_1_1$summary.hyperpar
model_1_1$summary.linear.predictor
model_1_1$summary.fitted.values

# `x` represents the value of the parameter, and `y` is the density
model_1_1$marginals.fixed

# Explore posterior marginals
alpha <- model_1_1$marginals.fixed[[1]]

ggplot(data.frame(inla.smarginal(alpha)), aes(x, y)) +
  geom_line() +
  labs(x = "", y = "", title = "Posterior of alpha") + 
  custom_theme

qq <- inla.qmarginal(0.05, alpha)
qq

inla.pmarginal(qq, alpha)

inla.dmarginal(0, alpha)

sigma_2 <- inla.tmarginal(fun = function(x) { 1 / x }, marginal = model_1_1$marginals.hyperpar$`Precision for the Gaussian observations`)

ggplot(data.frame(inla.smarginal(sigma_2))) +
  geom_line(aes(x, y)) +
  labs(x = "", y = "", title = "Posterior of the variance") + 
  custom_theme

post_fitted <- model_1_1$marginals.fitted.values

post_margin <- data.frame(do.call(rbind, post_fitted))
post_margin$cement <- rep(names(post_fitted), times = sapply(post_fitted, nrow))

ggplot(post_margin) + 
  geom_line(aes(x, y)) +
  facet_wrap(~ cement, ncol = 5) +
  labs(x = "", y = "Density") +
  custom_theme

