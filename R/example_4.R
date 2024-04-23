source("R/header.R")

data_USA <- readRDS(file = "data/example_4/data_USA.rds")
USA <- readRDS(file = "data/example_4/USA_filtered.rds")

ggplot() + 
  geom_sf(data = USA, fill = "white") +
  geom_sf(data = data_USA, aes(fill = mean), color = "black", size = 3, shape = 21) +
  scale_fill_gradientn(name = "PM2.5 level", colors = rainbow(9, start = 0.1, end = 0.9)) + 
  labs(x = "", y = "", title = "") +
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks  = element_blank())


# Boundary coordinates
USA_coor <- sf::st_coordinates(USA)
USA_coor <- matrix(c(USA_coor[, 1], USA_coor[, 2]), ncol = 2)
colnames(USA_coor) <- c("lon", "lat")

# Points coordinates
data_coor <- sf::st_coordinates(data_USA)
data_USA  <- bind_cols(data_USA, as_tibble(data_coor))
data_USA  <- data_USA %>% rename(lon = X, lat = Y) %>% dplyr::select(mean, sd, lon, lat, geometry)

# `max.edge`: the largest allowed triangle edge length. One or two values.
# `offset`: the automatic extension distance. One or two values, for an inner and an optional outer extension.

mesh <- inla.mesh.2d(loc.domain = USA_coor, max.edge = c(300, 3000), offset = c(300, 1500))
saveRDS(object = mesh, file = "data/example_4/mesh.RDS")

mesh$n # Number of nodes

# Plot `mesh`
{ 
  plot(mesh)
  plot(USA$geometry, lwd = 2, border = "red", add = TRUE)
  points(data_USA$lon, data_USA$lat, pch = 1, col = "green") 
}


# Building the SPDE model
# alpha = nu + d / 2 = 1 + 1 = 2, for nu = 1

# Flexible parameterization (with a default one)
# spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
# Parameterization for PC priors
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            alpha = 2,
                            prior.range = c(1e3, 0.90), # P(range < 1e3) = 0.90
                            prior.sigma = c(1.0, 0.01)) # P(sigma > 1.0) = 0.01

# Indices
indxs <- inla.spde.make.index("s", spde$n.spde)
# Projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = data_coor)

# Prediction
coord_pred <- create_prediction_grid(USA)
saveRDS(object = coord_pred, file = "data/example_4/coord_pred.RDS")
{
  coord_pred_cp <- as.data.frame(coord_pred)
  coordinates(coord_pred_cp) <- ~ x + y
  plot(coord_pred_cp)
}

# Projection matrix for prediction
Ap <- inla.spde.make.A(mesh = mesh, loc = coord_pred)

# Create stacks
# Stack for estimation 
stk_e <- inla.stack(tag = "est",
                    data = list(y = data_USA$mean),
                    A = list(1, A),
                    effects = list(data.frame(b0 = rep(1, nrow(data_USA))), s = indxs))

# Stack for prediction
stk_p <- inla.stack(tag = "pred",
                    data = list(y = NA),
                    A = list(1, Ap),
                    effects = list(data.frame(b0 = rep(1, nrow(coord_pred))), s = indxs))

# Full stack
stk_full <- inla.stack(stk_e, stk_p)

# Fit model
formula_1 <- y ~ 0 + b0 + f(s, model = spde)

model_4_1 <- inla(formula = formula_1,
                  family  = "gaussian", 
                  data = inla.stack.data(stk_full), 
                  control.predictor = list(compute = TRUE,
                                           A = inla.stack.A(stk_full))) # Matrix of predictors

saveRDS(object = model_4_1, file = "models/model_4_1.RDS")

summary(model_4_1)

# Posteriors

ss <- inla.tmarginal(fun = function(x) { 1 / sqrt(x) }, marginal = model_4_1$marginals.hyperpar$`Precision for the Gaussian observations`)

ggplot(data.frame(inla.smarginal(ss))) +
  geom_line(aes(x, y)) +
  labs(x = "", y = "", title = "Posterior of the standard deviation for the observations") + 
  custom_theme

ggplot(data.frame(inla.smarginal(model_4_1$marginals.hyperpar$`Range for s`))) +
  geom_line(aes(x, y)) +
  labs(x = "", y = "", title = "Posterior of the range in the Matérn model") + 
  custom_theme

ggplot(data.frame(inla.smarginal(model_4_1$marginals.hyperpar$`Stdev for s`))) +
  geom_line(aes(x, y)) +
  labs(x = "", y = "", title = "Posterior of the standard deviation in the Matérn model") + 
  custom_theme

# According to this parameterization: https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html#sec:matern
# range = sqrt(8 * smoothness) / scale

# Fitted values and prediction
idxs_pred <- inla.stack.index(stk_full, tag = "pred")$data

pred_mm <- as.data.frame(cbind(coord_pred, model_4_1$summary.fitted.values[idxs_pred, "mean"]))
pred_ll <- as.data.frame(cbind(coord_pred, model_4_1$summary.fitted.values[idxs_pred, "0.025quant"]))
pred_uu <- as.data.frame(cbind(coord_pred, model_4_1$summary.fitted.values[idxs_pred, "0.975quant"]))

r_mm <- pred_mm$V3 %>% range()
r_ll <- pred_ll$V3 %>% range()
r_uu <- pred_uu$V3 %>% range()
r <- c(min(r_mm[1], r_ll[1], r_uu[1]), max(r_mm[2], r_ll[2], r_uu[2]))

pp_mm <- plot_pred_USA(fitted_values = pred_mm, USA = USA, r = r, tt = "Mean")
pp_ll <- plot_pred_USA(fitted_values = pred_ll, USA = USA, r = r, tt = "2.5th")
pp_uu <- plot_pred_USA(fitted_values = pred_uu, USA = USA, r = r, tt = "97.5th")

(pp_ll + pp_mm + pp_uu) + plot_layout(guides = "collect") & theme(legend.position = "right")
