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

mesh <- readRDS(file = "data/example_4/mesh.RDS")
dual_mesh <- readRDS(file = "data/example_5/dual_mesh.RDS")
wgt <- readRDS(file = "data/example_5/wgt.RDS")

# Priors might differ
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            alpha = 2,
                            prior.range = c(1e3, 0.90), # P(range < 1e3) = 0.90
                            prior.sigma = c(1.0, 0.01)) # P(sigma > 1.0) = 0.01

#########################
#########################

# As before...

# Points coordinates
data_coor <- sf::st_coordinates(data_USA)
data_USA  <- bind_cols(data_USA, as_tibble(data_coor))
data_USA  <- data_USA %>% rename(lon = X, lat = Y) %>% dplyr::select(mean, sd, lon, lat, geometry)


# Boundary coordinates
USA_coor <- sf::st_coordinates(USA)
USA_coor <- matrix(c(USA_coor[, 1], USA_coor[, 2]), ncol = 2)
colnames(USA_coor) <- c("lon", "lat")

#########################
#########################

# Data and projection matrices

n_vtx <- mesh$n
n_pts <- nrow(data_USA)

coord_pred <- readRDS(file = "data/example_4/coord_pred.RDS")
n_pts_pred <- nrow(coord_pred)

indxs <- inla.spde.make.index("s", spde$n.spde)
indxv <- inla.spde.make.index("v", spde$n.spde)

y_pp <- rep(0:1, c(n_vtx, n_pts))
e_pp <- c(wgt, rep(0, n_pts)) 
imat <- Diagonal(n_vtx, rep(1, n_vtx))
ymat <- inla.spde.make.A(mesh, data_coor)
A_pp <- rbind(imat, ymat)
# Prediction
A_pp_p <- inla.spde.make.A(mesh = mesh, loc = coord_pred) 

# Create stacks

# `data` has two columns, one for each likelihood

stk_y_e <- inla.stack(tag = "est_y",
                      data = list(y = cbind(data_USA$mean, NA), e = rep(NA, n_pts)),
                      A = list(1, ymat),
                      effects = list(mu = rep(1, n_pts), s = indxs)) 

stk_y_p <- inla.stack(tag = "pred_y",
                      data = list(y = cbind(rep(NA, n_pts_pred), NA), e = rep(NA, n_pts_pred)),
                      A = list(1, A_pp_p),
                      effects = list(mu = rep(1, n_pts_pred), s = indxs)) 

stk_pp_e <- inla.stack(tag = "est_pp",
                       data = list(y = cbind(NA, y_pp), e = e_pp),
                       A = list(1, A_pp),
                       effects = list(alpha_pp = rep(1, n_vtx + n_pts), v = indxv))

stk_pp_p <- inla.stack(tag = "pred_pp",
                       data = list(y = cbind(NA, rep(NA, n_pts_pred)), e = rep(1, n_pts_pred)),
                       A = list(1, A_pp_p),
                       effects = list(alpha_pp = rep(1, n_pts_pred), v = indxv))


# Full stack
stk_full_pp_y <- inla.stack(stk_y_e, stk_y_p, stk_pp_e, stk_pp_p)

# Fitting the model

re_prior <- list(prior = "gaussian", param = c(0, 1))

formula_1 <- y ~ 0 + mu + alpha_pp + f(s, model = spde) + f(v, copy = "s", fixed = FALSE, hyper = list(beta = re_prior))

model_6_1 <- inla(formula = formula_1,
                  family  = c("gaussian", "poisson"), 
                  E = inla.stack.data(stk_full_pp_y)$e,
                  data = inla.stack.data(stk_full_pp_y), 
                  control.predictor = list(link = rep(c(1, 2), c((n_pts + n_pts_pred), (n_vtx + n_pts + n_pts_pred))),
                                           compute = TRUE,
                                           A = inla.stack.A(stk_full_pp_y)))

saveRDS(object = model_6_1, file = "models/model_6_1.RDS")

summary(model_6_1)

ggplot(data.frame(inla.smarginal(model_6_1$marginals.hyperpar$`Beta for v`))) +
  geom_line(aes(x, y)) +
  labs(x = "", y = "", title = "Posterior of the degree of preferentiality (gamma)") + 
  custom_theme

# Fitted values and prediction

## Response

idx_y  <- inla.stack.index(stk_full_pp_y, tag = "pred_y" )$data

pred_y_mm  <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_y,  "mean"]))
pred_y_ll  <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_y,  "0.025quant"]))
pred_y_uu  <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_y,  "0.975quant"]))

r_mm <- pred_y_mm$V3 %>% range()
r_ll <- pred_y_ll$V3 %>% range()
r_uu <- pred_y_uu$V3 %>% range()
r <- c(min(r_mm[1], r_ll[1], r_uu[1]), max(r_mm[2], r_ll[2], r_uu[2]))

pp_y_mm <- plot_pred_USA(fitted_values = pred_y_mm, USA = USA, r = r, tt = "Mean")
pp_y_ll <- plot_pred_USA(fitted_values = pred_y_ll, USA = USA, r = r, tt = "2.5th")
pp_y_uu <- plot_pred_USA(fitted_values = pred_y_uu, USA = USA, r = r, tt = "97.5th")

(pp_y_ll + pp_y_mm + pp_y_uu) + plot_layout(guides = "collect") & theme(legend.position = "right")

## Intensity process

idx_pp <- inla.stack.index(stk_full_pp_y, tag = "pred_pp")$data

pred_pp_mm <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_pp, "mean"]))
pred_pp_ll <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_pp, "0.025quant"]))
pred_pp_uu <- as.data.frame(cbind(coord_pred, model_6_1$summary.fitted.values[idx_pp, "0.975quant"]))

# Expected number of observations
sum(pred_pp_mm$V3 * (25 ** 2))

r_mm <- pred_pp_mm$V3 %>% range()
r_ll <- pred_pp_ll$V3 %>% range()
r_uu <- pred_pp_uu$V3 %>% range()
r <- c(min(r_mm[1], r_ll[1], r_uu[1]), max(r_mm[2], r_ll[2], r_uu[2]))

pp_pp_mm <- plot_pred_USA(fitted_values = pred_pp_mm, USA = USA, r = r, tt = "Mean",   should_round = FALSE)
pp_pp_ll <- plot_pred_USA(fitted_values = pred_pp_ll, USA = USA, r = r, tt = "2.5th",  should_round = FALSE)
pp_pp_uu <- plot_pred_USA(fitted_values = pred_pp_uu, USA = USA, r = r, tt = "97.5th", should_round = FALSE)

(pp_pp_ll + pp_pp_mm + pp_pp_uu) + plot_layout(guides = "collect") & theme(legend.position = "right")









