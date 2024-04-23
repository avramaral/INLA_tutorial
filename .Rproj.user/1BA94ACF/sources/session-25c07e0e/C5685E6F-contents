source("R/header.R")

data_USA <- readRDS(file = "data/example_4/data_USA.rds")
USA <- readRDS(file = "data/example_4/USA_filtered.rds")

ggplot() + 
  geom_sf(data = USA, fill = "white") +
  geom_sf(data = data_USA, color = "black", size = 3, shape = 3) +
  labs(x = "", y = "", title = "") +
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks  = element_blank())

mesh <- readRDS(file = "data/example_4/mesh.RDS")

# Priors might differ
spde <- inla.spde2.pcmatern(mesh = mesh, 
                            alpha = 2,
                            prior.range = c(1e3, 0.90), # P(range < 1e3) = 0.90
                            prior.sigma = c(1.0, 0.01)) # P(sigma > 1.0) = 0.01


dual_mesh <- inla.dual.mesh(mesh)

saveRDS(object = dual_mesh, file = "data/example_5/dual_mesh.RDS")

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

USA_poly <- SpatialPolygons(list(Polygons(list(Polygon(USA_coor)), ID = "1")))
wgt <- sapply(1:length(dual_mesh), function (i) {
  if (gIntersects(dual_mesh[i, ], USA_poly)) {
    return(gArea(gIntersection(dual_mesh[i, ], USA_poly)))
  } else {
    return(0)
  }
})
saveRDS(object = wgt, file = "data/example_5/wgt.RDS")
sum(wgt) # Area of the study region

{
  plot(mesh$loc, asp = 1, col = (wgt == 0) + 1, pch = 19, xlab = "", ylab = "", axes = F) 
  plot(dual_mesh, add = TRUE)
  plot(USA$geometry, add = TRUE, border = "green", lwd = 2)
}


# Data and projection matrices

n_vtx <- mesh$n
n_pts <- nrow(data_USA)

# Indices
indxs <- inla.spde.make.index("s", spde$n.spde)

# augmented data: `0` for the mesh nodes and `1` for the observations
y_pp <- rep(0:1, c(n_vtx, n_pts))
# Exposure vector
e_pp <- c(wgt, rep(0, n_pts)) 
# Projection matrix (in two parts)
# (1) For the integration points, this is just a diagonal matrix—as these locations are just the mesh vertices
imat <- Diagonal(n_vtx, rep(1, n_vtx))
# (2) For the observed points, the projection matrix is defined with `inla.spde.make.A`
ymat <- inla.spde.make.A(mesh, data_coor)
# (1) + (2)
A_pp <- rbind(imat, ymat)

coord_pred <- readRDS(file = "data/example_4/coord_pred.RDS")
n_pts_pred <- nrow(coord_pred)

# Prediction
A_pp_p <- inla.spde.make.A(mesh = mesh, loc = coord_pred) # Projection matrix for the prediction points

# Create stacks
# Stack for estimation 
stk_pp_e <- inla.stack(tag = "est_pp",
                       data = list(y = y_pp, e = e_pp),
                       A = list(1, A_pp),
                       effects = list(alpha_pp = rep(1, n_vtx + n_pts), s = indxs))

stk_pp_p <- inla.stack(tag = "pred_pp",
                       data = list(y = rep(NA, n_pts_pred), e = rep(1, n_pts_pred)),
                       A = list(1, A_pp_p),
                       effects = list(alpha_pp = rep(1, n_pts_pred), s = indxs))

# Full stack
stk_full_pp <- inla.stack(stk_pp_e, stk_pp_p)

# Fitting the model

formula_1 <- y ~ 0 + alpha_pp + f(s, model = spde)

model_5_1 <- inla(formula = formula_1,
                  family  = "poisson", 
                  E = inla.stack.data(stk_full_pp)$e,
                  data = inla.stack.data(stk_full_pp), 
                  control.predictor = list(link = 1,
                                           compute = TRUE,
                                           A = inla.stack.A(stk_full_pp)))

saveRDS(object = model_5_1, file = "models/model_5_1.RDS")

summary(model_5_1)

# Fitted values and prediction

idx_pp <- inla.stack.index(stk_full_pp, tag = "pred_pp")$data

pred_pp_mm <- as.data.frame(cbind(coord_pred, model_5_1$summary.fitted.values[idx_pp, "mean"]))
pred_pp_ll <- as.data.frame(cbind(coord_pred, model_5_1$summary.fitted.values[idx_pp, "0.025quant"]))
pred_pp_uu <- as.data.frame(cbind(coord_pred, model_5_1$summary.fitted.values[idx_pp, "0.975quant"]))

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


