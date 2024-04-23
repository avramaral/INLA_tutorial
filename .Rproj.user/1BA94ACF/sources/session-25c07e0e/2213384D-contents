source("R/header.R")

data_USA <- readRDS(file = "data/example_4/data_USA.rds")
USA <- readRDS(file = "data/example_4/USA_filtered.rds")

mesh <- readRDS(file = "data/example_4/mesh.RDS")
dual_mesh <- readRDS(file = "data/example_5/dual_mesh.RDS")
wgt <- readRDS(file = "data/example_5/wgt.RDS")

coord_pred <- readRDS(file = "data/example_4/coord_pred.RDS")

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

# Basis function

bandwidth <- 2500
smoothing_kernel <- "Wendland"
# center_pts <- as.matrix(rbind(c(-1600, 5307), c(  300, 3800), c( 1114, 4229),
#                               c( -800, 4300), c(  263, 4987), c(-2000, 4500)))

center_pts <- as.matrix(rbind(c(-2372, 4346), c(-1984, 5363), c( -918, 5143),
                              c( -518, 3817), c(  263, 4987), c(  677, 4096),
                              c( 1494, 4720), c( -198, 4291), c(-1710, 4513),
                              c(  875, 3602), c(-1202, 3742), c(-1310, 4946),
                              c( -190, 5189), c(  920, 4550), c( -233, 4800)))

bfs      <- basis_functions(center_pts = center_pts, loct = data_coor,  mesh = mesh, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth)
bfs_pred <- basis_functions(center_pts = center_pts, loct = coord_pred, mesh = NULL, smoothing_kernel = smoothing_kernel, bandwidth = bandwidth)

pps <- list()
for (i in 1:nrow(center_pts)) {
  pps[[i]] <- plot_pred_USA(fitted_values = as.data.frame(cbind(coord_pred, bfs_pred[[i]])), USA = USA, r = c(0, 1), tt = paste("Wendland ", i, sep = ""))
}

((pps[[ 1]] + pps[[ 2]] + pps[[ 3]]) / 
 (pps[[ 4]] + pps[[ 5]] + pps[[ 6]]) /
 (pps[[ 7]] + pps[[ 8]] + pps[[ 9]]) / 
 (pps[[10]] + pps[[11]] + pps[[12]]) /
 (pps[[13]] + pps[[14]] + pps[[15]])) + plot_layout(guides = "collect") & theme(legend.position = "right")

# Data and projection matrices

n_vtx <- mesh$n
n_pts <- nrow(data_USA)
n_pts_pred <- nrow(coord_pred)

indx <- 1:spde$n.spde

y_pp <- rep(0:1, c(n_vtx, n_pts))
e_pp <- c(wgt, rep(0, n_pts)) 
imat <- Diagonal(n_vtx, rep(1, n_vtx))
ymat <- inla.spde.make.A(mesh, data_coor)

# Estimation
A_pp_base <- rbind(imat, ymat)
n_basis_functions <- nrow(center_pts)
A_pp <- list()
for (i in 1:n_basis_functions) {
  A_pp[[i]] <- A_pp_base * bfs[[i]] # Multiply each column by the basis function evaluated at the corresponding locations
}

# Prediction
A_pp_p_base <- inla.spde.make.A(mesh = mesh, loc = as.matrix(coord_pred))
A_pp_p <- list()
for (i in 1:n_basis_functions) {
  A_pp_p[[i]] <- A_pp_p_base * bfs_pred[[i]]
}

# Create stacks

#########################
# Create a list of effects based on the number of basis function
#########################

fx_effects <- list(rep(1, n_vtx + n_pts))
rd_effects <- list()
for (i in 1:length(A_pp)) { rd_effects <- append(rd_effects, list(indx)) } 
effects <- append(fx_effects, rd_effects)
names(effects) <- c("alpha_pp", paste("v_", 1:length(A_pp), sep = ""))

effects_pred <- effects
effects_pred[[1]] <- rep(1, times = n_pts_pred)

#########################
#########################

# Stacks for `y` are similar to before; however, stacks for `pp` have different structures for `A` and `effects`

stk_y_e <- inla.stack(tag = "est_y",
                      data = list(y = cbind(data_USA$mean, NA), e = rep(NA, n_pts)),
                      A = list(1, ymat),
                      effects = list(mu = rep(1, n_pts), s = indx)) 

stk_y_p <- inla.stack(tag = "pred_y",
                      data = list(y = cbind(rep(NA, n_pts_pred), NA), e = rep(NA, n_pts_pred)),
                      A = list(1, A_pp_p_base),
                      effects = list(mu = rep(1, n_pts_pred), s = indx)) 

stk_pp_e <- inla.stack(tag = "est_pp",
                       data = list(y = cbind(NA, y_pp), e = e_pp),
                       A = append(1, A_pp),
                       effects = effects)

stk_pp_p <- inla.stack(tag = "pred_pp",
                       data = list(y = cbind(NA, rep(NA, n_pts_pred)), e = rep(1, n_pts_pred)),
                       A = append(1, A_pp_p),
                       effects = effects_pred)

# Full stack
stk_full_pp_y <- inla.stack(stk_y_e, stk_y_p, stk_pp_e, stk_pp_p)

re_prior <- list(prior = "gaussian", param = c(0, 1))

random_effects <- paste(paste("f(v_",  1:length(A_pp), ", copy = \"s\", fixed = FALSE, hyper = list(beta = re_prior))", sep = ""), collapse = " + ")
formula_1 <- paste("y ~ 0 + mu + alpha_pp + f(s, model = spde) + ", random_effects, sep = "")
formula_1 <- eval(parse(text = formula_1))
formula_1

model_7_1 <- inla(formula = formula_1,
                  family  = c("gaussian", "poisson"), 
                  E = inla.stack.data(stk_full_pp_y)$e,
                  data = inla.stack.data(stk_full_pp_y), 
                  control.predictor = list(link = rep(c(1, 2), c((n_pts + n_pts_pred), (n_vtx + n_pts + n_pts_pred))),
                                           compute = TRUE,
                                           A = inla.stack.A(stk_full_pp_y)))

saveRDS(object = model_7_1, file = "models/model_7_1.RDS")

summary(model_7_1)

# Fitted values and prediction

## Response

idx_y  <- inla.stack.index(stk_full_pp_y, tag = "pred_y" )$data

pred_y_mm  <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_y,  "mean"]))
pred_y_ll  <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_y,  "0.025quant"]))
pred_y_uu  <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_y,  "0.975quant"]))

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

pred_pp_mm <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_pp, "mean"]))
pred_pp_ll <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_pp, "0.025quant"]))
pred_pp_uu <- as.data.frame(cbind(coord_pred, model_7_1$summary.fitted.values[idx_pp, "0.975quant"]))

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

## Preferentiality surface


n_basis_functions
n_hyperparameter <- nrow(model_7_1$summary.hyperpar)

coeff <- model_7_1$summary.hyperpar[((n_hyperparameter - n_basis_functions + 1):n_hyperparameter), "mean"]
value <- rep(0, n_pts_pred) 
for (i in 1:n_basis_functions) {
  bfs_pred_tmp <- bfs_pred[[i]]
  value <- value + (bfs_pred_tmp * coeff[i])
}

pref_mm <- as.data.frame(cbind(coord_pred, value))

plot_pred_USA(fitted_values = pref_mm, USA = USA, r = range(pref_mm$value), tt = "Preferentiality surface", should_round = FALSE)


compute_preferentiatility <- function (fit, coord_pred, bfs_pred, xlim, by, n_coord_pr = NA, ...) {
  
  if (is.na(n_coord_pr)) { n_coord_pr <- length(bfs_pred[[1]]) }
  
  n_base_functions <- length(bfs_pred)
  n_hyperparameter <- nrow(fit$summary.hyperpar)
  
  coeff <- fit$summary.hyperpar[((n_hyperparameter - n_base_functions + 1):n_hyperparameter), "mean"][1:n_coord_pr]
  value <- rep(x = 0, times = n_coord_pr) 
  for (i in 1:n_base_functions) {
    # bfs_pred_tmp <- c(matrix(data = bfs_pred[[i]], nrow = length(seq(from = xlim[1] + (by / 2), to = xlim[2] - (by / 2), by = by)), byrow = TRUE))
    bfs_pred_tmp <- bfs_pred[[i]][1:n_coord_pr]
    value <- value + (bfs_pred_tmp * coeff[i])
  }
  
  pref <- data.frame(x = coord_pred$x[1:n_coord_pr], y = coord_pred$y[1:n_coord_pr], z = value)
  
  # Create a gridded spatial object from "pref"
  coordinates(pref) <- ~ x + y
  gridded(pref) <- TRUE
  
  raster(pref)
}


  fitted_preferentiality <- compute_preferentiatility(fit = result$fit, coord_pred = coord_pred, bfs_pred = bfs_pred, xlim = xlim, by = by, n_coord_pr = n_coord_pr)
  crs(fitted_preferentiality) = "+init=epsg:6345 +units=km +no_defs"
  
  fitted_preferentiality    <- as(fitted_preferentiality, "SpatialPixelsDataFrame")
  fitted_preferentiality_df <- as.data.frame(fitted_preferentiality)
  colnames(fitted_preferentiality_df) <- c("estimated", "x", "y")
  
  pal <- jet.col(n = 100, alpha = 0.9)
  labs <- seq(round(min(fitted_preferentiality_df$estimated) - 0.001, 3), round(max(fitted_preferentiality_df$estimated) + 0.001, 3), length.out = 6)
  
  p_3 <- ggplot() +
    geom_tile(data = fitted_preferentiality_df, mapping = aes(x = x, y = y, fill = estimated)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = "Deg. Pref.", colors = pal, breaks = labs, labels = as.character(format(labs, nsmall = 3)), limits = c(labs[1], tail(labs, 1))) +
    labs(x = "Longitude", y = "Latitude", title = "") + 
    theme_bw() +
    theme(text = element_text(size = 24, family = "LM Roman 10"), 
          legend.key.height = unit(2.57, "cm"),
          legend.key.width  = unit(0.75, "cm"),
          plot.margin = margin(0, 0, 0, 0, "cm"))








