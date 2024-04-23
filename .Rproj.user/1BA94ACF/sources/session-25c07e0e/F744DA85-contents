library("INLA")
library("tidyverse")
library("spatstat")
library("raster")
library("sf")
library("rgeos")
library("MASS")
library("spdep")
library("patchwork")

pal <- c("#00008FFF", "#0000F2FF", "#0063FFFF", "#00D4FFFF", "#46FFB8FF", "#B8FF46FF", "#FFD400FF", "#FF6300FF", "#F00000FF", "#800000FF")

custom_theme <-  theme_bw() + theme(legend.position = "right", 
                                    text = element_text(size = 14),
                                    plot.title = element_text(size = 16),
                                    legend.title = element_text(size = 12))

##################################################
##################################################

# Extracted from `SpatialEpi`
expected <- function (population, cases, n.strata, ...) {
  
  n <- length(population) / n.strata
  E <- rep(0, n)
  qNum <- rep(0, n.strata)
  qDenom <- rep(0, n.strata)
  q <- rep(0, n.strata)
  
  # Compute strata-specific rates
  for (i in 1:n.strata) {
    indices <- rep(i, n) + seq(0, (n - 1)) * n.strata
    qNum[i] <- qNum[i] + sum(cases[indices])
    qDenom[i] <- qDenom[i] + sum(population[indices])
  }
  q <- qNum / qDenom
  
  # Compute expected counts
  for (i in 1:n) {
    indices <- (1:n.strata) + (i - 1) * n.strata
    E[i] <- sum(population[indices] * q)
  }
  
  E
}

table_model_comparison <- function (models, ...) {
  
  n_models <- length(models)
  df <- as.data.frame(matrix(data = 0, nrow = n_models, ncol = 3))
  rownames(df) <- names(models)
  colnames(df) <- c("CPO", "DIC", "WAIC")
  
  for (i in 1:n_models) {
    tmp_CPO  <- sum(log(models[[i]]$cpo$cpo)) * -1
    tmp_DIC  <- models[[i]]$dic$dic
    tmp_WAIC <- models[[i]]$waic$waic
    
    df[i, ] <- c(tmp_CPO, tmp_DIC, tmp_WAIC)
  }
  
  df
}

plot_rr_years <- function (d, col_name = "fitted.values", temporal_name = "year", tt = "", ...) {
  ggplot(d) + 
    geom_sf(aes(fill = .data[[col_name]])) +
    facet_wrap(~ .data[[temporal_name]], dir = "h", ncol = 7) +
    scale_fill_gradient2(name = "Relative risk", midpoint = 1, low = "blue", mid = "white", high = "red") + 
    labs(x = "", y = "", title = tt) +
    custom_theme + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank())
}

create_prediction_grid <- function (loc, resolution = 25, ...) {
  pts_bdy <- loc$geometry[[1]][[1]]
  pts_bdy_x <- range(pts_bdy[, 1])
  pts_bdy_y <- range(pts_bdy[, 2])
  coord_pred <- expand.grid(x = seq(pts_bdy_x[1], pts_bdy_x[2], by = resolution), y = seq(pts_bdy_y[1], pts_bdy_y[2], by = resolution))
  coordinates(coord_pred) <- ~ x + y
  
  xx <- as(st_as_sf(coord_pred), "sf");   st_crs(xx) <- st_crs(loc$geometry)
  yy <- as(st_as_sf(loc$geometry), "sf"); st_crs(yy) <- st_crs(loc$geometry)
  
  # Compute intersection between grid and `loc` borders
  pppts <- st_intersection(x = xx, y = yy)
  
  coord_pred <- matrix(data = NA, nrow = length(pppts$geometry), ncol = 2)
  colnames(coord_pred) <- c("x", "y")
  for (p in 1:length(pppts$geometry)) { coord_pred[p, ] <- sf::st_coordinates(pppts$geometry[[p]]) }
  coord_pred <- data.frame(x = coord_pred[, 1], y = coord_pred[, 2])
  
  as.matrix(coord_pred)
}

plot_pred_USA <- function (fitted_values, USA, r, tt = "", should_round = TRUE, ...) {
  
  coordinates(fitted_values) <- ~ x + y
  gridded(fitted_values) <- TRUE
  fitted_values <- raster(fitted_values)
  crs(fitted_values) <- "+init=epsg:6345 +units=km +no_defs"
  
  fitted_values    <- as(fitted_values, "SpatialPixelsDataFrame")
  fitted_values_df <- as.data.frame(fitted_values)
  colnames(fitted_values_df) <- c("pred", "x", "y")
  
  if (should_round) {
    breaks <- seq(floor(r[1]), ceiling(r[2]), length.out = 5)
  } else {
    breaks <- seq(r[1], r[2], length.out = 5)
  }
  
  pal <- c("#00008FFF", "#0000F2FF", "#0063FFFF", "#00D4FFFF", "#46FFB8FF", "#B8FF46FF", "#FFD400FF", "#FF6300FF", "#F00000FF", "#800000FF")
  
  pp <- ggplot() +
    geom_tile(data = fitted_values_df, mapping = aes(x = x, y = y, fill = pred)) + 
    geom_sf(data = USA, color = "black", fill = NA, lwd = 0.5) +
    scale_fill_gradientn(name = "PM2.5", colors = pal, breaks = breaks, limits = c(breaks[1], tail(breaks, 1))) + 
    labs(x = "", y = "", title = tt) + 
    custom_theme + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank())
  
  pp
} 

# Extracted from "Advanced Spatial Modeling with Stochastic Partial Differential Equations Using R and INLA (Krainski et al., 2019)"
inla.dual.mesh <- function (mesh, ...) {
  if (mesh$manifold == "R2") {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function (i) { colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2]) }))
    library("parallel")
    pls <- mclapply(1:mesh$n, function (i) {
      p <- unique(Reduce("rbind", lapply(1:3, function (k) {
        j <- which(mesh$graph$tv[,k] == i)
        if (length(j) > 0) {
          return(rbind(ce[j, , drop = FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] + mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] + mesh$loc[mesh$graph$tv[j, c(2:3, 1)[k]], 2]) / 2))
        } else {
          return(ce[j, , drop = FALSE])
        }
      })))
      j1 <- which(mesh$segm$bnd$idx[, 1] == i)
      j2 <- which(mesh$segm$bnd$idx[, 2] == i)
      if ((length(j1) > 0) | (length(j2) > 0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2] / 2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2] / 2 + mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2] / 2))
        yy <- p[, 2] - mean(p[, 2]) / 2 - mesh$loc[i, 2] / 2
        xx <- p[, 1] - mean(p[, 1]) / 2 - mesh$loc[i, 1] / 2
      } else {
        yy <- p[, 2] - mesh$loc[i, 2]
        xx <- p[, 1] - mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy, xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function (i) { Polygons(list(pls[[i]]), i)} )))
  }
  else { stop("It only works for R2.") }
}

basis_functions <- function (center_pts, loct, mesh = NULL, bandwidth, smoothing_kernel = "Wendland", ...) {
  if (is.null(mesh)) {
    all_pts <- loct
  } else {
    all_pts <- rbind(mesh$loc[, 1:2], loct)
  }
  
  f <- get(smoothing_kernel)
  b <- list()
  n_basis_functions <- nrow(center_pts)
  
  for (i in 1:n_basis_functions) { b[[i]] <- f(x = all_pts[, 1], y = all_pts[, 2], center_x = center_pts[i, 1], center_y = center_pts[i, 2], bandwidth = bandwidth) }
  
  b
}

Wendland <- function (x, y, center_x, center_y, bandwidth = 0.5, ...) {
  pts <- cbind((x - center_x), (y - center_y))
  d <- apply(X = pts, MARGIN = 1, FUN = function (x) { sqrt(sum(x ** 2)) }) / bandwidth
  k <- (1 - d) ^ 6 * (35 * d ^ 2 + 18 * d + 3) / 3 * as.numeric(d <= 1)
  k <- k / max(k)
  k
}

