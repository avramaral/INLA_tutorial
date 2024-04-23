source("R/header.R")

data <- readRDS(file = "data/example_2/m_ohio.RDS")

g <- inla.read.graph(filename = "data/example_2/map.adj")

# Replicated effects only share the hyperparameters. This means that the values of the random effects in the different replicates can be different.

formula_1 <- Y ~ 0 + 1 + f(id_area, model = "besag", graph = g, replicate = id_time)

model_3_1 <- inla(formula = formula_1,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE))

saveRDS(object = model_3_1, file = "models/model_3_1.RDS")

summary(model_3_1)

data_3_1 <- bind_cols(data[, c("county", "year", "geometry")], fitted.values = model_3_1$summary.fitted.values$mean)

plot_rr_years(d = data_3_1, tt = "Replicates")

# Additive model

formula_2 <- Y ~ 0 + 1 + f(id_time, model = "rw1") + f(id_area, model = "besag", graph = g)

model_3_2 <- inla(formula = formula_2,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE))

saveRDS(object = model_3_2, file = "models/model_3_2.RDS")

summary(model_3_2)

data_3_2 <- bind_cols(data[, c("county", "year", "geometry")], fitted.values = model_3_2$summary.fitted.values$mean)

plot_rr_years(d = data_3_2, tt = "Additive model")

# Separable spatio-temporal models with grouping (Kronecker product models)
# REF: https://github.com/hrue/r-inla/blob/devel/internal-doc/group/group-models.pdf

# [INLA for Spatial Statistics: Grouped models](http://faculty.washington.edu/jonno/SISMIDmaterial/8-Groupedmodels.pdf)
# Within group correlation
# Between group correlation

formula_3 <- Y ~ 0 + 1 + f(id_area, model = "besag", graph = g, 
                           group = id_time, control.group = list(model = "rw1"))

model_3_3 <- inla(formula = formula_3,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE))

saveRDS(object = model_3_3, file = "models/model_3_3.RDS")

summary(model_3_3)

data_3_3 <- bind_cols(data[, c("county", "year", "geometry")], fitted.values = model_3_3$summary.fitted.values$mean)
plot_rr_years(d = data_3_3, tt = "Kronecker product model (grouped by \"time\")")

##########

# Bernardinelli model [Bayesian analysis of space-time variation in disease risk](https://onlinelibrary.wiley.com/doi/10.1002/sim.4780142112)

data$id_area_cp <- data$id_area

formula_4 <- Y ~ 0 + 1 + f(id_area, model = "bym", graph = g) + f(id_area_cp, id_time, model = "iid") + id_time

model_3_4 <- inla(formula = formula_4,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE))

saveRDS(object = model_3_4, file = "models/model_3_4.RDS")

summary(model_3_4)

data_3_4 <- bind_cols(data[, c("county", "year", "geometry")], fitted.values = model_3_4$summary.fitted.values$mean)
plot_rr_years(d = data_3_4, tt = "Bernardinelli model")

