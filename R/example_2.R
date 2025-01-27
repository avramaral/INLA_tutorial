source("R/header.R")

# Load data
d_ohio <- read_csv(file = "data/example_2/data_ohio.csv")
m_ohio <- read_sf(dsn = "data/example_2/ohio_shapefile/", layer = "fe_2007_39_county")

# Process data
d <- d_ohio %>% group_by(NAME, year) %>% summarise(Y = sum(y)) %>% ungroup() %>% rename(county = NAME) %>% arrange(county, year)

## Expected cases 
d_ohio <- d_ohio %>% arrange(county, year, gender, race)
n.strata <- 4 # 2 genders x 2 races 
E <- expected(population = d_ohio$n, cases = d_ohio$y, n.strata = n.strata)

## Compute SIR
n_years    <- length(unique(d_ohio$year))
n_counties <- length(unique(d_ohio$NAME))

counties_E <- rep(unique(d_ohio$NAME), each = n_years)
years_E    <- rep(unique(d_ohio$year), times = n_counties)

d_E <- data.frame(county = counties_E, year = years_E, E = E)

d <- d %>% left_join(y = d_E, by = c("county", "year"))
d <- d %>% mutate(SIR = Y / E)

## Link map information
d_wide <- d %>% pivot_wider(names_from = year, values_from = c(Y, E, SIR), names_glue = "{.value}.{year}") # `wide` format
m_ohio <- m_ohio %>% rename(county = NAME) %>% left_join(y = d_wide, by = c("county")) %>% dplyr::select(c(colnames(d_wide), "geometry"))

s_cols <- setdiff(colnames(m_ohio), c("county", "geometry")) # to `long` format
m_ohio <- m_ohio %>% 
          pivot_longer(cols = all_of(s_cols), names_to = c(".value", "year"), names_sep = "\\.") %>% 
          dplyr::select(county, year, Y, E, SIR, geometry) %>%
          mutate(year = as.numeric(year)) %>% 
          arrange(county, year)

# Create indices for random effects
m_ohio <- m_ohio %>% mutate(id_area = as.numeric(as.factor(county)),
                            id_time = (1 + year - min(year)))

saveRDS(object = m_ohio, file = "data/example_2/m_ohio.RDS")

ggplot(m_ohio) + 
  geom_sf(aes(fill = SIR)) +
  facet_wrap(~ year, dir = "h", ncol = 4) +
  scale_fill_gradient2(name = "SIR", midpoint = 1, low = "blue", mid = "white", high = "red") + 
  labs(x = "", y = "") +
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks  = element_blank())
  

#############
# Inference #
#############

data <- m_ohio %>% filter(year == 1988)

formula_1 <- Y ~ 0 + 1

model_2_1 <- inla(formula = formula_1,
                  family  = "poisson", 
                  data = data, 
                  E = E, # Known component in the mean for the Poisson likelihoods
                  control.predictor = list(compute = TRUE),
                  control.compute   = list(cpo  = TRUE, 
                                           dic  = TRUE, 
                                           waic = TRUE)) # For model comparison

saveRDS(object = model_2_1, file = "models/model_2_1.RDS")

summary(model_2_1)

formula_2 <- Y ~ 0 + 1 + f(id_area, model = "iid")

model_2_2 <- inla(formula = formula_2,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE),
                  control.compute   = list(cpo  = TRUE, 
                                           dic  = TRUE, 
                                           waic = TRUE))

saveRDS(object = model_2_2, file = "models/model_2_2.RDS")

summary(model_2_2)

# Model comparison
table_model_comparison(models = list(model_2_1 = model_2_1, model_2_2 = model_2_2))

###############################
# Including spatial structure #
###############################

nb <- poly2nb(data$geometry)
head(nb)

# Save matrix for later use within INLA
nb2INLA("data/example_2/map.adj", nb)
g <- inla.read.graph(filename = "data/example_2/map.adj")

data$id_area_cp <- data$id_area

formula_3 <- Y ~ 0 + 1 + f(id_area, model = "iid") + f(id_area_cp, model = "besag", graph = g, scale.model = TRUE)

model_2_3 <- inla(formula = formula_3,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE),
                  control.compute   = list(cpo  = TRUE, 
                                           dic  = TRUE, 
                                           waic = TRUE))

saveRDS(object = model_2_3, file = "models/model_2_3.RDS")

summary(model_2_3)

table_model_comparison(models = list(model_2_1 = model_2_1, model_2_2 = model_2_2, model_2_3 = model_2_3))

formula_4 <- Y ~ 0 + 1 + f(id_area, model = "bym2", graph = g)

model_2_4 <- inla(formula = formula_4,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE),
                  control.compute   = list(cpo  = TRUE, 
                                           dic  = TRUE, 
                                           waic = TRUE))

saveRDS(object = model_2_4, file = "models/model_2_4.RDS")

summary(model_2_4)

table_model_comparison(models = list(model_2_1 = model_2_1, model_2_2 = model_2_2, model_2_3 = model_2_3, model_2_4 = model_2_4))

# Setting Penalized Complexity (PC) priors

pc_prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(1, 0.01)), # P(1 / sqrt(prec) > 1) = 0.01
  phi = list(
    prior = "pc",
    param = c(0.5, 0.75)) # P(phi < 0.5) = 0.75
)

formula_5 <- Y ~ 0 + 1 + f(id_area, model = "bym2", graph = g, hyper = pc_prior)

model_2_5 <- inla(formula = formula_5,
                  family  = "poisson", 
                  data = data, 
                  E = E,
                  control.predictor = list(compute = TRUE),
                  control.compute   = list(cpo  = TRUE, 
                                           dic  = TRUE, 
                                           waic = TRUE))

saveRDS(object = model_2_5, file = "models/model_2_5.RDS")

summary(model_2_5)

table_model_comparison(models = list(model_2_1 = model_2_1, model_2_2 = model_2_2, model_2_3 = model_2_3, model_2_4 = model_2_4, model_2_5 = model_2_5))

#################
# Fitted values #
#################

# Selected model
model_2_3

fitted_values <- model_2_3$summary.fitted.values[, c("mean", "0.025quant", "0.975quant")] %>% 
                 as_tibble() %>% 
                 mutate(id_area = 1:nrow(data)) %>% 
                 rename(Mean = mean, `0.025` = `0.025quant`, `0.975` = `0.975quant`) %>% 
                 left_join(y = data[, c("county", "geometry", "id_area")], by = "id_area") %>% 
                 dplyr::select(county, Mean, `0.025`, `0.975`, geometry) %>% 
                 pivot_longer(cols = c(Mean, `0.025`, `0.975`)) %>% 
                 mutate(name = factor(name, levels = c("0.025", "Mean", "0.975"))) %>% 
                 st_as_sf()

ggplot(fitted_values) + 
  geom_sf(aes(fill = value)) +
  facet_wrap(~ name, dir = "h", ncol = 3) +
  scale_fill_gradient2(name = "Relative risk", midpoint = 1, low = "blue", mid = "white", high = "red") + 
  labs(x = "", y = "", title = "Estimated \"Relative Risk\" in 1988") +
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks  = element_blank())


# Random effects

random_effects <- data[, c("county", "geometry")] %>% 
                  mutate("IID" = model_2_3$summary.random$id_area$mean,
                         "BESAG" = model_2_3$summary.random$id_area_cp$mean) %>% 
                  pivot_longer(cols = c(IID, BESAG)) %>% 
                  mutate(name = factor(name, levels = c("IID", "BESAG"))) %>% 
                  st_as_sf()

ggplot(random_effects) + 
  geom_sf(aes(fill = value)) +
  facet_wrap(~ name, dir = "h", ncol = 2) +
  scale_fill_gradient2(name = "RE", midpoint = 0, low = "blue", mid = "white", high = "red") + 
  labs(x = "", y = "", title = "Posterior mean of the random effects") +
  custom_theme + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks  = element_blank())
