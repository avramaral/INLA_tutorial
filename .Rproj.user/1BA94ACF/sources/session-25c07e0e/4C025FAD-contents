## ----setting, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source('R/initial_setup.R')
opts_chunk$set(
  fig.path = 'figs/cond-sampl-'
)
library(splancs)

## ----halfpr0-------------------------------------------------------------
# loading the Parana border and precipitation data
data(PRborder)
data(PRprec)

# Which data points fall inside the left 50% of Parana
mid.long <- mean(range(PRborder[, 1]))
sel.loc <- which(PRprec[, 1] < mid.long)
sel.bor <- which(PRborder[, 1] < mid.long)

## ----halfpr, echo = FALSE, fig.height = 5, out.width = "80%", fig.cap = "Problem setting: available data over half of the domain."----
# visualizing it
par(mar =c(0, 0, 0, 0))
plot(PRborder, type = 'l', asp = 1, axes = FALSE)
points(PRprec[sel.loc, 1:2], cex = 0.7)
abline(v = mid.long, lty = 2)

## ----mesh1---------------------------------------------------------------
#  Building mesh1
mesh1 <- inla.mesh.2d(loc = PRprec[sel.loc, 1:2], max.edge = 1,
  cutoff = 0.1, offset = 1.2) 

## ----mesh2---------------------------------------------------------------
# define a boundary for the Parana border and an auxiliary mesh 
ibound <- inla.nonconvex.hull(PRborder, 0.05, 2, resol = 250) 
mesh2a <- inla.mesh.2d(mesh1$loc, boundary = ibound, 
  max.edge = 0.2, cutoff = 0.1) 
# Building mesh2 considering wider a boundary 
bound <- inla.nonconvex.hull(PRborder, 2) 
mesh2 <- inla.mesh.2d(loc = rbind(mesh1$loc, mesh2a$loc),
  boundary = bound, max.edge = 1, cutoff = 0.1) 

## ----plotmesh, echo=FALSE, fig.width=11, fig.height=5, fig.cap = '(ref:plotmesh)'----

par(mfrow = c(1, 2), mar = c(0, 0, 0, 0), xaxs = 'i', yaxs = 'i')
plot(mesh2$loc, type = 'n', asp = 1, main = "", axes = FALSE)
plot(mesh1, add = TRUE)
points(PRprec[sel.loc, 1:2], asp = 1, col = "blue", cex = 0.5)
lines(PRborder, col = 4)
plot(mesh2, asp = 1, main = "")
points(mesh1$loc, asp = 1,col = "red", cex = 0.5)
lines(PRborder, col = 4)

## ----spdes---------------------------------------------------------------
# SPDE model parameters 
range <- 3
std.u <- 1
# Define the SPDE models for mesh1 and mesh2
spde1 = inla.spde2.pcmatern(mesh1, prior.range = c(1, 0.1),
  prior.sigma = c(1, 0.1)) 
spde2 = inla.spde2.pcmatern(mesh2, prior.range = c(1, 0.1),
  prior.sigma = c(1, 0.1)) 

## ----Qmatrices-----------------------------------------------------------
# Obtain the precision matrix for spde1 and spde2
Q1 = inla.spde2.precision(spde1,
  theta = c(log(range), log(std.u))) 
Q2 = inla.spde2.precision(spde2,
  theta = c(log(range), log(std.u))) 

## ----usample-------------------------------------------------------------
u <- as.vector(inla.qsample(n = 1, Q = Q2, seed = 1))

## ----A-------------------------------------------------------------------
A1 <- inla.spde.make.A(mesh1,
  loc = as.matrix(PRprec[sel.loc, 1:2]))
A2 <- inla.spde.make.A(mesh2,
  loc = as.matrix(PRprec[sel.loc, 1:2]))

## ----ysample-------------------------------------------------------------
std.epsilon = 0.1
y <- drop(A2 %*% u) + rnorm(nrow(A2), sd = std.epsilon)

## ----datastack-----------------------------------------------------------
stk <- inla.stack(
  data = list(resp = y),
  A = list(A1, 1),
  effects = list(i = 1:spde1$n.spde,m = rep(1, length(y))),
  tag = 'est')

## ----fitting-------------------------------------------------------------
res <- inla(resp ~ 0 + m + f(i, model = spde1),
  data = inla.stack.data(stk),
  control.compute = list(config = TRUE),
  control.predictor = list(A = inla.stack.A(stk)))

## ----hyperpar------------------------------------------------------------
# Marginal for standard deviation of Gaussian likelihood
p.s.eps <- inla.tmarginal(function(x) 1 / sqrt(exp(x)), 
  res$internal.marginals.hyperpar[[1]])
# Summary of post. marg. of st. dev.
s.std <- unlist(inla.zmarginal(p.s.eps, silent = TRUE))[c(1:3, 7)] 

hy <- cbind(True = c(std.epsilon, range, std.u), 
 rbind(s.std, res$summary.hyperpar[2:3, c(1:3, 5)]))
rownames(hy) <- c('Std epsilon', 'Range field', 'Std field')
hy

## ----samples-------------------------------------------------------------
nn <- 100
s <- inla.posterior.sample(n = nn, res, intern = TRUE,
  seed = 1, add.names = FALSE)

## ----idxsamples----------------------------------------------------------
## Find the values of latent field "i" in samples from mesh1
contents <- res$misc$configs$contents
effect <- "i"
id.effect <- which(contents$tag == effect)
ind.effect <- contents$start[id.effect] - 1 + 
  (1:contents$length[id.effect])

## ----condsim-------------------------------------------------------------
# Obtain predictions at the nodes of mesh2
loc1 = mesh1$loc[,1:2]
loc2 = mesh2$loc[,1:2]
n = mesh2$n

mtch = match(data.frame(t(loc2)), data.frame(t(loc1)))
idx.c = which(!is.na(mtch))
idx.u = setdiff(1:mesh2$n, idx.c)
p = c(idx.u, idx.c)

ypred.mesh2 = matrix(c(NA), mesh2$n, nn)

m <- n - length(idx.c)
iperm <- numeric(m)

t0 <- Sys.time()
for(ind in 1:nn){
  
  Q.tmp = inla.spde2.precision(spde2, 
    theta = s[[ind]]$hyperpar[2:3])
  
  Q = Q.tmp[p, p]
  Q.AA = Q[1:m, 1:m]
  Q.BB = Q[(m + 1):n, (m + 1):n]
  Q.AB = t(Q[(m + 1):n, 1:m])
  Q.AA.sf = Cholesky(Q.AA,  perm = TRUE,  LDL = FALSE)
  perm = Q.AA.sf@perm + 1
  iperm[perm] = 1:m
  x = solve(Q.AA.sf, rnorm(m), system = "Lt")
  xc = s[[ind]]$latent[ind.effect]
  xx = solve(Q.AA.sf, -Q.AB %*% xc,  system = "A")
  
  x = rep(NA, n)
  x[idx.u] = c(as.matrix(xx))
  x[idx.c] = xc
  
  ypred.mesh2[, ind] = x
}
Sys.time() - t0

## ----projgrid------------------------------------------------------------
# Projection from the mesh nodes to a fine grid
projgrid  <- inla.mesh.projector(mesh2,
  xlim = range(PRborder[, 1]), ylim = range(PRborder[, 2]),
  dims = c(300, 300))

## ----visualize, out.width="70%", fig.height=10, fig.cap="The simulated field (top), the estimated posterior mean (middle) and the posterior marginal standard deviation (bottom)."----

# Find points inside the state of Parana
xy.in <- inout(projgrid$lattice$loc, PRborder[, 1:2])
# True field
r1 <- inla.mesh.project(projgrid , field = u)
r1[!xy.in] <- NA

# Mean predicted random field
r2 <- inla.mesh.project(projgrid , field = rowMeans(ypred.mesh2))
r2[!xy.in] <- NA

# sd of the predicted random field
sd.r2 <- inla.mesh.project(
  projgrid, field=apply(ypred.mesh2, 1, sd, na.rm = TRUE))
sd.r2[!xy.in] <- NA

# plotting
par(mfrow = c(3, 1), mar = c(0, 0, 0, 0))
zlm <- range(c(r1, r2), na.rm = TRUE)

# Map of the true field
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = r1),
  zlim = zlm)
points(PRprec[sel.loc, 1:2], col = "black", asp = 1, cex = 0.3)

## Map of the mean of the mean predicted random field
book.plot.field(list(x = projgrid$x, y = projgrid$y, z = r2),
  zlim = zlm)
points(PRprec[sel.loc, 1:2], col = "black", asp = 1, cex = 0.3)

book.plot.field(list(x = projgrid$x, y = projgrid$y, z = sd.r2))
points(PRprec[sel.loc, 1:2], col = "black", asp = 1, cex = 0.3)

