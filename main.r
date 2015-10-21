## Hierarchichal spatial model for the presence of fish with INLA
##
## This is companion code of 
## F. Muñoz, M. G. Pennino, D. Conesa, A. López-Quílez, J. M. Bellido (2013).
## Estimation and prediction of the spatial occurrence of fish species using
## Bayesian latent Gaussian models. Stochastich Environmental Research and Risk
## Assessment 27(5):1171-1180. +init=epsg:23030DOI:10.1007/s00477-012-0652-3
##
## NOTE: The code has been revised to meet current practices, methods and
## syntaxis in INLA. Therefore, the results can differ slightly from those in
## the paper.
## 
## LICENCSE: GPL-3 (https://www.gnu.org/copyleft/gpl.html)
## (c) 2015 Facundo Muñoz


library(INLA)       # Tested with INLA version 0.0-1441315467 (3 Sep 2015)
library(parallel)   # Parallel computation
library(doParallel) # ...
library(foreach)    # ...
library(knitr)      # kable()
library(lattice)    # Plots of maps

## Auxiliar functions 
source("functions.r")


## Data and cartography
## http://dx.doi.org/10.1007/s00477-012-0652-3
obs <- read.csv("https://zenodo.org/record/32494/files/MediterraneusT.csv", 
                 header = TRUE)
load(url('https://zenodo.org/record/32494/files/spatial_data.RData'))


## Derive the log-Bathymetry as an alternative covariate
logDepth <- Depth
slot(logDepth, 'data')[[1]] <- 
  log(1+abs(slot(logDepth, 'data'))[[1]])


## Covariates at observation spots
obs.sp <- SpatialPointsDataFrame(
  coords  = cbind(obs$x, obs$y),
  data    = obs,
  proj4string = CRS(proj4string(Depth))
)

dat <- transform(
  obs,
  Depth    = over(obs.sp, Depth)[[1]],     # Bathymetry
  logDepth = over(obs.sp, logDepth)[[1]],  # log-Bathymetry
  Chl     = over(obs.sp, Chl)[[1]],        # Chlorophyll a
  U       = seq_len(nrow(obs))             # Heterogeneity effect (nugget)
)


## Build the list of competing GLM models:
## List all combinations of potential variables and then
## remove those models with both Depth and logDepth
sel.terms <- c('Depth', 
               'logDepth',
               'Chl',
               'f(spatial, model=spde)'
               # , 'f(Year, model = "iid")'   # Numerical issues: excluded
               )
f.list <- unlist(sapply(0:length(sel.terms), comb.terms))
f.list = f.list[-grep('Depth \\+ logDepth', f.list)]


## Prediction mesh
## set the maximum length for the triangulation
## as 10% of the diameter of the locations
max.edge <- max(dist(coordinates(obs)))*.1
cutoff <- 500   # Merge close vertices into one
# Triangulate the region, without restrictions (mesh/default)
mesh <- inla.mesh.create(
  coordinates(obs.sp),
  cutoff = cutoff,
  extend=list(n=8, offset=-0.1),  # Encapsulate the data region with a relative margin
  refine=list(min.angle=26,       # Refined triagulation
              max.edge=max.edge)
)


## Prediction covariates
## Values of variables in triangulation vertices
pred.loc <- SpatialPoints(mesh$loc[,-3],
                          proj4string = CRS(proj4string(Depth)))
covar.pred <- data.frame(
  mu = 1,
  spatial = seq_len(mesh$n),
  sapply(list(Depth = Depth, logDepth = logDepth, Chl = Chl),
         function(x)
           over(pred.loc, x)[[1]])
)

## Observation and Prediction stacks
A.obs  <- inla.spde.make.A(mesh, loc = coordinates(obs.sp))
stack.obs <- 
  inla.stack(data    = list(Y = dat$Presence),
             A       = list(A.obs, 1),
             effects = list(data.frame(mu      = rep(1, mesh$n),
                                       spatial = seq_len(mesh$n)),
                            dat),
             tag="obs")
stack.pred <- 
  inla.stack(data    = list(Y = NA),
             A       = 1,
             effects = covar.pred,
             tag     = "pred")
stack  <-  inla.stack(stack.obs, stack.pred)






## Prior specifications
## Default vague priors for the fixed effects: N(0, 1e-05)
# plot(mesh)
size = max(apply(mesh$loc[,1:2], 2, function(x) diff(range(x))))
sigma0 = 1      # prior median sd
range0 = size/5 # prior median range (~ 20% diameter of region)
kappa0 = sqrt(8)/range0
tau0   = 1/(sqrt(4*pi)*kappa0*sigma0)
## Code for tuning the prior range
# ## the range can be between about size/100 and size, so:
# sqrt(8)/size*c(1, 100)  # kappa0 \in 4.4e(-05, -03), and
# log(sqrt(8)/size*c(1, 100))  # theta2 \in (-10, -5.4), with median
# log(kappa0)   # -8.4
# ## This gives a semi-support of:
# max(abs(log(sqrt(8)/size*c(1, 100)) - log(kappa0) )) # 3 units in the log scale
# ## so we use a gaussian prior for theta2 with mean log(kappa0)
# ## and sd = 3/3 (as 3 sd is the effective support of a Normal dist.)
# ## Being generous, we use an sd = 1.5
rho0.mar <- 
  inla.tmarginal(function(x) sqrt(8)/exp(x),
                 transform(data.frame(x = log(kappa0) + seq(-2.5, 10,
                                                            length=2021)),
                           y = dnorm(x, mean = log(kappa0),
                                     sd = 1.5)))
# plot(rho0.mar, type = 'l', xlab = 'distance', ylab = '')

sigma20.mar <- 
  inla.tmarginal(function(x) 1/(4*pi*kappa0^2*exp(2*x)),
                 transform(data.frame(x = log(tau0) + seq(-1.7, 10,
                                                          length=2021)),
                           y = dnorm(x, mean = log(tau0),
                                     sd = 1)))
# plot(sigma20.mar, type = 'l', xlab = 'variance', ylab = '')

sigma_gt1 <- function(x, kappa) 1/sqrt(4*pi)/kappa/exp(x)
# sigma_gt1(log(tau0)+3, kappa0)

spde = inla.spde2.matern(mesh,
                         B.tau = cbind(log(tau0), -1, 1),
                         B.kappa = cbind(log(kappa0), 0, -1),
                         theta.prior.mean = c(0, 0),
                         theta.prior.prec = c(1, 1/1.5**2),
                         constr = TRUE)


## Parallel fit of models
cl <- makeCluster(8)   # number of cores available
registerDoParallel(cl)
# Go for it!
rr <- foreach(i=seq_along(f.list)) %dopar% 
  runinla(formula = eval(parse(text=f.list[i])), 
          stack   = stack, 
          mesh    = mesh,
          spde    = spde,
          cf      = list())

stopCluster(cl)


## Results

## Selected model
r <- rr[[12]]

## Region limits (extended by 25%)
## (make sure we have covariates in the full prediction area)
xl <- bbox(obs.sp)[1,] + .25*c(-1,1)*diff(bbox(obs.sp)[1,])  # ~ 80km
yl <- bbox(obs.sp)[2,] + .25*c(-1,1)*diff(bbox(obs.sp)[2,])  # ~ 33km

## Mapping between triangulation vertices and grid points:
proj  <-  inla.mesh.projector(mesh, 
                              dims=c(200, 200))

## Stack indices
obs.idx <- inla.stack.index(stack, 'obs')
pred.idx <- inla.stack.index(stack, 'pred')

## Coordinates of the coast for plotting purposes
coord_coast <- coordinates(slot(slot(coast, "polygons")[[1]], 
                                "Polygons")[[1]])


## Fig. 1 Visualize data and region
par(mai=rep(0.1,4))
plot.new()
plot.window(xlim = range(mesh$loc[, 1]), ylim = range(mesh$loc[, 2]), "", asp=1)
plot(mesh, col='lightgray', add=T, xlab='Easting')
plot(coast, col='lightgray', add=T)
plot(obs.sp, pch=ifelse(obs$Presence, 20, 1), col=as.numeric(obs$Presence)+1, add=T)
SpatialPolygonsRescale(layout.north.arrow(1), offset= c(573849,4072084), scale = 6000,
                       plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(578667,4051436), scale= 10000, fill=c("transparent", "black"), plot.grid= F)
# Annotations
#    z <- locator()
text(583667,4050436, "10KM", cex= 1)
text(553203,4046000, "Mediterranean Sea", cex= 1)
z <- list(x=549564.8, y=4077379)
points(z$x, z$y, cex=3, col='black', pch=19)
text(z$x, z$y, pos=3, offset=.8, "Almería")
box()


## Fig. 2 Maps of covariates
spplot(Depth, col.regions=topo.colors, sp.layout=coast, xlim=xl, ylim=yl)
spplot(Chl, col.regions=topo.colors, sp.layout=coast, xlim=xl, ylim=yl)



## Table 1 Model comparison 
f.list.labels <- gsub('Y ~ ', '', f.list)
f.list.labels <- gsub('f\\(spatial, model=spde\\)', '\\$\\\\theta\\$', f.list.labels)
f.list.labels <- gsub('f\\(Year, model=\\"iid\\"\\)', '\\$Y\\$', f.list.labels)
kable(
  data.frame(
    Model = f.list.labels, 
    LCPO = sapply(1:length(rr), function(x) -mean(log(rr[[x]]$cpo$cpo[obs.idx$data]))),
    BS = sapply(1:length(rr), function(x) mean((1-rr[[x]]$cpo$cpo[obs.idx$data])^2)),
    DIC = unlist(sapply(1:length(rr), function(x) ifelse(is.null(rr[[x]]$dic), NA, rr[[x]]$dic$dic))),
    P.eff = unlist(sapply(1:length(rr), function(x) ifelse(is.null(rr[[x]]$dic), NA, rr[[x]]$dic$p.eff))),
    Mean.dev = unlist(sapply(1:length(rr), function(x) ifelse(is.null(rr[[x]]$dic), NA, rr[[x]]$dic$mean.deviance))),
    Dev.mean = unlist(sapply(1:length(rr), function(x) ifelse(is.null(rr[[x]]$dic), NA, rr[[x]]$dic$deviance.mean)))
  ), 
  digits = 2
)
  

## Table 2 Summary of posterior distributions
kable(r$summary.fixed[, c(1:5)], digits = 2)


## Fig. 3 Posterior distributions of fixed effects
posterior.plot(r$marginals.fixed[['mu']])
posterior.plot(r$marginals.fixed[['logDepth']])
posterior.plot(r$marginals.fixed[['Chl']])


## Fig. 4 Posterior mean and sd of the spatial effect
map.plot(inla.mesh.project(proj, r$summary.random$spatial[,"mean"]))
map.plot(inla.mesh.project(proj, r$summary.random$spatial[,"sd"]))


## Fig. 5 Posterior mean of the linear predictor (latent scale)
map.plot(inla.mesh.project(proj, 
                           r$summary.linear.predictor[pred.idx$data, "mean"]))


## Fig. 6 Posterior median predicted prob. of ocurrence (probability scale)
map.plot(inla.mesh.project(proj, 
                           r$summary.fitted.values[pred.idx$data, "0.5quant"]))


## Fig. 7 Posterior quartiles of prob. of ocurrence (probability scale)
map.plot(inla.mesh.project(proj, 
                           r$summary.fitted.values[pred.idx$data, "0.025quant"]))
map.plot(inla.mesh.project(proj, 
                           r$summary.fitted.values[pred.idx$data, "0.975quant"]))

