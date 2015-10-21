# 
# # Load auxiliar functions from Finn
# source(file.path('code', 'leuk-demo', 'utils.R'))
# 
# # SPDE/GMRF model, (kappa^2-Delta)(tau x) = W:
# spde <- function(mesh, hyper.default=NULL){
#   s <- inla.spde2.matern(mesh, alpha = 2)
#   
#   # If a custom prior specification is given, overwrite the default
#   if(!is.null(hyper.default))
#     s$f$hyper.default <- hyper.default
#     
#   return(s)
# } 
# 
# ## Calculate mapping between triangulation vertices and grid points:
# meshproj <- function(mesh, xlm=xl, ylm=yl) UseMethod("meshproj")
# 
# meshproj.inla.mesh <- function(mesh, xlm=xl, ylm=yl) inla.mesh.projector(mesh, dims=c(200,200), xlim=xlm, ylim=ylm)
# 
# meshproj.grid <- function(mesh, xlm=xl, ylm=yl) return(mesh)
# 
# # Parameters of the SPDE model in terms of range and variance
# kappa <- function(range) sqrt(8)/range  # these relations are valid in dim=2.
# tau <- function(range, var=1.0) 1/sqrt(4*pi*kappa(range)^2*var)
# 
# # export the boundary of a mesh as a segment object
# # for use in the creation of another mesh
# get.bnd <- function(mesh){
#   return(inla.mesh.segment(mesh$loc[mesh$segm$bnd$idx[,1],]))
# }

# Returns a formula string with response variable Y, and intercept term
# and all the combinations of m elements in v
comb.terms <- function(m, v=sel.terms) {
  if(m==0) return('Y ~ 0 + mu')
  else {
    combis <- apply(combn(v, m), 2, paste, collapse=' + ')
    return(paste('Y ~ 0 + mu', combis, sep=' + '))
  }
}


runinla <- function(formula, stack, mesh, spde, cf, keep = FALSE)
{
  library(INLA)

  ## I need to change the environment of the formula to that of the present function,
  ## since INLA will look for spde in it.
  environment(formula) <- environment()

  ## Run INLA:
  r  <- inla(
    formula, 
    family="binomial",
    # Ntrials = 1, # default value - Bernoulli response
    data = inla.stack.data(stack),
    control.compute = list(
      return.marginals = TRUE, 
      dic = TRUE, 
      cpo = TRUE),
    control.predictor = list(
      compute = TRUE,
      A = inla.stack.A(stack),
      link = 1),
    ## Correction of the Laplace Approximation
    ## for Binomial data with many zeros
    ## http://arxiv.org/abs/1503.07307
    control.inla = list(
      correct = TRUE,
      correct.factor = 10,
      tolerance = 1e-5,  # We don't need to overoptimise:
      numint.maxfeval = 10e6),
    ## Verbose output:
    verbose = FALSE,
    keep = keep, 
    control.fixed = cf
  )
  return(r)
}



## Auxiliar functions for plotting

# Posterior distributions
posterior.plot <- function(mar) {
  par(mai=c(0.5, 0.1, 0.1, 0.1))
  #     curve(dnorm(x, mean=sum['mean'], sd=sum['sd']), xlim=sum[c(3,5)]+c(-1,1)*sum['sd'], xlab='', ylab='', yaxt='n')
  # I better use inla functions to plot marginals
  # but I don't like them to show so much unnecessary range
  crop <- .15    # percent to crop from each side
  shortmar <- lapply(inla.smarginal(mar), function(x) x[ceiling(length(x)*crop+1):floor(length(x)*(1-crop))])
  plot(inla.smarginal(shortmar), type='l', xlab='', ylab='', yaxt='n')
  segments(0, 0, y1=inla.dmarginal(0, mar), lwd=2, col='gray')
}


# Maps of projected results
map.plot <- function(plotdata, p=proj, palette=topo.colors, ...)
{
  bbb = (levelplot(row.values=p$x, column.values=p$y, x=plotdata,
                   mm=coord_coast, pp=obs, panel=levelplotmap,
                   col.regions=palette, 
                   xlim=range(p$x), ylim=range(p$y), aspect="iso",
                   contour=TRUE, cuts=11, labels=FALSE, pretty=TRUE,
                   scales=list(draw=FALSE),
                   xlab='',ylab='', ...))
  print(bbb)
}

# Panel function
levelplotmap = function(..., mm, pp, colland='white') {
  panel.levelplot(...)
  panel.polygon(mm, lwd=3, col=colland, add=T)
  #    panel.points(x=coordinates(pp), pch=20, col=pp$Presencia+1)
  #    legend("topleft", c("Absense", "Presence"), col=1:2, pch=20)
}
