library(INLA)
inla.setOption(pardiso.license=NULL)
library(knitr)

source('R/spde-book-functions.R')

# knitr options
opts_chunk$set( 
  tidy = FALSE,
  collapse = TRUE, # colapse chunks 
  fig.width = 7,
  fig.height = 7,
  out.width = "97%",
  fig.align = "center",
  message = FALSE, 
  warning = FALSE,
#  background = '#99E6FF', ## from rgb(.6,9,1) 
  tidy.opts = list(blank = FALSE, width.cutoff = 65)
)

knit_hooks$set(par.args = function(before, options, envir) {
  if (before) graphics::par(options$par.args)
})

options(width = 63, prompt = " ", continue = "   ", digits = 4)

