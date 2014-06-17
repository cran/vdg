#' Sampler Function
#' 
#' This is a wrapper for the sampling funcions of the \pkg{vdg} package. It extracts design properties from the 
#' design passed to it to take appropriate samples.
#' 
#' @param n number of points to sample
#' @param design  design for which the sample is required (either a matrix or data frame)
#' @param type type of design region/sampling method. One of "spherical", "cuboidal", "lhs", "mlhs", "slhs", or "rslhs"
#' @param at logical; should sampling be done on the surface of hyperspheres or hypercubes? Not used for LHS methods.
#' @param \dots other arguments passed to the underlying sampling functions.
#' @return Matrix with samples as rows, with S3 class \code{smpl}
#' @seealso \code{\link{runif_sphere}}, \code{\link{runif_cube}}, \code{\link{LHS}}, 
#' \code{\link{MLHS}}, \code{\link{SLHS}}, \code{\link{RSLHS}}
#' @author Pieter C. Schoonees
#' @examples
#' 
#' set.seed(1896)
#' sampler(n = 10, design = expand.grid(x = -1:1, y = -1:1))
sampler <-
function(n, design, type = "spherical", at = FALSE, ...){
  m <- ncol(design) 
  samp <- switch(type, spherical = runif_sphere(n = n, m = m, at = at, ...), 
                 cuboidal = runif_cube(n = n, m = m, at = at, ...), 
                 lhs = LHS(n = n, m = m, ...), 
                 mlhs = MLHS(n = n, m = m, ...),
                 slhs = SLHS(n = n, m = m, ...), 
                 rslhs = RSLHS(n = n, m = m, ...))
  colnames(samp) <- colnames(design)
  class(samp) <- c("smpl", "matrix")
  samp
}
