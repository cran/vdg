% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runif_cube.R, R/runif_sphere.R
\name{runif_cube}
\alias{runif_cube}
\alias{runif_sphere}
\title{Sampling for hyperspheres/hypercubes}
\usage{
runif_cube(n, m = 2, max.dist = 1, at = FALSE, nr.dist = 21)

runif_sphere(n, m = 2, max.radius = sqrt(m), at = FALSE, nr.rad = 21)
}
\arguments{
\item{n}{number of points to sample}

\item{m}{number of design factors}

\item{max.dist}{maximum distance from origin (L-infinity norm/supremum distance) for the hypercuboidal 
design region (enveloping hypercube)}

\item{at}{logical indicating whether to sample on concentric hyperspheres/hypercubes or not. With this
option \code{n} is distributed proportionally across radii / supremum distances so that the density
of samples on each concentric hypercube / hypersphere are uniform across the different hyperspheres / hypercubes..}

\item{nr.dist}{the number of concentric hypercubes to use in case at is \code{TRUE}}

\item{max.radius}{maximum radius of the hyperspherical design region (enveloping hypersphere)}

\item{nr.rad}{number of concentric hyperspheres to sample on in case of \code{at} being \code{TRUE}}
}
\description{
Sample uniformly in or on a hyperspheres or hypercubes.
}
\examples{

set.seed(1234)
runif_sphere(n = 10)

set.seed(1234)
samp <- runif_sphere(n = 500, at = TRUE)
plot(samp, asp = 1)
}
\references{
Pieter C. Schoonees, Niel J. le Roux, Roelof L.J. Coetzer (2016). Flexible Graphical Assessment of 
Experimental Designs in R: The vdg Package. \emph{Journal of Statistical Software}, 74(3), 1-22. 
\doi{10.18637/jss.v074.i03}.
}
\author{
Pieter C. Schoonees
}
