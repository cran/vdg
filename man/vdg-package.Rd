% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/_vdg.R
\name{vdg-package}
\alias{vdg-package}
\title{Variance Dispersion Graphs, Fraction-of-Design-Space Plots and Variants}
\description{
This package provides functionality for producing variance dispersion graphs (VDGs), 
fraction-of-design (FDS) plots and related graphics for assessing the prediction variance
properties of experimental designs. Random sampling is used to assess the distribution of the 
prediction variance throughout the design region. Multiple design and/or model formulae
can be assessed at the same time. Graphics are produced by the \pkg{ggplot2} package.
}
\details{
The workhorse function in the package is \code{\link{spv}}, which takes lists of 
experimental designs and / or model formulae and produces samples throughout the design region
at which the prediction variance is evaluated. Depending on the type of input for the 
\code{design} and \code{formula} arguments, \code{\link{spv}} creates output objects of S3 classes
\code{spv}, \code{spvlist}, \code{spvforlist} or \code{spvlistforlist}. The graphical output are
obtained with the \code{\link{plot}} methods of these classes, and the \code{which} argument
can be used to control the type of plots produced. 

The design regions allowed for are typically spherical or cuboidal in nature, but the 
\code{keepfun} argument to \code{\link{spv}} can be used for rejection sampling. In this way 
nonstandard design regions can be allowed for. See also the \code{type} argument of \code{\link{spv}}.
The output from the \code{\link{plot}} methods for objects created by \code{\link{spv}} are 
typically named lists of graphical objects created by \pkg{ggplot2}. These are best stored in an 
object and recreated by printing the required plot. Storing such graphical objects also enables 
post-hoc manipulation of the plots, such as changing the background colour by using 
\pkg{ggplot2}'s \code{\link{theme}} function.
}
\references{
Pieter C. Schoonees, Niel J. le Roux, Roelof L.J. Coetzer (2016). Flexible Graphical Assessment of 
Experimental Designs in R: The vdg Package. \emph{Journal of Statistical Software}, 74(3), 1-22. 
\doi{10.18637/jss.v074.i03}.
}
\seealso{
\code{\link{spv}}, \code{\link{plot.spv}}, and \code{vignette(topic = "vdg")}.
}
\keyword{package}
