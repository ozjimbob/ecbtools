\name{ecbtools-package}
\alias{ecbtools-package}
\alias{ecbtools}
\docType{package}
\title{
Environmental Change Biology tools.
}
\description{
This package provides functions developed by the Environmental Change Biology lab at the University of Tasmania, including functions to determine canopy cover from hemispherical/fish eye photos, and stratified sample plot placement functions.
}
\details{
\tabular{ll}{
Package: \tab ecbtools\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-09-05\cr
License: \tab BSD 2-clause\cr
}
This package provides two functions useful for processing hemispherical canopy cover photos. \code{\link{darkpix}}, which processes a single photographic image file, extracts the central circle representing the hemispherical horizon, and calculates the percentage of dark (cover) pixels, and \code{\link{darkpix.dir}}, which applies the \code{\link{darkpix}} function to an entire directory of images, and returns a data frame.
}
\author{
Grant Williamson
Maintainer: Grant Williamson <grant.williamson@utas.edu.au>
}

\keyword{ package }
\keyword{ cover }
\keyword{ hemispherical }
\keyword{ photo }
\keyword{ canopy }


