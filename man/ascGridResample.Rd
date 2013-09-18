\name{ascGridResample}
\alias{ascGridResample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crops, resamples and registers an input ASCII Grid raster to a template.
}
\description{
This function is useful for generating a set of precisely matching ASCII Grid raster files, for input to software such as MaxEnt. The function takes an input raster to adjust, a template raster (which must be an ASCII Grid .asc file), and a filename for output, and crops, resamples and aligns the input raster to exactly match the template. 
}
\usage{
ascGridResample(x, template, filename = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The input raster to be cropped.  Either a filename, or a RasterLayer object.
}
  \item{template}{
The template raster.  Either a filename or a RasterLayer object.  Must be an ASCII Grid (.asc) file.
}
  \item{filename}{
The filename the output raster will be written to.
}
}
\details{
 The standard crop/extent functions in the raster package seem to have issues with coordinate precision, so software such as MaxEnt often complains of projection errors even if the number of cells in the raster input files are the same.  In order to overcome this, the function manually copies the header data from the template file to the output file.  All rasters should be in WGS84 projection.  
}
\value{
The function writes the output ASCII Grid raster to disk, and returns the filename.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Grant Williamson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ raster }
\keyword{ maxent }% __ONLY ONE__ keyword per line
\keyword{ crop }% __ONLY ONE__ keyword per line
\keyword{ resample }% __ONLY ONE__ keyword per line