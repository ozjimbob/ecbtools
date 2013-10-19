\name{raster.kmeans}
\alias{raster.kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
K-means Classification on a RasterStack
}
\description{
This function performs a k-means unsupervised classification of a stack of rasters.  Number of clusters can be specified, as well as number of iterations and starting sets.  An optional geographic weighting system can be turned on that constrains clusters to a geographic area, by including coordinates in the clustering.  All variables are normalized before clustering is performed.
}
\usage{
raster.kmeans(x, k = 12, iter.max = 100, nstart = 10, geo = T, geo.weight = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A RasterStack object, or a string pointing to the directory where all raster layers are stored.  All rasters should have the same extent, resolution and coordinate system.
}
  \item{k}{
Number of clusters to classify to.
}
  \item{iter.max}{
Maximum number of iterations allowed
}
  \item{nstart}{
Number of random sets to be chosen.
}
  \item{geo}{
True/False - should geographic weighting be used?
}
  \item{geo.weight}{
A weighting multiplier indicating the strength of geographic weighting relative to the other variables.  A value of 1 gives equal weight.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A raster object, of the same extent and coordinates as the input data, containing the cluster classification.
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
\keyword{ kmeans }
\keyword{ classify }
\keyword{ cluster }
\keyword{ raster }% __ONLY ONE__ keyword per line

