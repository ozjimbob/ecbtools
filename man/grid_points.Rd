\name{grid_points}
\alias{grid_points}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a regular grid of points around a SpatialPointsDataFrame's points.
}
\description{
This function takes as input a set of geographic points, as a SpatialPointsDataFrame. For each point in the input object, it genereates a regular grid of points surrounding it, of chosen dimensions.  This is useful for setting up an array of repeat sampling locations, for example, for assessing vegetation condition or canopy cover, for a sequence of geographic locations.
}
\usage{
grid_points(points, spacing = 25, gridsize = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
A SpatialPointsDataFrame object, which can be loaded from a shapefile using the \code{\link{readShapePoints}} command of the package.
}
  \item{spacing}{
The spacing between lines of the grid. The default is 25.  This value is expressed in the same coordinate system as the input spatial data, so be aware that for data in geographic (latitude/longitude) projections, this value represents degrees.
}
  \item{gridsize}{
The number of rows/columns in the grid, for example, the default value of 3 will produce a 3x3 grid of nine points surrounding the centre point.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A SpatialPointsDataFrame object containing the grid of points.
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
data(SplashDams)
plot(SplashDams)
gridded = grid_points(SplashDams,0.03,3)
plot(gridded)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sampling }
\keyword{ gridded }% __ONLY ONE__ keyword per line
