\name{pointDistSample}
\alias{pointDistSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subsamples spatial points with bias against dense points.
}
\description{
This function takes a subsample of points in a SpatialPointsDataFrame object.  The subsample is biased towards removing points that are close to other points, therefore reducing the density of highly dense areas.  This may be useful in bias adjusting species distribution patterns away from densely sampled areas.  
}
\usage{
pointDistSample(points, sample.size,id=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{
A SpatialPointsDataFrame object.
}
  \item{sample.size}{
The number of points to retain in the sample.
}
 \item{id}{
Square-root weighting.
}
}
\details{
Distance from each point to each other point is calculated, and a mean neighbour distance for each point is determined.  Distances are normalized, and points with a greater mean distance to other points have a greater probability of being sampled.
}
\value{
A SpatialPointsData Frame.
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
\keyword{ sample }
\keyword{ distribution }% __ONLY ONE__ keyword per line
