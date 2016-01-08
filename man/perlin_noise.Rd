\name{perlin_noise}
\alias{perlin_noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a matrix of Perlin noise
}
\description{
Generates a matrix of Perlin noise of defined dimensions
}
\usage{
perlin_noise(vx = 7, vy = 8, ix = 100, iy = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vx}{
X dimension of field size for vector noise
}
  \item{vy}{
Y dimension of field size for vector noise
}
  \item{ix}{
X dimension of output matrix
}
  \item{iy}{
Y dimension of output matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

A matrix of values of dimension ix,iy.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Grant Williamson 
}
\note{
Based on code from \url{http://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
i = perlin_noise()
image(i)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~noise }
\keyword{ ~perlin }% __ONLY ONE__ keyword per line
