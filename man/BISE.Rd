\name{BISE}
\alias{BISE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Best Index Slope Extraction
}
\description{
Best Index Slope Extraction (BISE) function for removing noise from satellite NDVI data.
}
\usage{
BISE(x, slide_period = 20, slope_threshold=0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A sequence of equally-spaced values, for example, from daily NDVI or meteorological data.
}
  \item{slide_period}{
Days to look ahead for decrease.
}
  \item{slope_threshold}{
Percentage decrease threshold.
}
}
\details{
This is a simple implementation of the BISE algorithm for reducing noise in time-series data where sudden drops in value (for example, due to cloud in satellite imagery) should be ignored. It is based on Viovy et al. 1992.
}
\value{
A vector of corrected data.
}
\references{
Viovy, N., Arino, O., & Belwards, A.S. (1992) The Best Index Slope Extraction (BISE): A method for reducing noise in NDVI time-series. International Journal of Remote Sensing. 13:8, 1585-1590.
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
h=c(seq(10,90,2),seq(90,10,-2))
h=h+rnorm(82,0,10)
o=BISE(h)
plot(h)
lines(o)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~BISE}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
