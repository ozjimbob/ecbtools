\name{index_SDI}
\alias{index_SDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate SDI Soil Dryness Index
}
\description{
Calculate Mount's Soil Dryness Index
}
\usage{
index_SDI(Rain, Temperature)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rain}{
Vector of daily precipitation in mm.
}
  \item{Temperature}{
Vector of daily maximum temperature in Celcius.
}
}
\details{
This function calculates the daily SDI (in mm) from rainfall and temperature.  Evapotranspiration calculation is currently based on  a lookup table for Hobart, Tasmania, based on Mount's original report. 
}
\value{
A vector of Soil Dryness Index values.
}
\references{
Mount, A. B. 1972. The derivation and testing of a soil dryness index using run-off data. Bulletin No. 4, Tasmanian Forestry Commission.
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
