\name{get_et}
\alias{get_et}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SDI Evapotranspiration Lookup Function
}
\description{
SDI Evapotranspiration Lookup Function
}
\usage{
get_et(Temp, SDI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temp}{
Daily maximum temperature in Celcius.
}
  \item{SDI}{
Soil Dryness Index in mm.
}
}
\details{
This is a lookup table for Hobart, Tasmania, based on Mount 1972.
}
\value{
A daily evapotranspiration amount (mm).
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
