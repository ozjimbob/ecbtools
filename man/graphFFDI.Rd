\name{graphFFDI}
\alias{graphFFDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graph of McArthur FFDI for selected Australian Met Station
}
\description{
This function retrieves current and recent observations for a given Australian meteorological station, plots a graph of Forest Fire Danger Index against time, and returns a data frame containing the calculated values.
}
\usage{
graphFFDI(state = "N", station = "94757", DF = 9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{state}{
A single character indicating the state the station is in. T=Tasmania, V=Victoria, N=New South Wales, Q=Queensland, S=South Australia, W=Western Australia, D=Northern Territory, A=Australian Capital Territory
}
  \item{station}{
The WMO station number for the station.
}
  \item{DF}{
The drought factor value, from 0-10.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots a graph with fire danger ratings indicated by colour, and returns a data frame of observations.
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
FFDI=graphFFDI("T",94975,7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FFDI }
\keyword{ fire }% __ONLY ONE__ keyword per line
