\name{daylengths}
\alias{daylengths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate daily lengths of day for a given latitude.
}
\description{
Calculates a vector of length 366 containing the number of hours the sun was above the horizon for each day of the year.
}
\usage{
daylengths(latitude)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{latitude}{
Latitude of site, southern latitudes take negative values.
}
}
\details{
This is an approximate function, and day lengths are calculated to a resolution of approximately 2 minutes.
}
\value{
A vector of day lengths in hours.
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
\code{\link{dayintegral}},\code{\link{hourwatt}},\code{\link{sunaltitude}}
}
\examples{
latitude = -42
l=daylengths(latitude)
plot(l,type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ solar }
\keyword{ hours }% __ONLY ONE__ keyword per line
\keyword{ day }% __ONLY ONE__ keyword per line