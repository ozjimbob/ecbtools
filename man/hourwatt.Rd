\name{hourwatt}
\alias{hourwatt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the energy (Watts) of sunlight reaching the earth for a given day of year, latitude and hour of day.
}
\description{
Calculates the energy (Watts) of sunlight reaching the earth for a given day of year, latitude and hour of day.
}
\usage{
hourwatt(dayofyear, latitude, hour)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dayofyear}{
Day of year, julian day.
}
  \item{latitude}{
Latitude of location, southern locations take negative latitudes.
}
  \item{hour}{
Hour of day, from midnight.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Watts per square meter of radiation received.
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
\code{\link{dayintegral}},\code{\link{daylengths}},\code{\link{sunaltitude}}
}
\examples{
latitude = -42
hour=11
doy=20
hourwatt(doy,latitude,hour)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ energy }
\keyword{ solar }% __ONLY ONE__ keyword per line
\keyword{ radiation }% __ONLY ONE__ keyword per line\
\keyword{ watts }% __ONLY ONE__ keyword per line