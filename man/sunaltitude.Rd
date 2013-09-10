\name{sunaltitude}
\alias{sunaltitude}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Altitude of sun above horizon.
}
\description{
Calculates the altitude of the sun above the horizon for a given latitude, day of year and hour of day.
}
\usage{
sunaltitude(latitude, dayofyear, hour)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{latitude}{
Latitude of location, southern locations take negative latitude.
}
  \item{dayofyear}{
Day of year, julan day.
}
  \item{hour}{
Hour of day, from midnight.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Altitude (angle) of sun above the horizon.
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
\code{\link{dayintegral}},\code{\link{daylengths}},\code{\link{hourwatt}}
}
\examples{
latitude = -42
hour=11
doy=45
sunaltitude(latitude,doy,hour)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ altitude}
\keyword{ solar }% __ONLY ONE__ keyword per line
