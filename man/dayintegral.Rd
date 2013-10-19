\name{dayintegral}
\alias{dayintegral}

\title{
Time sun is above horizon for a day and latitude.
}
\description{
A coarse but quick function to calculate the number of hours (down to approximately 2-minute resolution) that the sun is above the horizon for a given day of the year, and latitude.
}
\usage{
dayintegral(dayofyear, latitude)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dayofyear}{
Day of year (Julian day).
}
  \item{latitude}{
Latitude, southern latitudes take negative values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The number of hours, as a floating point value, the sun was above the horizon.
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
\code{\link{daylengths}},\code{\link{hourwatt}},\code{\link{sunaltitude}}
}
\examples{
latitude = -42
doy = 23
dayintegral(doy,latitude)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ solar }
\keyword{ sunlight }% __ONLY ONE__ keyword per line
\keyword{ hours }