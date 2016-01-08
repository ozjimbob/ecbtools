\name{calc_humid}
\alias{calc_humid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate relative humidity from temperature and dew point.
}
\description{
Calculates a relative humidity value, based on input temperature and dew point values.
}
\usage{
calc_humid(T, Dp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
Temperature (C)
}
  \item{Dp}{
Dew point (C)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A variable or vector containing relative humidity (\%) values.
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
\code{\link{calc_dew}}
}
\examples{
Temp=23.5
DP=14.3
RH=calc_humid(Temp,DP)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Temperature }
\keyword{ Humidity }% __ONLY ONE__ keyword per line
\keyword{ Dew }% __ONLY ONE__ keyword per line