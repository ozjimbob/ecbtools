\name{calc_dew}
\alias{calc_dew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the dew point, based on temperature and relative humidity.
}
\description{
A function to convert temperature and relative humidity data into a dew point value.
}
\usage{
calc_dew(T, RH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
Temperature (°C)
}
  \item{RH}{
Relative humidity (\%)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A variable or vector containing dew point values in °C.
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
\code{\link{calc_humid}}
}
\examples{

Temp=34
RH=56
dp = calc_dew(Temp,RH)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Temperature }
\keyword{ Humidity }% __ONLY ONE__ keyword per line
\keyword{ Dew }