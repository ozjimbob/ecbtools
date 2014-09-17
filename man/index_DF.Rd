\name{index_DF}
\alias{index_DF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the FFDI Drought Factor from input meteorological data.
}
\description{
Calculates the FFDI Drought Factor from input meteorological data.
}
\usage{
index_DF(Temperature, Rain, KBDI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temperature}{
A vector of daily maximum temperatures in degrees C.
}
  \item{Rain}{
A vector of daily rainfall amounts in mm.
}
  \item{KBDI}{
A daily record of Keech-Bryam Drought Index values.
}
}
\details{
Calculates drought factor.
}
\value{
Returns a vector of Drought Factor values, from 0 to 10.
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
index_DF(temp,rain,k)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ drought }

