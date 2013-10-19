\name{normalize}
\alias{normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize a vector or matrix within a range.
}
\description{
Normalizes a vector or matrix values to spread between a minimum and maximum value.  By default, it normalizes between zero and one.
}
\usage{
normalize(x, low = 0, high = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector or matrix of numeric values.
}
  \item{low}{
The minimum value desired in the output data.
}
  \item{high}{
The maximum value desired in the output data.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector or matrix of normalized values.
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
k=matrix(runif(20),4,5)
normalize(k)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalize }
\keyword{ range }% __ONLY ONE__ keyword per line
