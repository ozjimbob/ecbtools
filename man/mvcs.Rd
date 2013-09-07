\name{mvcs}
\alias{mvcs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Continuous Stratification
}
\description{
Selects a subsample of a data frame where the distribution of values of a set of variables in the subset match those of the input data frame.
}
\usage{
mvcs(data, number, variables, iter = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data frame to be subsampled.
}
  \item{number}{
The number of rows of data to be returned, ie. the sample size.
}
  \item{variables}{
A vector containing the names of variables to be used for the stratification.
}
  \item{iter}{
The number of iterations to be performed to find the optimal subset.
}
}
\details{
This function uses a Cramer test to select rows of an input data frame, where the distribution of values across multiple variables closely matches the distribution of input data.  For example, you may have a large input data set of geographic points, where elevation is skewed towards low values, and rainfall is skewed towards high values.  This function returns a subset of the data, of chosen size, where the rainfall and elevation have distributions matching those of the input data.

This function works only on continuous variables, for which a statistical distribution can be calculated.  For stratification within factors, use \code{\link{mvs}}.
}
\value{
A data frame.
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
\code{\link{mvs}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stratification }
\keyword{ stratify }% __ONLY ONE__ keyword per line
\keyword{ subsample }% __ONLY ONE__ keyword per line
\keyword{ cramer }% __ONLY ONE__ keyword per line