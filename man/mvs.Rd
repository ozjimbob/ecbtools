\name{mvs}
\alias{mvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Stratification
}
\description{
Selects a subsample of a data frame where the proportion of points in each class of a set of factor variables in the subset match those of the input data frame.
}
\usage{
mvs(data, number, variables, iter = 20)
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
A vector containing the names of factor variables to be used in the stratification
}
  \item{iter}{
The number of iterations to try to find the optimal solution.
}
}
\details{
This function uses a chi square test to select rows of an input data frame, where the proportion of values falling in classes of multiple factors closely matches the propotion in the input data.  For example, you may have a large input data set of geographic points with variables describing soil type or land use.  This function returns a subset of the data, of chosen size, the proportion of soil types and land uses matches those of the input data.

This function works only on factors for which a proportion of counts in each class can be calculated.  For stratification with continuous variables, use \code{\link{mvcs}}.
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
\code{\link{mvcs}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stratification }
\keyword{ stratify }% __ONLY ONE__ keyword per line
\keyword{ factor }% __ONLY ONE__ keyword per line
