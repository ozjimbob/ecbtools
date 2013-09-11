\name{timespline}
\alias{timespline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
2D smoothed plots of data by year and day of year.
}
\description{
This function produces a 2D field plot of a value, with the x-axis representing year, and the y-axis representing day of year (julian day).  A raw plot, and a plot with a 2D spline smooth are produced.  These plots are useful for visually examining seasonal trends and cycles in meteorological variables over long time periods.
}
\usage{
timespline(data, datevar, plotvar, k = 40, hinge = F, col = div_pal(256))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the input data.
}
  \item{datevar}{
The name of the column in the data frame containing the date of the observation.  This should contain a Date-type value, as formatted by the \code{\link{as.Date}} function.
}
  \item{plotvar}{
The name of column in the data frame containing the numeric variable to plot.
}
  \item{k}{
The degree of spline smoothing to be performed.  Higher values give more detailed plots, but may overfit.
}
  \item{hinge}{
True/False - should the data be hinged over the middle of the year?  If false, the y-axis of the plot represents day of year from 1 to 365.  If true, values range from 0 to 1 representing "winterness" and "summerness".
}
  \item{col}{
A palette specification.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs a plot with two panels, one with raw data, and one with 2D spline-smoothed values.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spline }
\keyword{ smooth }% __ONLY ONE__ keyword per line
