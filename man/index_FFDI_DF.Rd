\name{index_FFDI_DF}
\alias{index_FFDI_DF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
McArthur Forest Fire Danger Index with precalculated drought factor.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_FFDI_DF(Temperature, DF, Humidity, Wind)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temperature}{
%%     ~~Describe \code{Temperature} here~~
}
  \item{DF}{
%%     ~~Describe \code{DF} here~~
}
  \item{Humidity}{
%%     ~~Describe \code{Humidity} here~~
}
  \item{Wind}{
%%     ~~Describe \code{Wind} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Temperature, DF, Humidity, Wind) 
{
    McA <- 0
    McA_list <- c()
    llength <- length(Temperature)
    for (day in 1:llength) {
        F <- 2 * exp(-0.45 + 0.987 * log(DF[day] + 0.001) - 0.0345 * 
            Humidity[day] + 0.0338 * Temperature[day] + 0.0234 * 
            Wind[day])
        McA_list <- append(McA_list, F)
    }
    McA_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
