\name{index_FFDI_G4}
\alias{index_FFDI_G4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
McArthur Grassland Fire Index Mk4.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_FFDI_G4(Temperature, DewPoint, Wind)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temperature}{
%%     ~~Describe \code{Temperature} here~~
}
  \item{DewPoint}{
%%     ~~Describe \code{DewPoint} here~~
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
function (Temperature, DewPoint, Wind) 
{
    McA <- 0
    McA_list <- c()
    llength <- length(Temperature)
    for (day in 1:llength) {
        satvp <- 6.11 * 10^((7.5 * Temperature[day]/(237.7 + 
            Temperature[day])))
        actvp = 6.11 * 10^((7.5 * DewPoint[day]/(237.7 + DewPoint[day])))
        Humid <- (actvp/satvp) * 100
        H <- Humid
        T <- Temperature[day]
        U <- Wind[day]
        C = 100
        Q = 4.5
        F = exp(-1.523 + 1.027 * log(Q) - 0.009432 * (100 - C)^1.535 + 
            0.02764 * T - 0.2205 * sqrt(H) + 0.6422 * sqrt(U))
        McA_list <- append(McA_list, F)
    }
    McA_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
