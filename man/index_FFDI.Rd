\name{index_FFDI}
\alias{index_FFDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
McArthur Forest Fire Danger Index Mk5
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_FFDI(Temperature, Rain, DewPoint, Wind, KBDI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temperature}{
%%     ~~Describe \code{Temperature} here~~
}
  \item{Rain}{
%%     ~~Describe \code{Rain} here~~
}
  \item{DewPoint}{
%%     ~~Describe \code{DewPoint} here~~
}
  \item{Wind}{
%%     ~~Describe \code{Wind} here~~
}
  \item{KBDI}{
%%     ~~Describe \code{KBDI} here~~
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
function (Temperature, Rain, DewPoint, Wind, KBDI) 
{
    McA <- 0
    McA_list <- c()
    llength <- length(Temperature)
    DSinceRain <- 0
    LastRain <- 0
    for (day in 1:llength) {
        if (Rain[day] > 0) {
            LastRain <- Rain[day]
            DSinceRain <- 0
        }
        else {
            DSinceRain <- DSinceRain + 1
        }
        satvp <- 6.11 * 10^((7.5 * Temperature[day]/(237.7 + 
            Temperature[day])))
        actvp = 6.11 * 10^((7.5 * DewPoint[day]/(237.7 + DewPoint[day])))
        Humid <- (actvp/satvp) * 100
        DF <- (0.191 * (KBDI[day] + 104) * (DSinceRain + 1)^1.5)/((3.52 * 
            (DSinceRain + 1)^1.5) + LastRain - 1)
        if (DF > 10) {
            DF <- 10
        }
        F <- 2 * exp(-0.45 + 0.987 * log(DF + 0.001) - 0.0345 * 
            Humid + 0.0338 * Temperature[day] + 0.0234 * Wind[day])
        McA_list <- append(McA_list, F)
    }
    McA_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
