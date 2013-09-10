\name{index_KBDI}
\alias{index_KBDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Keech-Byram Drought Index.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_KBDI(Temperature, Rain, MAP = sum(Rain)/length(Rain/365.25))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temperature}{
%%     ~~Describe \code{Temperature} here~~
}
  \item{Rain}{
%%     ~~Describe \code{Rain} here~~
}
  \item{MAP}{
%%     ~~Describe \code{MAP} here~~
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
function (Temperature, Rain, MAP = sum(Rain)/length(Rain/365.25)) 
{
    KBDI <- 200
    KBDI_list <- c()
    llength <- length(Temperature)
    for (day in 1:llength) {
        pKBDI <- KBDI
        KBDI <- KBDI - Rain[day]
        if (KBDI < 0) {
            KBDI <- 0
        }
        todayET <- lookupET(pKBDI, Temperature[day], MAP)
        KBDI <- KBDI + todayET
        KBDI_list <- append(KBDI_list, KBDI)
    }
    KBDI_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
