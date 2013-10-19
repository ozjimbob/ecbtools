\name{index_NI}
\alias{index_NI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nesterov Fire Danger Index.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_NI(DewPoint, Temperature, Rain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DewPoint}{
%%     ~~Describe \code{DewPoint} here~~
}
  \item{Temperature}{
%%     ~~Describe \code{Temperature} here~~
}
  \item{Rain}{
%%     ~~Describe \code{Rain} here~~
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
function (DewPoint, Temperature, Rain) 
{
    Deficit <- Temperature - DewPoint
    G <- 0
    G_list <- c()
    llength <- length(DewPoint)
    for (day in 1:llength) {
        if (Rain[day] > 3) {
            G <- 0
            G_list <- append(G_list, G)
            next
        }
        if (Temperature[day] > 0) {
            G <- G + (Temperature[day] * Deficit[day])
        }
        G_list <- append(G_list, G)
    }
    G_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
