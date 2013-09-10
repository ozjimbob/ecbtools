\name{index_Sharples}
\alias{index_Sharples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sharples Fire Danger Index.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_Sharples(Temp, Wind, Dew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Temp}{
%%     ~~Describe \code{Temp} here~~
}
  \item{Wind}{
%%     ~~Describe \code{Wind} here~~
}
  \item{Dew}{
%%     ~~Describe \code{Dew} here~~
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
function (Temp, Wind, Dew) 
{
    humid <- calc_humid(Temp, Dew)
    FFI_list <- c()
    llength <- length(Temp)
    for (day in 1:llength) {
        FMI = 10 - 0.25 * (Temp[day] - humid[day])
        f = max(Wind[day], 1)/FMI
        FFI_list <- append(FFI_list, f)
    }
    FFI_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
