\name{index_Fosberg}
\alias{index_Fosberg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fosberg Fire Weather Index.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_Fosberg(Temp, Wind, Dew)
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
    Wind <- Wind/1.6
    Temp = Temp * (9/5) + 32
    FFI_list <- c()
    llength <- length(Temp)
    for (day in 1:llength) {
        if (is.na(humid[day])) {
            FFI_list <- append(FFI_list, NA)
        }
        else {
            if (humid[day] < 10) {
                m = 0.03229 + 0.281073 * humid[day] - 0.000578 * 
                  humid[day] * Temp[day]
            }
            if (humid[day] >= 10 & humid[day] <= 50) {
                m = 2.22749 + 0.160107 * humid[day] - 0.01478 * 
                  Temp[day]
            }
            if (humid[day] > 50) {
                m = 21.0606 + 0.005565 * humid[day]^2 - 0.00035 * 
                  humid[day] * Temp[day] - 0.483199 * humid[day]
            }
            n = 1 - 2 * (m/30) + 1.5 * (m/30)^2 - 0.5 * (m/30)^3
            FFWI = n * ((1 + Wind[day]^2)^0.5)/0.3002
            FFI_list <- append(FFI_list, FFWI)
        }
    }
    FFI_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
