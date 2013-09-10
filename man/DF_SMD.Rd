\name{DF_SMD}
\alias{DF_SMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a drought factor vector from rain and soil moisture deficit.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DF_SMD(rain, SMD_vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rain}{
%%     ~~Describe \code{rain} here~~
}
  \item{SMD_vec}{
%%     ~~Describe \code{SMD_vec} here~~
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
function (rain, SMD_vec) 
{
    conseq_days_rain = rep(0, length(rain))
    conseq_sum_rain = rep(0, length(rain))
    size_of_largest = rep(0, length(rain))
    days_since_largest = rep(0, length(rain))
    rain[is.na(rain)] = 0
    conseq = 0
    conseq_sum = 0
    rfvec <- rain
    for (d in 2:length(rfvec)) {
        prev_rain = rfvec[d - 1]
        cur_rain = rfvec[d]
        if (cur_rain > 2) {
            conseq = conseq + 1
            conseq_sum = conseq_sum + cur_rain
        }
        if (cur_rain <= 2 & conseq != 0) {
            conseq = 0
            conseq_sum = 0
        }
        biggest = max(rfvec[(d - conseq):d])
        which_day = which.max(rev(rfvec[(d - conseq):d])) - 1
        conseq_days_rain[d] = conseq
        size_of_largest[d] = biggest
        if (conseq > 0) {
            days_since_largest[d] = which_day
        }
        if (conseq == 0) {
            days_since_largest[d] = days_since_largest[d - 1] + 
                1
        }
        if (conseq > 0) {
            conseq_sum_rain[d] = conseq_sum
        }
        if (conseq == 0) {
            conseq_sum_rain[d] = conseq_sum_rain[d - 1]
        }
    }
    look_frame = data.frame(rfvec, conseq_days_rain, conseq_sum_rain, 
        size_of_largest, days_since_largest)
    look_frame$conseq_sum_rain[look_frame$days_since_larges > 
        20] = 0
    N_vec = look_frame$days_since_largest
    P_vec = look_frame$conseq_sum_rain
    x_vec = rep(0, length(P_vec))
    for (a in 1:length(N_vec)) {
        N = N_vec[a]
        P = P_vec[a]
        if (N >= 1 & P >= 2) {
            x = N^1.3/(N^1.3 + P - 2)
        }
        if (N == 0 & P >= 2) {
            x = 0.7481988/(0.7481988 + P - 2)
        }
        if (P < 2) {
            x = 1
        }
        x_vec[a] = x
    }
    xlim_vec = rep(0, length(SMD_vec))
    for (a in 1:length(SMD_vec)) {
        SMD = SMD_vec[a]
        if (is.na(SMD)) {
            next
        }
        if (SMD < 20) {
            xlim = 1/(1 + 0.1135 * SMD)
        }
        if (SMD >= 20) {
            xlim = 75/(270.525 - 1.267 * SMD)
        }
        if (x_vec[a] > xlim) {
            x_vec[a] = xlim
        }
    }
    DF_vec = rep(0, length(SMD_vec))
    for (a in 1:length(DF_vec)) {
        SMD = SMD_vec[a]
        x = x_vec[a]
        if (is.na(SMD)) {
            next
        }
        DF = 10.5 * (1 - exp(-(SMD + 30)/40)) * (41 * x^2 + x)/(40 * 
            x^2 + x + 1)
        DF_vec[a] = DF
    }
    DF_vec
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
