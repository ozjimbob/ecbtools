\name{index_CFWI}
\alias{index_CFWI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Canadian Fire Weather Index
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_CFWI(Month, Days, Temp, Dew, WS, Rain, daylist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Month}{
Vector of decimal months of observations.
}
  \item{Days}{
Vector of Julian days of observations.
}
  \item{Temp}{
Vector of temperature (C)
}
  \item{Dew}{
Vector of dew Point (C)
}
  \item{WS}{
Wind speed (km/h)
}
  \item{Rain}{
Vector of rain (mm)
}
  \item{daylist}{
Day length vector for each day of year.
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
Van Wagner, C. E., & Pickett, T. L. (1985). Equations and FORTRAN Program for the Canadian Forest Fire Weather Index System (Vol. 33).
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Month, Days, Temp, Dew, WS, Rain, daylist) 
{
    Humid <- calc_humid(Temp, Dew)
    DLA <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    CFWI <- 0
    F <- 85
    P <- 6
    D <- 15
    CFWI_list <- c()
    llength <- length(Temp)
    for (day in 1:llength) {
        T = Temp[day]
        H = Humid[day]
        if (H > 100) {
            H = 100
        }
        W = WS[day]
        if (W < 0) {
            W = 0
        }
        ro = Rain[day]
        if (ro < 0) {
            ro = 0
        }
        Le = daylist[Days[day]]
        Fo = F
        mo = 147.2 * (101 - Fo)/(59.5 + Fo)
        if (ro > 0.5) {
            rf = ro - 0.5
            if (mo <= 150) {
                mr = mo + 42.5 * rf * (exp(-100/(251 - mo))) * 
                  (1 - exp(-6.93/rf))
            }
            else {
                mr = mo + 42.5 * rf * (exp(-100/(251 - mo))) * 
                  (1 - exp(-6.93/rf)) + 0.0015 * ((mo - 150)^2) * 
                  (rf^0.5)
            }
            if (mr > 250) {
                mr = 250
            }
            mo = mr
        }
        Ed = 0.942 * H^0.679 + 11 * exp((H - 100)/10) + 0.18 * 
            (21.1 - T) * (1 - exp(-0.115 * H))
        if (mo > Ed) {
            ko = 0.424 * (1 - (H/100)^1.7) + 0.0694 * W^0.5 * 
                (1 - (H/100)^8)
            kd = ko * 0.581 * exp(0.0365 * T)
            m = Ed + (mo - Ed) * 10^(-kd)
        }
        Ew = 0.618 * H^0.753 + 10 * exp((H - 100)/10) + 0.18 * 
            (21.1 - T) * (1 - exp(-0.115 * H))
        if (mo < Ew & mo < Ed) {
            kl = 0.424 * (1 - ((100 - H)/100)^1.7) + 0.0694 * 
                W^0.5 * (1 - ((100/H)/100)^8)
            kw = kl * 0.581 * exp(0.0365 * T)
            m = Ew - (Ew - mo) * 10^(-kw)
        }
        if (Ed >= mo & mo >= Ew) {
            m = mo
        }
        F = 59.5 * (250 - m)/(147.2 + m)
        Po = P
        if (ro > 1.5) {
            re = 0.92 * ro - 1.27
            Mo = 20 + exp(5.6348 - Po/43.43)
            if (Po <= 33) {
                b = 100/(0.5 + 0.3 * Po)
            }
            if (Po > 33 & Po <= 65) {
                b = 14 - 1.3 * log(Po)
            }
            if (Po > 65) {
                b = 6.2 * log(Po) - 17.2
            }
            Mr = Mo + 1000 * re/(48.77 + b * re)
            Pr = 244.72 - 43.43 * log(Mr - 20)
            if (Pr < 0) {
                Pr = 0
            }
            Po = Pr
        }
        if (T < -1.1) {
            TT = -1.1
        }
        else {
            TT = T
        }
        K = 1.894 * (TT + 1.1) * (100 - H) * Le * 1e-05
        P = Po + 100 * K
        Do = D
        if (ro > 2.8) {
            rd = 0.83 * ro - 1.27
            Qo = 800 * exp(-Do/400)
            Qr = Qo + 3.937 * rd
            Dr = 400 * log(800/Qr)
            if (Dr < 0) {
                Dr = 0
            }
            Do = Dr
        }
        if (T < 2.8) {
            TT = 2.8
        }
        else {
            TT = T
        }
        V = 0.36 * (TT + 2.8)
        D = Do + 0.5 * V
        fW = exp(0.05039 * W)
        fF = 91.9 * exp(-0.1386 * m) * (1 + m^5.31/(4.93 * 1e+08))
        R = 0.208 * fW * fF
        if (P <= 0.4 * D) {
            U = 0.8 * P * D/(P + 0.4 * D)
        }
        if (P > 0.4 * D) {
            U = P - (1 - 0.8 * D/(P + 0.4 * D)) * (0.92 + (0.0114 * 
                P)^1.7)
        }
        if (U <= 80) {
            fD = 0.626 * U^0.809 + 2
        }
        if (U > 80) {
            fD = 1000/(25 + 108.64 * exp(-0.023 * U))
        }
        B = 0.1 * R * fD
        if (is.na(B)) {
            S = 0
        }
        else {
            if (B > 1) {
                S = exp(2.72 * (0.434 * log(B))^0.647)
            }
            if (B <= 1) {
                S = B
            }
        }
        CFWI_list = append(CFWI_list, S)
    }
    CFWI_list
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~fire }
\keyword{ ~cfwi }% __ONLY ONE__ keyword per line
