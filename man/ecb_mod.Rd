\name{ecb_mod}
\alias{ecb_mod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ECB Linear Modelling Wizard
}
\description{
Simple wrapper for running lm, glm, lmer, glmer and nb.glmer models on a full suite of variable combinations, outputing an AIC table and w+ variable importance table.
}
\usage{
ecb_mod(data, term.names, response, random.def, model.family = gaussian(), model.type = "lm")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame containing the model input.
}
  \item{term.names}{
Vector containing a list of the variables to be included in the model.
}
  \item{response}{
Character string containing the response variable of the model.
}
  \item{random.def}{
If running a mixed model, a character string containing the random effect part of the model formula, eg. "(1|Site)".
}
  \item{model.family}{
Model family function, see \link{family}.
}
  \item{model.type}{
String defining the model type, can be lm, glm, lmer, glmer or nb.glmer.
}
}
\details{
Details
}
\value{
Returns a list containing the AIC table and w+ table.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Grant Williamson (grant.williamson@utas.edu.au)
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

data=read.csv("myfile.csv")

ecb_mod(data, term.names=c("MAP","Slope","Soil"), response="Height", random.def="(1|Site)", model.family = poisson(), model.type = "lmer")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
