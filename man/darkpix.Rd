\name{darkpix}
\alias{darkpix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates percentage of dark pixels in an image.
}
\description{
This function takes a colour image file, typically a JPEG, and counts the percentage of pixels below a certain darkness threshold.  By default, only a central circle representing the hemispherical horizon is counted.  The threshold at which darkness is detected can be altered.
}
\usage{
darkpix(file, cutoff = 128, clip.circle = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Character string for the filename of the image file to be processed.
}
  \item{cutoff}{
The cutoff value, from 0-255.  Pixels with brighness values below this threshold are considered dark.
}
  \item{clip.circle}{
True/fale - if true, central circle of image is used. If false, whole image is used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A value representing the percentage of dark pixels in the image.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Grant Williamson <grant.williamson@utas.edu.au>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{darkpix.dir}}
}
\examples{
#rlogo=system.file("pictures/Rlogo.jpg", package="rgdal")[1]
#darkpix(rlogo,cutoff=100,clip.circle=T)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ hemispherical }
\keyword{ photo }
\keyword{ image }
\keyword{ cover }
