\name{leafa}
\alias{leafa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates area in cm of leaves in an image.
}
\description{
This function takes a scanned colour image file, typically a JPEG, counts the dark pixels and converts them into an area in cm2.
}
\usage{
leafa(file, cutoff = 100, dpi=300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
Character string for the filename of the image file to be processed.
}
  \item{cutoff}{
The cutoff value, from 0-255.  Pixels with brighness values below this threshold are considered dark.
}
  \item{dpi}{
DPI (dots per inch) at which the image was scanned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An area, in cm, of dark pixels in the image.
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
rlogo=system.file("pictures/Rlogo.jpg", package="rgdal")[1]
darkpix(rlogo,cutoff=100,dpi=50)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ leaf }
\keyword{ photo }
\keyword{ image }
\keyword{ area }
