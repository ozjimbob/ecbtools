\name{leafa.dir}
\alias{leafa.dir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process leaf area image scan directory.
}
\description{
Processes a directory of scanned image files, calculates area of dark pixels in each and outputs a table.
}
\usage{
leafa.dir(directory, pattern = "jpg", cutoff = 100, dpi=300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{
String with path to directory to be processed.
}
  \item{pattern}{
A pattern, suitable for the \code{\link{list.files}} function, indicating which image files to load. "jpg" is the most common, but "JPG" or "JPEG" may be produed by some cameras.
}
  \item{cutoff}{
Cutoff brighness level, 0-255.
}
  \item{dpi}{
DPI (dots per inch) of scanned images.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing the file name, and area of dark pixels in cm2.
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
\code{\link{darkpix}}
}
\examples{
#leafa.dir("fieldtrips/photo12/",pattern="jpg")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ leaf }
\keyword{ photo }
\keyword{ image }
\keyword{ area }
