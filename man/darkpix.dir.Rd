\name{darkpix.dir}
\alias{darkpix.dir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process hemispherical photo directory.
}
\description{
Processes a directory if image files, counts dark percentage and produces data table.
}
\usage{
darkpix.dir(directory, pattern = "jpg", cutoff = 100, clip.circle = T)
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
  \item{clip.circle}{
Clip processing to a central circle in the image, otherwise use entire image.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing the file name, and the percentage of dark pixels.
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
#darkpix.dir("fieldtrips/photo12/",pattern="jpg")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{hemispherical}
\keyword{photo}
\keyword{image}
\keyword{cover}