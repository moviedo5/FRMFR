\name{Pork}
\alias{Pork}
\docType{data}
\title{Pork fat samples}
\description{The spectrum curves of Raman spectroscopy and the percentages of the unsaturated fatty
  acid obtained from 105 samples of the fat taken from the daily
  production stock of a slaughterhouse.}
\usage{data(Pork)}
\format{A list of two matrix:
 \describe{
  \item{X}{a matrix with 105 rows and 566 columns with each row
    representing the spectrum curve of Raman spectroscopy from one fat sample.}
  \item{Y}{a vector of length 105, giving the percentages of the
    unsaturated fatty acid of the 105 samples.}
  }
}
\keyword{datasets}
\references{
  Lotte Boge Lyndgaard, Klavs Martin Sorensen, Frans van den Berg and
  Soren Balling Engelsen (2012): Depth profiling of porcine adipose
  tissue by Raman spectroscopy, Journal of Raman Spectroscopy, 43,
  482-489.

  Xin Qi and Ruiyan Luo. (Manuscript) Functional regression for highly
  densely observed functional data with novel regularity.
  }
\source{ \url{http://www.models.life.ku.dk/RAMANporkfat}}
\seealso{ \code{\link{cv.sof.spike}}}
\examples{
 data(Pork)
 str(Pork)
}

