\name{corn}
\alias{corn}
\docType{data}
\title{NIR of corn samples}
\description{This data set consists of measurements on 80 corn samples. For each sample, its moisture, oil, protein
and starch values, together with its near-infrared (NIR) spectrum curve, were measured. The
wavelength range of the NIR curve is from 1100nm to 2498nm, and the curve was measured at every 2 nm.}
\usage{data(corn)}
\format{A list of 2 matrices:
 \describe{
  \item{X}{an 80*700 matrix with each row giving the NIR spectrum of a
    corn sample.}
  \item{Y}{an 80*4 matrix with the 4 columns giving the moisture, oil,
    protein and starch values, respectively.}
  }
}
\source{ \href{http://www.eigenvector.com/data/Corn/index.html}{Data Link}}
\seealso{ \code{\link{cv.msof}}}
\examples{
 data(corn)
 str(corn)
}

