\name{goats}
\alias{goats}
\docType{data}
\title{
Zoometric measurements of goats
}
\description{
Zoometric measurements of 27 week old creole goats collected by Dorantes-Coronado (2013).

}
\usage{data(goats)}
\format{
  A data frame with 52 rows and 7 columns containing measurements (in kilograms and centimeters) on the following variables.
  \describe{
    \item{\code{body.weight}}{ }
    \item{\code{body.length}}{ }
    \item{\code{trunk.length}}{}
    \item{\code{withers.height}}{}
    \item{\code{thoracic.perimeter}}{}
    \item{\code{hip.length}}{}
    \item{\code{ear.length}}{}
  }
}

\source{
Dorantes-Coronado (2013).
}
\references{
Dorantes-Coronado, E.J. (2013). Estudio preliminar para el establecimiento de un programa de mejoramiento genetico de cabras en el Estado de Mexico. Ph.D. Thesis. Colegio de Postgraduados, Mexico.
}
\examples{
data(goats)
plot(goats)
}
\keyword{datasets}
