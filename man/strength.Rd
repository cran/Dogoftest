\name{strength}
\alias{strength}
\docType{data}
\title{
Compressive strength of maize seeds
}
\description{
Compressive strength and strain of maize seeds.
}
\usage{data("strength")}
\format{
  A data frame with 90 observations on the following 2 variables.
  \describe{
    \item{\code{strain}}{a numeric vector giving the relative change in length under compression stress  in millimeters.}
    \item{\code{cstrength}}{a numeric vector giving the compressive strength in Newtons.}
  }
}
\details{
These data correspond to maize seeds with floury endosperm and 8\% of moisture.  
}
\source{
Mancera-Rico, A. (2014).
}
\references{
Mancera-Rico, A. (2014). Contenido de humedad y tipo de endospermo en la resistencia a compresion en semillas de maiz. Ph.D. Thesis. Colegio de Postgraduados, Mexico.
}
\examples{
data(strength)
plot(strength)     # plot of "strain" versus "cstrength"
}
\keyword{datasets}
