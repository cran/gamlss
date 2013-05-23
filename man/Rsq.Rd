\name{Rsq}
\alias{Rsq}

\title{Generalised R-squared for GAMLSS models}
\description{
This function gives the generalised R-squared of Nagelkerke (1991) for a GAMLSS model
}
\usage{
Rsq(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a GAMLSS object}
}
\details{
This uses the generalised definition for R-squared defined as 
\deqn{R^2=1- \left(\frac{L(0)}{L(\hat{\theta})}\right)^(2/n)}{R-squared =1-(L(null)/L(fitted))^(2/n)}
 where \eqn{L(0)} is the null model (only a constant is  fitted to all parameters) and 
 \eqn{L(\hat{\theta})}{L(fitted)} is the current fitted model.
}
\value{A single value}
\references{
Nagelkerke, N. J. (1991). A note on a general definition of the coefficient of determination. Biometrika, 78(3), 691-692.
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} }
\note{The null model is fitted using the function gamlssML() which can create warning messages}

\seealso{\code{\link{GAIC}}}
\examples{
data(aids)
m1 <- gamlss(y~x+qrt, data=aids, family=NBI)
Rsq(m1)
rm(m1)
}

\keyword{regression}

