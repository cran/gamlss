\name{lms}
\alias{lms}

\title{A function to fit LMS curves for centile estimation}
\description{
This function is design to help the user to construct centile estimation.
It is only applicable when only "one" explanatory variable is available (usually age).
}
\usage{
lms(y, x, families = LMS, data = NULL, k = 2, 
     cent = 100*pnorm((-4:4)*2/3), 
     calibration = TRUE, legent = FALSE, mu.df = NULL, 
     sigma.df = NULL, nu.df = NULL, tau.df = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The response variable}
  \item{x}{The unique explanatory variable
}
  \item{families}{a list of \code{gamlss.families} with default \code{LMS=c("BCCGo",  "BCPEo", "BCTo")}
}
  \item{data}{the data frame}
  \item{k}{the penanlty to be used in the GAIC}
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{calibration}{whether calibrartion is required with default \code{TRUE}}
  \item{legent}{whether a legent is required with default \code{FALSE}}
  \item{mu.df}{\code{mu} effective degrees of freedon if required otherwise are estimated}
  \item{sigma.df}{\code{sigma} effective degrees of freedon if required otherwise are estimated}
  \item{nu.df}{\code{nu} effective degrees of freedon if required otherwise are estimated}
  \item{tau.df}{\code{tau} effective degrees of freedon if required otherwise are estimated}
  \item{\dots}{extra argument which cen be passed to \code{gamlss}}
}
\details{
This function should be used if only one explanatory variable is available.
The function starts by fitting a normal error distribution and a smooth function for \code{mu}. Then it proceeds by fitting several "appropriate" distributions for the response variable.
The set of \code{gamlss.family} distributions to fit is specified by the argument \code{families}. 
The default \code{families} arguments is \code{LMS=c("BCCGo",  "BCPEo", "BCTo")} that is the LMS class of distributions, Cole and Green (1992).
Note that this class is only appropriate when y is positive (with no zeros). If the response variable contains negative values and zeros then use the argument \code{families=theSHASH} where \code{theSHASH <-  c("NO", "SHASHo")} or add any other list of distributions which you may think is appropriate.
Justification of using the specific centile (0.38  2.27  9.1211220 25.25, 50, 74.75, 90.88, 97.72, 99.62) is given in Cole (1994).
}
\value{It returns a \code{gamlss} fitted object}

\references{
Cole, T. J. (1994) Do growth chart centiles need a face lift? \emph{BMJ}, 308--641.

Cole,  T. J. and Green,  P. J. (1992) Smoothing reference centile curves: the LMS method and penalized likelihood, \emph{Statist. Med.} \bold{11}, 1305--1319

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} and Vlasios Voudouris \email{v.voudouris@londonmet.ac.uk}
}
\note{The function is fitting several models and for large data can be slow}

\seealso{\code{\link{gamlss}}, \code{\link{centiles}}, \code{\link{calibration}}}
\examples{
\dontrun{
data(abdom)
m1<-lms(y,x , data=abdom, n.cyc=30)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
