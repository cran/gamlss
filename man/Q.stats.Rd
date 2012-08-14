\name{Q.stats}
\alias{Q.stats}

\title{A function to calculate the Q-statistics}
\description{
 This function calculates and prints the Q-statistics which are useful to test normality of the residuals within a range 
 of an independent variable, for example age in centile estimation, see Royston and Wright (2000).
}
\usage{
Q.stats(obj, xvar = NULL, xcut.points = NULL, n.inter = 10, zvals = TRUE, 
        save = TRUE, ...)
}

\arguments{
  \item{obj}{a GAMLSS object or any other residual vector}
  \item{xvar}{a unique explanatory variable}
  \item{xcut.points}{the x-axis cut off points e.g. \code{c(20,30)}. If \code{xcut.points=NULL} then the \code{n.inter} argument is activated }
  \item{n.inter}{if \code{xcut.points=NULL} this argument gives the number of intervals in which the x-variable will be split, with default 4}
  \item{zvals}{if \code{TRUE} the output matrix contains the individual z's rather that  Q statistics}
  \item{save}{whether to save the Q-statistics or not with default equal to \code{TRUE}. 
            In this case the functions produce a matrix giving individual Q (or z) statistics and the final aggregate Q's}
   \item{\dots}{for extra arguments}
}

\value{
  A matrix containing the individual Q's and the aggregate Q-statistics
}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Royston P. and Wright E. M. (2000) Goodness of fit statistics for the age-specific reference intervals. 
\emph{Statistics in Medicine}, 19, pp 2943-2962.  

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}, with contributions from Elaine Borghie}

\seealso{  \code{\link{gamlss}}, \code{\link{centiles.split}},  \code{\link{wp}} }
\examples{
data(abdom)
h<-gamlss(y~pb(x), sigma.formula=~pb(x), family=BCT, data=abdom) 
Q.stats(h,xvar=abdom$x,n.inter=8)
Q.stats(h,xvar=abdom$x,n.inter=8,zvals=FALSE)
rm(h)
}
\keyword{regression}% 
