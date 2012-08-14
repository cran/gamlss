\name{wp}
\alias{wp}
\title{ Worm plot }
\description{
Provides single or multiple worm plots for GAMLSS fitted objects. This is 
a diagnostic tool for checking the residuals within different ranges (by default not overlapping) of the explanatory variable
          }
\usage{
wp(object=NULL, xvar = NULL,  resid = NULL, n.inter = 4,
    xcut.points = NULL, 
   overlap = 0, xlim.all = 4, xlim.worm = 3.5, 
   show.given = TRUE, line = TRUE, 
   ylim.all = 12 * sqrt(1/length(resid)), 
   ylim.worm = 12 * sqrt(n.inter/length(resid)), 
   cex = 1, pch = 21, ...)
}

\arguments{
  \item{object}{a GAMLSS fitted object or any other fitted model where the resid() method works (preferably it should be quantile residuals)}
  \item{xvar}{ the explanatory variable against which the worm plots will be plotted }
  \item{resid}{if object is missing this argument can be used to specify the residual vector (again it should a quantile residuals or it be assumed to come from a normal distribution)}
  \item{n.inter}{the number of intervals in which the explanatory variable \code{xvar} will be cut}
  \item{xcut.points}{the x-axis cut off points e.g. \code{c(20,30)}. If \code{xcut.points=NULL} then the \code{n.inter} argument is activated }
  \item{overlap}{how much overlapping in the \code{xvar} intervals. Default value is \code{overlap=0} for non overlapping intervals}
  \item{xlim.all}{for the single plot, this value is the x-variable limit, default is \code{xlim.all=4}}
  \item{xlim.worm}{for multiple plots, this value is the x-variable limit, default is \code{xlim.worm=3.5}}
  \item{show.given}{whether to show the x-variable intervals in the top of the graph, default is \code{show.given=TRUE} }
  \item{line}{whether to plot the polynomial line in the worm plot, default value is \code{line=TRUE}}
  \item{ylim.all}{for the single plot, this value is the y-variable limit, default value is \code{ylim.all=12*sqrt(1/length(fitted(object)))}}
  \item{ylim.worm}{for multiple plots, this values is the y-variable limit, default value is \code{ylim.worm=12*sqrt(n.inter/length(fitted(object)))}}
  \item{cex}{ the cex plotting parameter with default \code{cex=1}}
  \item{pch}{ the pch plotting parameter with default \code{pch=21} }
   \item{\dots}{for extra arguments}
 }
\details{
 If the \code{xvar} argument is not specified then a single worm plot is used. In this case a worm plot is a detrended normal QQ-plot 
 so departure from normality is highlighted.
 If the  \code{xvar} is specified then we have as many worm plot as \code{n.iter}. In this case the x-variable is cut into \code{n.iter}
 intervals with an equal number observations and detrended normal QQ (i.e. worm) plots for each interval are plotted. 
 This is a way of highlighting failures of the model within different ranges of the 
 explanatory variable. The fitted coefficients from fitting cubic polynomials to the residuals (within each x-variable interval) can be obtain by e.g. 
 \code{coeffs<-wp(model1,xvar=x,n.iner=9)}.  van Buuren \emph{et al.} (2001) used these residuals to identify regions (intervals) of the 
 explanatory variable within which the model does not fit adequately the data (called "model violation")  
 
 
}
\value{
  For multiple plots the \code{xvar} intervals and the coefficients of the fitted cubic polynomials to the residuals (within each \code{xvar} interval) are returned.   
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), \emph{Appl. Statist.}, \bold{54}, part 3,
1-38. 
 
Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
           
van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for modelling growth reference curves. 
            \emph{Statistics in Medicine}, \bold{20}, 1259--1277
            }
\author{Mikis Stasinopoulos and Bob Rigby}


\seealso{  \code{\link{gamlss}}, \code{\link{plot.gamlss} }}

\examples{
data(abdom)
a<-gamlss(y~pb(x),sigma.fo=~pb(x,1),family=LO,data=abdom)
wp(a)
coeff1<-wp(a,abdom$x)
coeff1
rm(a,a1)
}
\keyword{regression}% 
