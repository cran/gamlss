\name{rqres.plot}
\alias{rqres.plot}

\title{Plotting Randomized Quantile Residuals }
\description{
 This function plots QQ-plots of the normalized randomized quantile residuals (see Dunn and Smyth, 1996) for a model using a discrete GAMLSS family distribution.   
 
}
\usage{
rqres.plot(obj = NULL, howmany = 6, all = TRUE, save = FALSE, ...)
}

\arguments{
  \item{obj}{a fitted GAMLSS model object from a "discrete" type of family }
  \item{howmany}{The number of QQ-plots required up to ten i.e. \code{howmany=6}}
  \item{all}{if TRUE QQ-plots from \code{howmany} realizations are plotted. 
             If FALSE then a single qq-plot of the median of the \code{howmany} realizations is plotted}
  \item{save}{If TRUE the median residuals can be saved }
   \item{\dots}{for extra arguments}
}
\details{For discrete family distributions, the \code{\link{gamlss}()} function saves on exit one realization of randomized quantile residuals which 
can be plotted using the generic function \code{plot} which calls the \code{plot.gamlss}. Looking at only one realization can be misleading, so the 
current function creates QQ-plots for several 
realizations. The function allows up to 10 QQ-plots to be plotted. Occasionally  one wishes to create a lot of realizations 
and then take a median of them (separately for each ordered value) to create a single median realization. The option \code{all} in combinations 
with the option \code{howmany} creates a 
QQ-plot of the medians of the normalized randomized quantile residuals. These 'median' randomized quantile residuals can be saved using the option
(\code{save=TRUE}).  
 
}
\value{
  
  If \code{save} it is TRUE then the vector of the median residuals is saved. 
}
\references{ Dunn, P. K. and Smyth, G. K. (1996) Randomised quantile residuals,
    \emph{J. Comput. Graph. Statist.}, \bold{5}, 236--244
    
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}}



\section{Warning }{....}  
\seealso{ \code{\link{plot.gamlss}}, \code{\link{gamlss}}  }

\examples{
data(aids) # fitting a model from a discrete distribution 
h<-gamlss(y~pb(x)+qrt, family=NBI, data=aids) # 
plot(h)
# plot qq- plots from 6 realization of the randomized quantile residuals
rqres.plot(h) 
# a qq-plot from the medians from 40 realizations
rqres.plot(h,howmany=40,all=FALSE) # 
}
\keyword{regression}% 
