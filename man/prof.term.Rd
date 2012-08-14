\name{prof.term}
\alias{prof.term}

\title{Plotting the Profile: deviance or information criterion for one of the terms (or hyper-parameters) in a GAMLSS model }
\description{
This functions plots the profile deviance for a chosen parameter included in the linear predictor of any of the \code{mu}, 
\code{sigma}, \code{nu} or \code{tau} models so profile confidence intervals can be obtained. In can also be used to plot the profile 
of a specified information criterion for any 
hyperparameter.  
}
\usage{
prof.term(model = NULL, criterion = "GD", penalty = 2.5, other = NULL, 
         min = NULL, max = NULL, step = NULL, type = "o", xlabel = NULL, 
         plot = TRUE, term = TRUE, perc = 95, ... )
}

\arguments{
  \item{model}{this is a GAMLSS model, e.g. \cr
                     \code{model=gamlss(y~cs(x,df=this),sigma.fo=~cs(x,df=3),data=abdom)}, where \code{this} indicates the (hyper)parameter to be profiled}
  \item{criterion}{whether global deviance ("GD") or information criterion ("IC") is profiled. The default is \code{criterion="GD"}}
  \item{penalty}{ The penalty value if information criterion is used in \code{criterion}, default \code{penalty=2.5}  }
  \item{other}{this can be used to evaluate an expression before the actual fitting of the model}
  \item{min}{the minimum value for the parameter e.g. \code{min=1} }
  \item{max}{the maximum value for the parameter e.g.  \code{max=20}}
  \item{step}{how often to evaluate the global deviance (defines the step length of the grid for the parameter) e.g. \code{step=1}}
  \item{type}{what type of plot required. This is the same as in \code{type} for \code{plot}, default value is \code{type="o"}, that is, both line and points}
  \item{xlabel}{if a label for the axis is required }
  \item{plot}{whether to plot, \code{plot=TRUE} or save the results, \code{plot=FALSE}  }
  \item{term}{this has the value TRUE and it should be changed to FALSE if a profile global deviance is required for a hyperparameter so the IC are suppressed}
  \item{perc}{what \% confidence interval is required}
  \item{\dots}{for extra arguments}
}
\details{
 This function can be use to provide likelihood based confidence intervals for a parameter involved in terms in the linear predictor(s). 
 These confidence intervals are more accurate than the ones obtained from the parameters' standard errors. 
 The function can also be used to plot a profile information criterion (with a given penalty) against a hyperparameter. This can be used to check  
 the uniqueness in hyperparameter determination using for example \code{find.df}. 
}
\value{
  Return a profile plot (if the argument \code{plot=TRUE}) or the values of the parameters and the IC or GD values otherwise 
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{ Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk} and Bob Rigby \email{r.rigby@londonmet.ac.uk} }
\section{Warning }{A dense grid (i.e. small step) evaluation of the global deviance can take a long time, so start with a sparse grid (i.e. large step) 
and decrease  gradually the step length for more accuracy. }  
\seealso{\code{\link{gamlss}},  \code{\link{prof.dev}} }
\examples{
data(aids)
gamlss(y~x+qrt,family=NBI,data=aids)
mod<-quote(gamlss(y ~ offset(this * x) + qrt, data = aids, family = NBI))
prof.term(mod, min=0.06, max=0.11, step=0.001)
mod1<-quote(gamlss(y ~ cs(x,df=this) + qrt, data = aids, family = NBI))
prof.term(mod1, min=1, max=15, step=1, criterion="IC")
mod2 <- quote(gamlss(y ~ x+I((x>this)*(x-this))+qrt,family=NBI,data=aids))
prof.term(mod2, min=1, max=45, step=1, criterion="GD")
rm(mod,mod1,mod2)
}
\keyword{regression}% 
