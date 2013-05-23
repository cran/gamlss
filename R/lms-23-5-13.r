# Authors Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris
# created 11-04-12
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#This function is design to help the user to construct centile estimation.
#It is only applicable when only "one" explanatory variable is available (usually age).
#It stats by fitting a normal error distribution and smooth function for mu and then proceeds by fitting #several appropriate distributions.
#The set of gamlss.family distribution to fit are specified in the argument families. 
#The default families arguments is LMS=c("BCCGo",  "BCPEo", "BCTo") that is the LMS class of distributions.
#Note that this class is only appropriate when y is positive (with no zeros). If the response variable contains negative values and/or zeros then use 
#the argument theSHASH theSHASH <-  c("NO", "SHASHo") or add any other distribution which you think is appropriate
LMS <- c("BCCGo",  "BCPEo", "BCTo")
theSHASH <-  c("NO", "SHASHo")
#------------
lms <- function(y, x,
        families = LMS,
            data = NULL, 
               k = 2, # for the AIC
            cent = 100*pnorm((-4:4)*2/3),
     calibration = TRUE,
          legend = FALSE,
           mu.df = NULL,
        sigma.df = NULL,
           nu.df = NULL,
          tau.df = NULL,
          method.pb = c("ML", "GAIC"),
              ... 
                )
{
  require(gamlss)
  # the families to fit
      FAM <- families      
  # get the data
if (!is.null(data)) 
  {
  DaTa <- with(data, data.frame(y=y, x=x))
#  y <- with(data, y)
  x.x <- with(data, x)
  }
if (is.null(data))  
  {
  DaTa <- data.frame(y=y, x=x)
  }
  method.pb <- match.arg(method.pb)
  # starting value model (we assume that this will work). Note no sigma is fitted here
    switch(method.pb, 
         "ML"= {m0 <- gamlss(y~pb(x), sigma.formula=~1, data=DaTa, c.crit = 0.01)},
       "GAIC"= {m0 <- gamlss(y~pb(x, method="GAIC", k=k), sigma.formula=~1, data=DaTa, c.crit = 0.01)}
          ) # initial fit
  # in an ideal situation we should the fit tail function to decide which distrbutions we should fit here 
  # here we use the given list          
    failed <- list() 
      fits <- list()
       aic <- AIC(m0, k=k)
      fits <- c(fits, aic) 
  for (i in 1:length(FAM)) 
  {
    cat('*** Fitting', FAM[i], "***","\n")  
     switch(method.pb, 
         "ML"= { m1 <- try(gamlss(y ~ pb(x, df=mu.df),
              sigma.formula = ~pb(x, df=sigma.df),
                 nu.formula = ~pb(x, df=nu.df), 
                tau.formula = ~pb(x, df=tau.df), 
                family = FAM[i], data = DaTa,
               mu.start = fitted(m0), ...), silent=TRUE)},
        "GAIC"= { m1 <- try(gamlss(y ~ pb(x,  method="GAIC", k=k, df=mu.df),
              sigma.formula = ~pb(x,  method="GAIC", k=k, df=sigma.df),
                 nu.formula = ~pb(x,  method="GAIC", k=k, df=nu.df), 
                tau.formula = ~pb(x,  method="GAIC", k=k, df=tau.df), 
                family = FAM[i], data = DaTa,
               mu.start = fitted(m0), ...), silent=TRUE)
         })      
    if (any(class(m1)%in%"try-error")) # if fitting failed
    {
      cat(FAM[i], " failed", "\n")
          failed <- c(failed, FAM[i]) 
    }
    else
    {
             aic <- AIC(m1, k=k)
      names(aic) <- FAM[i]
            fits <- c(fits, aic)
      if (AIC(m1, k=k) < AIC(m0, k=k)) 
      {
        m0<-m1 
      }
    }
  }
  m0$failed <- failed
     fits <- unlist(fits)
  m0$fits <- fits[order(fits)] 
  m0$xvar <- with(DaTa,x)
     m0$y <- with(DaTa,y)
  if (!is.null(data)) m0$call$data  <- substitute(data)
  #if (m0$family[1]=="SHASH")  m0$sigma.fv <- m0$sigma.fv*MaxOilProduction
  # calibration
  if (calibration)
  {
    calibration(m0, xvar=x.x, cent=cent, pch = 15, cex = 0.5, col = gray(0.7), ylab=deparse(substitute(y)), xlab=deparse(substitute(x)), legend=legend)	
  } 
  else 
  {
    centiles(m0, xvar=x.x, cent=cent, pch = 15, cex = 0.5, 
             col = gray(0.7), ylab=deparse(substitute(y)), xlab=deparse(substitute(x)), legend=legend)		
  }
  m0  # save the last model
}
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# this function is appropriate to used when fitted model fails to c
calibration <- function(object, xvar, cent=100*pnorm((-4:4)*2/3), legend=FALSE, fan=FALSE,  ...)
{
  z   <-  quantile(resid(object), probs = cent/100)
  p   <-  pNO(z, mu=0, sigma=1)
  percent <- 100*p
  if (fan)
  {
    centiles.fan(object, xvar=xvar, cent=percent,   ...)  
  }
  else
  {
    centiles(object, xvar=xvar, cent=percent, legend=legend,  ...)
  }
}
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------