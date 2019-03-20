# library(gamlss)
# those are function for checking the tails of a distribution
# created by Mikis Stasinopoulos Bob Rigby and Vlasios Voudouris Dec 2011 
# rewriten MARCH 2018 by Mikis
# TO DO
# 
#-------------------------------------------------------------------------------
# log log Survival plots
# there are three functions for that plus one which connect them al
# i)   loglogSurv1()  : method 2 for for type I tails (log WEI)
# ii)  loglogSurv2()  : method 2 for type II tails (WEI)
# iii) loglogSurv3()  : method 2 for type III tails (GUMBEL)
# iv)  loglogSurv     : combines the above functions and selects the best
#-------------------------------------------------------------------------------
# log Survival plot
# v)    logSurv() : plots the empirical survival function log(1-ecdf) 
#        against log(y) for part of the data
# vi)  loglogplot() plots the empirical survival function log(1-ecdf) for all data
# v)  ECDF()
#-------------------------------------------------------------------------------
# those are in gamlss.tr
# vii)     fitTail : fits a truncated gamlss.family distribution to the tail of the data
# viii) fitTailAll : fits a (Hill type) series of Fit using the fitTail function 
#------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# TYPE I
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
loglogSurv1 <- function(y, 
                         prob = 0.9, 
                         print = TRUE,
                         title = NULL,
                         lcol = gray(.1),
                         ltype = 1,
                        ...)
{
  #-----------------
  Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
    cdF <- ECDF(y) # get ecdf for all data
     mY <- quantile(y, probs=prob)
      Y <- y[y>mY]
if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
   newY <- log(-log(1-cdF(Y)))
  place <- "bottomleft"
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
      x <- log(log(Y))
     m1 <- gamlss(newY ~ x, trace=FALSE)
    ess <- sum((newY-fitted(m1))^2)
   Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
    {
    title(paste("Log Log Survival plot (Type I) for", prob, 
                          "of the right tail of",  deparse(substitute(y))))
    } else title(title)
  invisible(m1)  
}

#--------------------------------------------------------------------------------
# Type II
#--------------------------------------------------------------------------------
loglogSurv2 <- function(y, 
                          prob = 0.9, 
                         print = TRUE,
                         title = NULL,
                          lcol = gray(.1),
                         ltype = 1,
                         ...)
{
  #-----------------
  Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
    cdF <- ECDF(y) # get ecdf for all data
     mY <- quantile(y, probs=prob)
      Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
   newY <- log(-log(1-cdF(Y)))
  place <- "bottomleft"
   Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
howmany <- length(Y)
      x <- log(Y)
     m1 <- gamlss(newY ~ x, trace=FALSE)
    ess <- sum((newY-fitted(m1))^2)
   Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
  {
    title(paste("Log Log Survival plot (Type II) for", prob, 
                "of the right tail of",  deparse(substitute(y))))
  } else title(title)
  invisible(m1)  
}
#--------------------------------------------------------------------------------
# Type III
#--------------------------------------------------------------------------------
loglogSurv3 <- function(y, 
                     prob = 0.9, 
                    print = TRUE,
                    title = NULL,
                     lcol = gray(.1),
                    ltype = 1,
                        ...)
{
  #-----------------
    Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
      cdF <- ECDF(y) # get ecdf for all data
       mY <- quantile(y, probs=prob)
        Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
     newY <- log(-log(1-cdF(Y)))
    place <- "bottomleft"
     Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
  howmany <- length(Y)
        x <- Y
       m1 <- gamlss(newY ~ x, trace=FALSE)
      ess <- sum((newY-fitted(m1))^2)
     Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")  
  if (print)
  {
    cat("coefficients",  coef(m1), "\n")
    cat("error sum of squares", ess, "\n")
  }
  plot(newY~x, xlab=Xlab1, ylab=Ylab, ...)
  lines(fitted(m1)~x, col=lcol, lty=ltype)  
  if (is.null(title))
  {
    title(paste("Log Log Survival plot (Type III) for", prob, 
                "of the right tail of",  deparse(substitute(y))))
  } else title(title)
  invisible(m1)  
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# Select the best from type I II and III
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
loglogSurv <- function(y, 
                 prob = 0.9, 
                print = TRUE,
                title = NULL,
                 lcol = gray(.1),
                ltype = 1,
                 plot = TRUE,
                      ...)
{
#-----------------
# body of the function starts here
#  Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
  cdF <- ECDF(y) # get ecdf for all data
   mY <- quantile(y, probs=prob)
    Y <- y[y>mY]
  if (mY<1) stop("For the method for Type I to work the minimum value of the tail must be greater than 1") 
 newY <- log(-log(1-cdF(Y)))
 Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
  howmany <- length(Y)
    x1 <- log(log(Y))
    m1 <- gamlss(newY ~ x1, trace=FALSE)
  ess1 <- sum((newY-fitted(m1))^2)
    x2 <- log(Y)
    m2 <- gamlss(newY ~ x2, trace=FALSE)
  ess2 <- sum((newY-fitted(m2))^2)
    x3 <- Y
    m3 <- gamlss(newY ~ x3, trace=FALSE)
  ess3 <- sum((newY-fitted(m3))^2) 
   ess <- c(ess1, ess2, ess3)
   num <- which.min(ess)
  matcoef <- rbind(coef(m1), coef(m2), coef(m3))
  matcoef <- cbind(matcoef, ess)
 dimnames(matcoef) <- list(c("type I", "type II", "type III"), c(" Intercept",
" slope", " Error SS"))
if  (print)
{
  cat("Linear regression coefficients", "\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  matk <- matcoef[,-3]
  matk[,1] <-exp(matk[,1])
  colnames(matk) <- c("k:2,4,6", "k:1,3,5")
  cat("Estimates for parameters k", "\n")
  printCoefmat(matk, digits = 3, signif.stars = TRUE) 
}
if (plot)
 {
  
  Ylab <- paste("log(-log( S(", deparse(substitute(y)), ")))", sep="")
  switch(num, 
         { Xlab1 <- paste("log(log(", deparse(substitute(y)), "))", sep="")
           plot(newY~x1, xlab=Xlab1, ylab=Ylab, ...)
           lines(fitted(m1)~x1, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type I) for", prob, 
                         "of the right tail of",  deparse(substitute(y))))
           } else title(title)
           },
         {
           Xlab1 <- paste("log(", deparse(substitute(y)), ")", sep="")
           plot(newY~x2, xlab=Xlab1, ylab=Ylab, ...)
           lines(fitted(m2)~x2, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type II) for", prob, 
                         "of the right tail of",  deparse(substitute(y))))
           } else title(title)
         },
         {
           Xlab1 <- paste(deparse(substitute(y)), sep="")
           plot(newY~x3, xlab=Xlab1, ylab=Ylab,...)
           lines(fitted(m3)~x3, col=lcol, lty=ltype)
           if (is.null(title))
           {
             title(paste("Log Log Survival plot (Type III) for", prob, 
                         "of the right tail of",  deparse(substitute(y))))  
           } else title(title)   
           
         }
         )   
 }
model <- switch(num, m1,m2,m3)   
return(model) 
}
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# this plots the empirical log(1-ecdf) agains log(y)
# The complementary cumulative distribution function (CCDF) plot
# fits also a linear and a quadraitic fit to the resulting plot yo help interpetation 
#----------------------------------------------------------------------------- 
logSurv <- function(y, 
                  prob = 0.9, 
                  tail = c("right", "left"), 
                  plot = TRUE, 
                 print = TRUE,
                 title = NULL,
                  lcol = c(gray(.1),gray(.2), gray(.3)), 
                 ltype = c(1,2,3),
...)
{
# body of the function starts here
#  require(gamlss)
  tail <- match.arg(tail)
  Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
   cdF <- ECDF(y) # get ecdf for all data
    mY <- quantile(y, probs=prob)
  if (tail=="right")
        {  Y <- y[y>mY]
        newY <-  log(1-cdF(Y))
       newY2 <- log(-log(1-cdF(Y)))
       place <- "bottomleft"
        Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")
     howmany <- length(Y)
        } else
        {
           Y <- y[y<mY]
        newY <-  log(cdF(Y))
       newY2 <- log(-log(cdF(Y)))
       place <- "topleft"
        Ylab <- paste("log(F(",deparse(substitute(y)), "))", sep="")
     howmany <- length(Y)
        }
# model fitting
     m1 <- gamlss(newY ~ log(Y), trace=FALSE) # fitting k1=1 model
     m2 <- gamlss(newY ~ log(Y)+I(log(Y)^2), trace=FALSE)  # fitting k1=2 model
     m3 <- gamlss(newY2 ~ log(Y), trace=FALSE) # fitting k
   fv3  <- -exp(fitted(m3))
  if (plot)
  {
    plot(newY[order(Y)]~log(Y)[order(Y)], xlab=Xlab, ylab=Ylab, ...)
    lines(fitted(m1)[order(Y)]~log(Y)[order(Y)], col=lcol[1], lty=ltype[1], lwd=2)
    lines(fitted(m2)[order(Y)]~log(Y)[order(Y)], col=lcol[2], lty=ltype[2] , lwd=2)
    lines(fv3[order(Y)]~log(Y)[order(Y)], col=lcol[3], lty=ltype[3], lwd=2 )
    if (is.null(title))
    {
      title(main=paste(paste(prob,"%",sep=""), "of the", tail, "tail", "of", 
                       deparse(substitute(y)),howmany,"obs"))
    } else title(title)
    
    legend(place, legend=c("linear", "quadratic", "exponential"), col=lcol, lty=ltype, lwd=2 )
  }
    M <-rbind(coef(m1), coef(m2)[-2], coef(m3))
    rownames(M) <- c("linear", "quadratic", "exponential")
    invisible(M) 
}
#-------------------------------------------------------------------------------
loglogplot <- function(y, nplus1=TRUE, ...)
{
  Xlab <- paste("log(", deparse(substitute(y)), ")", sep="")
  Ylab <- paste("log(1-F(",deparse(substitute(y)), "))", sep="")   
  if (any(y<=0)) {
     y <- y+abs(min(y))+1
    warning("negative values in y, it is shifted to y+abs(min(y))+1 ")
  }
     Y <- unique(sort(y))
     n <- length(y)
  ecdf <- if (nplus1)  cumsum(tabulate(match(y, Y)))/(n+1)
          else         cumsum(tabulate(match(y, Y)))/n
     y <- if (nplus1) 1-ecdf else  1-c(0, ecdf[-n])
  plot(Y, y, log="xy", xlab=Xlab,ylab=Ylab,...)
  invisible(data.frame(x=Y,y=y))
} 
#------------------------------------------------------------
loglogplot0 <- function(x,  ...)
{
  ff <- ecdf(x)
   x <- unique(sort(x))
   F <- ff(x)
  FF <- c(0, F[-length(F)])
   y <- 1-FF
 plot(x, y, log="xy",...)
   M <- data.frame(x,y)
  invisible(M)
} 
#-------------------------------------------------------------------------------
# a function for ecdf the the difference that it divide by n+1
ECDF <- function(y)
{
  ysort <- unique(sort(y))
      n <- length(y)
# ecdf0 <- cumsum(tabulate(match(y, ysort)))/n# standard definition
   ecdf <- cumsum(tabulate(match(y, ysort)))/(n+1)
  fun <- stepfun(ysort,c(0,ecdf))
  class(fun) <- c("ecdf", "stepfun")
  fun
}
#------------------------------------------------------------------------------