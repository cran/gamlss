# this is a general profile (deviance or GAIC) function 
# it can be used for profiling DF power transormation in x 
# or in fact any coeficient in a GAMLSS models
# created by MS Thursday, June 19, 2003 at 13:29 
# modified by MS Monday, August 25, 2003 at 13:25 to include intervals for profile GD 
prof.term <- function (model = NULL, 
                   criterion = "GD", 
                     penalty = 2.5, 
                       other = NULL,
                         min = NULL, 
                         max = NULL, 
                        step = NULL, 
                        type = "o", 
                      xlabel = NULL,
                        plot = TRUE,
                        term = TRUE,
                        perc = 95,
                        ...) 
{
if(is.null(model)) stop("you have not defined the model")
if(is.null(min)) stop("you have not defined the minimum")
if(is.null(max)) stop("you have not defined the maximum")
if(is.null(step)) stop("you have not defined the step")
if(!criterion%in%c("IC","GD")) stop("criterion should be IC or GD")
interval <- seq(from=min, to=max, by=step)
     I.C <- rep(0, length(interval)) 
    call <- model
 if (!is.null(model$data)) 
        {
        a<-model$data  
        attach(eval(substitute(a)))
        on.exit(detach(eval(substitute(a))))
        }
    for (i in 1:length(interval))
    {
    this<<- this <-interval[i] # mikis Thursday, March 27, 2008 
    if (!is.null(other)) eval(other)
      mod.1<-eval(call) 
      call<-mod.1$call
     call$mu.start<-fitted(mod.1,"mu")
     if ("sigma"%in%mod.1$parameters)  call$sigma.start <- fitted(mod.1,"sigma")
     if (   "nu"%in%mod.1$parameters)  call$nu.start    <- fitted(mod.1,"nu")
     if (  "tau"%in%mod.1$parameters)  call$tau.start   <- fitted(mod.1,"tau")
 I.C[i]<-  if(criterion=="GD")  deviance(mod.1) else  IC(mod.1,penalty)
    }

    xlab <- if(!is.null(xlabel)) xlabel else "parameter" 
    ylab <- if(criterion=="GD") "Global Deviances" else   paste("GAIC pen=",penalty)
    main <- if(criterion=="GD") "Profile Global Deviance" else  "Profile GAIC"
prof.out <- cbind(interval, I.C)
 if(criterion=="GD") 
    {
    Gmin <- min(I.C)
    mx <- which.min(I.C)
    min.parameter <-interval[mx]
    lim <- Gmin + qchisq((perc/100), 1)
    xl <- as.vector(interval)
     m <- length(interval)
    }
if (plot) 
   {
    #op <- par(pin=c(5,4),col.axis="blue4",col.main="blue4",col.lab="blue4") 
    # plot.new() 
    plot(interval, I.C, xlab=xlab, ylab=ylab, main=main, col="darkgreen", frame.plot = TRUE, type=type)
    #par(op)# points(par(col="blue4"))
    plims <- par("usr")
    if(criterion=="GD" & term==TRUE)
      { 
            if (lim < max(I.C))
               {
                 abline(h = lim, lty = 3)
                  y0 <- plims[3]
                  scal <- (1/10 * (plims[4] - y0))/par("pin")[2] #par("pin")[2]=the height of the plot in inches
                  scx <- (2/10 * (plims[2] - plims[1]))/par("pin")[1] #par("pin")[1]=the width of the plot in inches 
                  # MS change to 2/10, Sunday, December 9, 2007 at 23:28  
                  text(xl[1] + scx, lim + scal, paste(perc,"%") )
                  la <- xl[mx]
                if (mx > 1 && mx < m) 
                    segments(la, y0, la, Gmin, lty = 3)
                }
               # if (mx > 1 && mx < m) 
               #     segments(la, y0, la, Gmin, lty = 3)
                ind <- range((1:m)[I.C < lim]) #gets the x-values for the given range that their GD is in the range GD-or+3.84
                  #Defines the lower bound
             if (I.C[1] > lim) 
               {
                 i <- ind[1]
                 x <- xl[i - 1] + ((lim - I.C[i - 1]) * (xl[i] - 
                 xl[i - 1]))/(I.C[i] - I.C[i - 1])
                 min.ci.par <- x
                  segments(x, y0, x, lim, lty = 3)
                }
              #Defines the upper bound
            if (I.C[m] > lim) 
               {
                i <- ind[2] + 1
                x <- xl[i - 1] + ((lim - I.C[i - 1]) * (xl[i] - 
                xl[i - 1]))/(I.C[i] - I.C[i - 1])
                max.ci.par <- x
                segments(x, y0, x, lim, lty = 3)
                }
                cat("*******************************************************************", "\n")
                cat("Best estimate of the fixed parameter is " ,min.parameter, "\n")
                cat("with a Global Deviance equal to ", Gmin, " at position ", mx, "\n")
                if ((I.C[1] > lim) && (I.C[m] > lim))    
                {cat("A ", perc,"% Confidence interval is: (" ,min.ci.par, ",", max.ci.par, ") \n")}
                cat("*******************************************************************", "\n")             
                
      }           
            return(invisible(prof.out))
    }
    else {
           return(prof.out)
         }
}
