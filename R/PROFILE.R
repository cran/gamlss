prof.dev<- function ( object, 
                      which = NULL,  
                      min = NULL, 
                      max = NULL, 
                      step = NULL, 
                      startlastfit = TRUE, 
                      type = "o",
                      plot = TRUE,
                      perc = 95, 
                      ...) 
{
    cat("*******************************************************************", "\n")
   if (is.null(object$y) ) 
       stop(paste(deparse(substitute(object)), "does not have `y' components"))
   if (is.null(min)) 
       stop(paste("The min has not been specified"))
   if (is.null(max)) 
       stop(paste("The max has not been specified"))
   if (is.null(step)) 
       stop(paste("The step has not been specified"))
   if (is.null(which)) 
       stop(paste("The parameter has not been specified"))
   else 
    {
       if (!any(which == c("mu", "sigma", "nu", "tau")))  
           stop(paste("The parameter has not been specified correctly."))
    }     
    interval <- seq(from=min, to=max, by=step)
    interval <- ifelse(interval==0,0.0001,interval)
    fix.parameter <- rep(0, length(interval))
    G.deviances <- rep(0, length(interval))
    m <- length(interval)
    for (i in 1:length(interval))
    {
        if (which == "mu") 
        {   
            if ("mu"%in%object$parameters) 
            {
                mu.start <- rep(interval[i],length(object$y))
                object$call$mu.fix = TRUE
                object$call$mu.start <- mu.start
                cat("mu.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The mu is not a valid parameter of the current model."))
        }
        if (which == "sigma")
        {
            if ("sigma"%in%object$parameters)  
            {
                sigma.start <- rep(interval[i],length(object$y))
                object$call$sigma.fix = TRUE
                object$call$sigma.start <- sigma.start
                cat("sigma.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The sigma is not a valid parameter of the current model."))
        }
        if (which == "nu")
        {
            if ("nu"%in%object$parameters) 
            {
                nu.start <- rep(interval[i],length(object$y))
                object$call$nu.fix = TRUE
                object$call$nu.start <- nu.start
                cat("nu.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The nu is not a valid parameter of the current model."))
        }
        if (which == "tau") 
        {
            if ("tau"%in%object$parameters)
            {
                tau.start <- rep(interval[i],length(object$y))
                object$call$tau.fix = TRUE
                object$call$tau.start <- tau.start
                cat("tau.start=(" ,interval[i],")","\n")
            }
            else
                stop(paste("The tau is not a valid parameter of the current model."))
        }
        
        #                    FIT THE MODEL     
        object$call$control <- object$control
                        out <- eval(object$call)
           fix.parameter[i] <- interval[i]
             G.deviances[i] <- out$G.deviance
            if (startlastfit)
        {
            if ("mu"%in%object$parameters)       object$call$mu.start <- out$mu.fv
            if ("sigma"%in%object$parameters) object$call$sigma.start <- out$sigma.fv
            if ("nu"%in%object$parameters)       object$call$nu.start <- out$nu.fv
            if ("tau"%in%object$parameters)     object$call$tau.start <- out$tau.fv
        }
        cat("*******************************************************************", "\n")
    }
 
 prof.out <- cbind(interval, G.deviances)
if (plot)
  {   
      xl <- as.vector(interval)
  loglik <- G.deviances
  xlabel <- paste("Grid of the",which,"parameter")
    xlab <- xlabel
    ylab <- "Global Deviances"
    main <- "Profile Global Deviance"
    Gmin <- min(G.deviances)
    mx <- which.min(G.deviances)
    min.parameter <-fix.parameter[mx]
    lim <- Gmin + qchisq((perc/100), 1)
    op <- par(pin=c(5,4),col.axis="blue4",col.main="blue4",col.lab="blue4")
   # plot.new() 
    #plot(xl, loglik, type="n")
    plot(xl, loglik, xlab = xlab, ylab = ylab, main=main,
         col="darkgreen", frame.plot = TRUE,  type=type) # MS Tuesday, June 10, 2003 at 09:50
         # points(par(col="blue4")), taken out
    #lines(xl, loglik, col="red") # MS Tuesday, June 10, 2003 at 09:49
    par(op)
    plims <- par("usr")
    la <- xl[mx]
    y0 <- plims[3]
    if (lim < max(loglik))
    {
        abline(h = lim, lty = 3)
        y0 <- plims[3]
        scal <- (1/10 * (plims[4] - y0))/par("pin")[2] #par("pin")[2]=the height of the plot in inches
        scx <- (2/10 * (plims[2] - plims[1]))/par("pin")[1] #par("pin")[1]=the width of the plot in inches
        # MS change to 2/10 from 1/10 Sunday, December 9, 2007 at 23:29
        text(xl[1] + scx, lim + scal, paste(perc,"%"))
    }
    if (mx > 1 && mx < m) 
        segments(la, y0, la, Gmin, lty = 3)
    ind <- range((1:m)[loglik < lim]) #gets the x-values for the given range that their GD is in the range GD-or+3.84
    #Defines the lower bound
    if (loglik[1] > lim) {
        i <- ind[1]
        x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] - 
                xl[i - 1]))/(loglik[i] - loglik[i - 1])
        min.ci.par <- x
        segments(x, y0, x, lim, lty = 3)
    }
    #Defines the upper bound
    if (loglik[m] > lim) {
        i <- ind[2] + 1
        x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] - 
                xl[i - 1]))/(loglik[i] - loglik[i - 1])
        max.ci.par <- x
        segments(x, y0, x, lim, lty = 3)
    }
    cat("*******************************************************************", "\n")
    cat("Best estimate of the fixed parameter is " ,min.parameter, "\n")
    cat("with a Global Deviance equal to ", Gmin, " at position ", mx, "\n")
    if ((loglik[1] > lim) && (loglik[m] > lim))    
        {cat("A ", perc,"%  Confidence interval is: (" ,min.ci.par, ",", max.ci.par, ") \n")}
    cat("*******************************************************************", "\n")           
  return(invisible(prof.out))
  }
    else return(prof.out)



}
