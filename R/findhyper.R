# last modification Friday, June 13
find.hyper <-function (model = NULL, 
                      parameters = NULL, 
                      other = NULL, 
                      penalty = 2,  
                      steps = c(0.1), 
                      lower = -Inf, 
                      upper = Inf, 
                      method="L-BFGS-B", 
                      ... )# 
{
  fn <- function(p,pen=penalty)
    {p<<-p
    if (!is.null(other)) eval(other)
      mod.1<-eval(call) 
      call<<-mod.1$call
     call$mu.start<<-fitted(mod.1,"mu")
     if ("sigma"%in%mod.1$parameters)  call$sigma.start <<- fitted(mod.1,"sigma")
     if (   "nu"%in%mod.1$parameters)  call$nu.start    <<- fitted(mod.1,"nu")
     if (  "tau"%in%mod.1$parameters)  call$tau.start   <<- fitted(mod.1,"tau")
     cat("par",p,"crit=",IC(mod.1,pen),"with pen=",pen,"\n")
      IC(mod.1,pen)
    }
 if(is.null(model)) stop("you have not defined the model")
 if(is.null(parameters)) stop("you have not define the parameters")
  lp <- length(parameters)
 if (lp==length(steps)) ndeps <- steps else ndeps <-rep(steps[1],lp)
   p <- parameters
  assign("p", parameters, envir=globalenv())
  rm(p)
  if (!is.null(model$data)) 
  { 
       attach(eval(model$data))
    #attach(eval(substitute(model$data))
    on.exit(detach(eval(model$data)))
  }
  if (!is.null(other)) eval(other)
  mod.1<-eval(model)
  call<-mod.1$call
  call$mu.start<-fitted(mod.1,"mu")
 if ("sigma"%in%mod.1$parameters)  call$sigma.start <- fitted(mod.1,"sigma")
 if (   "nu"%in%mod.1$parameters)  call$nu.start    <- fitted(mod.1,"nu")
 if (  "tau"%in%mod.1$parameters)  call$tau.start   <- fitted(mod.1,"tau")
  o2<-optim(parameters, fn, lower = lower, upper = upper, 
            method = method, control = list(ndeps = ndeps , ...))
  on.exit(rm(p,envir=.GlobalEnv),add=TRUE)
  o2
}
