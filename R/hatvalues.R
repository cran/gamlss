######################################################################
######################################################################
# the hatvalues for GAMLSS
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
# the linear leverage of a gamlss object
hatvalues.gamlss <- function(model, ...)
{
#####################################################################
#####################################################################  
Formulae2one <- function(formula, sigma=~1, nu=~1, tau=~1, data )
  {
    form <- formula(formula)
    nform <- paste(paste(form[[2]],form[[1]]), deparse(form[[3]]), "+",  
                   deparse(sigma[[2]]),"+",
                   deparse(nu[[2]]),"+",
                   deparse(tau[[2]]))[1]
    ff<- formula(paste(nform, collapse = " "))
    environment(ff) <- globalenv()
    ff
  }
#####################################################################
#####################################################################
leverage <- function(formula = list(), 
                             data = NULL, 
                             weights = NULL, 
                             subset = NULL, 
                             na.action)
{
####################################################################  
# local
Formulae2data <- function(formula = list(), data=NULL, weights=NULL, subset=NULL, 
                            na.action, print = FALSE )
  {
    if (is(formula,"list"))
    {
      lenList <- length(formula)
      if (lenList==0) stop("no formula detected")
      if (lenList==1) 
      {
        ff <- deparse(formula[[1]])
      } else
      {
        # the first formula  
        form <- formula(formula[[1]])
        # create y~x+   
        f1 <- paste(paste(form[[2]],form[[1]]), deparse(form[[3]]), "+")
        # now add the of he formulae    
        for (i in 2:lenList)
        {
          ff <- if (i==lenList) paste(f1, deparse(formula[[i]][[2]]))
          else paste(f1, deparse(formula[[i]][[2]]),"+")
        } 
      }
    } else if (is(formula,"formula")) {ff  <- formula}
    else stop("The formula argument should be a formula or a list") 
    if (!is.null(weights)) 
    {
      # formula(paste(ff[[3]], collapse = " "))
      ff <- paste(ff, deparse(substitute(weights)), sep="+")
      # ff[[3]] <- paste(ff[[3]],deparse(substitute(weights)), sep="+")
    }
    environment(ff) <- globalenv()    # do I need this
    all.vars <- get_all_vars(ff, data=data)
    if (!is.null(data)&&class(data)!="data.frame") warning("data is not a data frame class attributes will be lost")
    M <- dim(all.vars)[1]
    ## subsetting             
    if (!is.null(subset)) {
      r <- if (!is.null(data))  eval(substitute(subset), data,  parent.frame())
      else eval(substitute(subset),  parent.frame())
      if (!is.logical(r)) stop("'subset' must be logical")
      all.vars <- all.vars[r,]
      M <- dim(all.vars)[1]
      if (print) cat( M, "observations left after subsetting \n" )           
    }
    # it need a futher warning here      N <- dim(all.vars)[1]  
    # na.omit   
    all.vars <- na.omit(all.vars)                             # clear NA's
    N <- dim(all.vars)[1]     
    if (print) {if (M-N > 0) cat(M-N, "rows with NAs are deleted", "\n" )}
    if (print) cat( N, "observations with", dim(all.vars)[2], "variables \n")    
    attr(all.vars, "formula") <- ff
    all.vars
  }
################################################################  
    d_f <- Formulae2data(formula=formula, data=data, weights=weights, subset=subset, 
                       na.action, print = FALSE)
  dimDF <- dim(d_f)
      p <-  dimDF[2]-1 # in general the number of x variables
      n <- dimDF[1]   # in general the number of cases, n=99 here
     m1 <- lm(attr(d_f, "formula"), data= d_f)
      h <- hatvalues(m1)
  return(h)
}
################################################################
################################################################
if (!is.null(model$call[["data"]]))
{
  DaTa <- eval(model$call[["data"]])
} else stop("the data argument is not used in the fitted model")
  lpar <- length(model$parameters)
  form <- switch(lpar, 
               Formulae2one(model$mu.formula),
               Formulae2one(model$mu.formula, model$sigma.formula),
               Formulae2one(model$mu.formula, model$sigma.formula, model$nu.formula),
               Formulae2one(model$mu.formula, model$sigma.formula, model$nu.formula, model$tau.formula)
               )
lev <- leverage(form, data=DaTa)
lev
}
##############################################################
##############################################################
##############################################################
##############################################################  
