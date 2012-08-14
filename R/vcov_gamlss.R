# the variance covariance function of GAMLSS
# rewriten for R 2.4.0 
vcov.gamlss <- function (object, type=c("vcov", "cor", "se", "coef", "all"), ...)
{  
#----------------------------------------------------------------------------------------
#****************************************************************************************
#----------------------------------------------------------------------------------------
##  the proper function here 
   type <- match.arg(type)
   if (!is.gamlss(object))  stop(paste("This is not an gamlss object", "\n", ""))
      call <- object$call
   if (is.null(call))   stop("need an object with call component")
 #  if  (grep("$", deparse(call)))  warning("need an object with no call component containing $")
        # the first thing is to get the formula from gamlss
   formula <-  call$formula
        # rename the call
 call[[1]] <- as.name("gamlssNonLinear")
        # get the y
      yvar <- formula[[2]]
        # delete the formula from call
        call[["formula"]] <- NULL
        # make sure that the gamlssNonLinear control is one iteration
      cont <- NonLinear.control(iterlim = 1)
        # create the list with the extra arguments
    extras <- list(y=yvar, mu.formula=formula, control=cont)
        # from  "mu to tau"
        coefBeta<-list()
        for (i in object$par)
         {
         if (i=="mu") 
           {
           if (!is.null(unlist(attr(terms(formula(object),specials=gamlss:::.gamlss.sm.list), "specials"))))# gamlss:::
                warning("addive terms exists in the mu formula standard errors for the linear terms maybe are not appropriate")
           }
          else
           {
            if (!is.null(unlist(attr(terms(formula(object,i),specials =gamlss:::.gamlss.sm.list), "specials")))) #gamlss:::
                warning(paste("addive terms exists in the ", i, " formula standard errors for the linear terms maybe are not appropriate"))
           } 
         parname <- paste(i,"start",sep=".")
      extras$ppp <- coef(object,i)
        names(extras)[[length(extras)]] <- parname
        coefBeta <- c(coefBeta,coef(object,i))
         }  
 existing <- !is.na(match(names(extras), names(call)))
   for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
       if (any(!existing)) 
          {
           call <- c(as.list(call), extras[!existing])
           call <- as.call(call)
          }
   suppressWarnings(        
       a <- eval(call)
                   )
           rownames(a$cov) <- attr(coefBeta,"names")
           colnames(a$cov) <- attr(coefBeta,"names")
           rownames(a$cor) <- attr(coefBeta,"names")
           colnames(a$cor) <- attr(coefBeta,"names")  
               names(a$se) <- attr(coefBeta,"names")                 
                  coefBeta <- unlist(coefBeta)
           names(coefBeta) <- attr(a$se, "names")
  switch(type,"vcov" = vcov(a),
               "cor" = a$cor, 
                "se" = a$se,
              "coef" = coefBeta,
               "all" = list(coef=coefBeta, se=a$se, vcov=vcov(a), cor=a$cor))
}
#****************************************************************************************
inteprFormula <- function (.z, ...) 
UseMethod("inteprFormula")
#----------------------------------------------------------------------------------------
inteprFormula.default <- function (.z, .envir = parent.frame(), .formula = FALSE, 
                               .vector = TRUE, .args = NULL, .start = 1, 
                               .name = NULL, .expand = TRUE, .intercept = TRUE, 
                               .old = NULL, .response = FALSE, ...) 
{
    if (!inherits(.z, "formula")) 
        return(NULL)
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    if (!is.environment(.envir)) {
        if (is.null(.name)) 
            .name <- paste(deparse(substitute(.envir)))
    #    if (inherits(.envir, "repeated")) 
    #        return(inteprFormula.repeated(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old, 
    #            .response))
    #    if (inherits(.envir, "tccov")) 
    #        return(inteprFormula.tccov(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old))
    #    if (inherits(.envir, "tvcov")) 
    #        return(inteprFormula.tvcov(.z, .envir, .formula, .vector, 
    #            .args, .start, .name, .expand, .intercept, .old))
        if (inherits(.envir, "data.frame")) 
            return(inteprFormula.data.frame(.z, .envir, .formula, .vector, 
                .args, .start, .name, .expand, .intercept, .old))
    }
    .pars <- .range <- NULL
    if (!is.null(.old)) {
        if (!is.list(.old)) 
            .old <- list(.old)
        for (.j in .old) {
            if (!inherits(.j, "formulafn")) 
                stop("objects in .old must have class, formulafn")
            .pars <- c(.pars, attr(.j, "parameters"))
            .range <- c(.range, attr(.j, "range")[1]:attr(.j, 
                "range")[2])
        }
        if (.start <= max(.range)) 
            warning("possible conflict in vector indexing - check .start")
    }
    if (!is.null(.args) && !is.character(.args)) 
        stop(".args must be a character string")
    .zz <- funMobJ(.z)
    .ch <- .zz$formula
    .mem <- .zz$objects
    .fcn <- .zz$functions
    if ("$" %in% .fcn) 
        stop("sublists not allowed (attach dataframes and use variable names)")
    .ex <- .zz$covariates
    .fac <- .zz$factors
    .local <- .zz$local
    rm(.zz)
    if (length(.mem) > 0) {
        .un <- unique(.mem[!.ex & !.fac & !.local])
        if (length(unique(.mem[.ex | .fac | .local])) == 0 && 
            length(.un) == 0) 
            warning("inteprFormula.default: no variables found")
    }
    if (length(.mem) == 0 || all(.ex | .fac | .local)) {
        if (.formula) 
            return(.z)
        else {
            if (any("offset" %in% .fcn)) 
                stop("offset not allowed")
            .mt <- terms(.z)
            if (is.numeric(.mt[[2]])) {
                .dm <- matrix(1)
                colnames(.dm) <- "(Intercept)"
            }
            else {
                .dm <- model.matrix(.mt, model.frame(.mt, .envir))
                if (!.intercept) 
                  .dm <- .dm[, -1, drop = FALSE]
            }
            .fna <- function(.p) as.vector(.dm %*% .p[attr(.fna, 
                "range")[1]:attr(.fna, "range")[2]])
            attributes(.fna) <- list(formula = .z, model = colnames(.dm), 
                covariates = if (length(.mem) > 0) unique(.mem[.ex | 
                  .fac]) else NULL, parameters = paste("p[", 
                  1:dim(.dm)[2], "]", sep = ""), range = c(.start, 
                  .start + dim(.dm)[2] - 1), class = "formulafn")
            .obj <- ls(all.names = TRUE)
            rm(list = .obj[.obj != ".fna" & .obj != ".dm"])
            rm(.obj)
            return(.fna)
        }
    }
    if (!is.null(.fac) && any(.fac)) 
        stop(paste("covariates in formulae with unknowns must not be factors\ncheck", 
            .mem[.fac]))
    .fna <- function(.p) eval(attr(.fna, "model"))
    if (.vector) {
        if (!is.null(.args)) {
            .tmp <- match(.args, .un)
            if (all(!is.na(.tmp))) 
                .un <- .un[-.tmp]
            .par <- "alist(.p="
            for (.j in 1:length(.args)) {
                .par <- paste(.par, ",", collapse = "")
                .par <- paste(.par, .args[.j], "=", collapse = "")
            }
            .par <- paste(.par, ")", collapse = "")
            formals(.fna) <- eval(parse(text = .par))
        }
        if (!is.null(.old)) {
            .j <- match(.pars, .un)
            .un <- .un[-.j]
            .pars <- .pars[!is.na(.j)]
            .range <- .range[!is.na(.j)]
            for (.j in 1:length(.pars)) .ch <- gsub(paste(" ", 
                .pars[.j], " ", sep = ""), paste(" .p[", .range[.j], 
                "] ", sep = ""), .ch)
        }
        if (length(.un) > 0) 
            for (.j in 1:length(.un)) .ch <- gsub(paste(" ", 
                .un[.j], " ", sep = ""), paste(" .p[", .start + 
                .j - 1, "] ", sep = ""), .ch)
    }
    else {
        .par <- "alist("
        for (.j in 1:length(.un)) {
            if (.j > 1) 
                .par <- paste(.par, ",", collapse = "")
            .par <- paste(.par, .un[.j], "=", collapse = "")
        }
        .par <- paste(.par, ")", collapse = "")
        formals(.fna) <- eval(parse(text = .par))
    }
    attributes(.fna) <- list(formula = .z, model = parse(text = .ch), 
        parameters = .un, common = .pars, covariates = unique(.mem[.ex]), 
        range = c(.start, .start + length(.un) - 1), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fna"])
    rm(.obj)
    return(.fna)
}
#----------------------------------------------------------------------------------------
inteprFormula.data.frame <- function (.z, .envir = NULL, .formula = FALSE, .vector = TRUE, 
    .args = NULL, .start = 1, .name = NULL, .expand = NULL, .intercept = TRUE, 
    .old = NULL, ...) 
{
    if (!inherits(.z, "formula")) 
        return(NULL)
    .pars <- .range <- NULL
    if (!is.null(.old)) {
        if (!is.list(.old)) 
            .old <- list(.old)
        for (.j in .old) {
            if (!inherits(.j, "formulafn")) 
                stop("objects in .old must have class, formulafn")
            .pars <- c(.pars, attr(.j, "parameters"))
            .range <- c(.range, attr(.j, "range")[1]:attr(.j, 
                "range")[2])
        }
        if (.start <= max(.range)) 
            warning("possible conflict in vector indexing - check .start")
    }
    if (!is.null(.args) && !is.character(.args)) 
        stop(".args must be a character string")
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    .ndata <- if (is.null(.name)) 
        paste(deparse(substitute(.envir)))
    else .name
    .cn <- colnames(.envir)
    .ex1 <- NULL
    .zz <- funMobJ(.z)
    .ch <- .zz$formula
    .mem <- .zz$objects
    .fcn <- .zz$functions
    .local <- .zz$local
    rm(.zz)
    if (length(.mem) > 0) {
        .ex1 <- match(.mem, .cn)
        .un <- unique(.mem[is.na(.ex1) & !.local])
        if (length(unique(.mem[!is.na(.ex1)])) == 0 && length(.un) == 
            0) 
            warning("inteprFormula.data.frame: no variables found")
    }
    .ex1a <- if (is.null(.ex1)) 
        NULL
    else .ex1[!is.na(.ex1)]
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) .ch <- gsub(paste(" ", .cn[.ex1a[.j]], 
            " ", sep = ""), paste(" ", .ndata, "$", .cn[.ex1a[.j]], 
            sep = ""), .ch)
    if (is.null(.ex1) || all(!is.na(.ex1))) {
        if (.formula) 
            return(.z)
        else {
            if (any("offset" %in% .fcn)) 
                stop("offset not allowed")
            .ch <- as.formula(paste("~", .ch))
            .mt <- terms(.ch)
            if (is.numeric(.mt[[2]])) {
                if (!.intercept) 
                  return(NULL)
                .n <- dim(.envir)[1]
                .dm <- matrix(1)
                colnames(.dm) <- "(Intercept)"
                .fna <- function(.p) rep(.p[attr(.fna, "range")[1]], 
                  .n)
            }
            else {
                .dm <- model.matrix(.mt, model.frame(.mt, data = .envir))
                if (!.intercept) 
                  .dm <- .dm[, -1, drop = FALSE]
                .fna <- function(.p) as.vector(.dm %*% .p[attr(.fna, 
                  "range")[1]:attr(.fna, "range")[2]])
            }
            attributes(.fna) <- list(formula = .z, model = colnames(.dm), 
                covariates = if (length(.mem) > 0) unique(.mem[!is.na(.ex1)]) else NULL, 
                parameters = paste("p[", 1:dim(.dm)[2], "]", 
                  sep = ""), range = c(.start, .start + dim(.dm)[2] - 
                  1), class = "formulafn")
            .obj <- ls(all.names = TRUE)
            rm(list = .obj[.obj != ".i" & .obj != ".fna" & .obj != 
                ".dm" & .obj != ".n"])
            rm(.obj)
            return(.fna)
        }
    }
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) if (is.factor(.envir[, .ex1a[.j]])) 
            stop(paste(colnames(.envir)[.ex1a[.j]], "is a factor variable"))
    .fna <- function(.p) eval(attr(.fna, "model"))
    if (!is.null(.args)) {
        .tmp <- match(.args, .un)
        if (all(!is.na(.tmp))) 
            .un <- .un[-.tmp]
        .par <- "alist(.p="
        for (.j in 1:length(.args)) {
            .par <- paste(.par, ",", collapse = "")
            .par <- paste(.par, .args[.j], "=", collapse = "")
        }
        .par <- paste(.par, ")", collapse = "")
        formals(.fna) <- eval(parse(text = .par))
    }
    if (!is.null(.old)) {
        .j <- match(.pars, .un)
        .un <- .un[-.j]
        .pars <- .pars[!is.na(.j)]
        .range <- .range[!is.na(.j)]
        for (.j in 1:length(.pars)) .ch <- gsub(paste(" ", .pars[.j], 
            " ", sep = ""), paste(" .p[", .range[.j], "] ", sep = ""), 
            .ch)
    }
    if (length(.un) > 0) 
        for (.j in 1:length(.un)) .ch <- gsub(paste(" ", .un[.j], 
            " ", sep = ""), paste(" .p[", .start + .j - 1, "] ", 
            sep = ""), .ch)
    attributes(.fna) <- list(formula = .z, model = parse(text = .ch), 
        parameters = .un, common = .pars, covariates = unique(.mem[!is.na(.ex1)]), 
        range = c(.start, .start + length(.un) - 1), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fna"])
    rm(.obj)
    return(.fna)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
"NonLinear.control" <- function(fscale = 1,   
                      typsize = NULL, #abs(p0), 
                      stepmax = NULL, # sqrt(p0 %*% p0),  
                      iterlim = 100,  
                       ndigit = 10,   
                      steptol = 1e-05, 
                      gradtol = 1e-05,
                  print.level = 0, 
            check.analyticals = TRUE,
                      hessian = TRUE )
{
##  Control iteration for GLIM
##  MS  Sunday, February 17, 2002 at 19:18
##
        if(fscale <= 0) 
        {
warning("the scale value supplied is zero or negative the default value of 1 was used instead")
          fscale <- 1
        }
 #       if(any(typsize <= 0)) 
 #       {
#warning("the value of typsize supplied is zero or negative the default value of abs(p0) was used instead")
#         typsize <- abs(p0)
#        }
#        if(stepmax < 0) 
#        {
#warning("the value of stepmax supplied is zero or negative the default value of  sqrt(p0 %*% p0) was used instead")
#         stepmax <-  sqrt(p0 %*% p0)
#        }
        if( iterlim <= 0) 
        {
warning("the value of  iterlim supplied is zero or negative the default value of 100 was used instead")
          iterlim <- 100
        } 
         if( ndigit < 0) 
        {
warning("the value of  ndigit supplied is zero or negative the default value of 10 was used instead")
           ndigit <- 10
        }   
         if(  steptol < 0) 
        {
warning("the value of  steptol supplied is zero or negative the default value of  1e-05 was used instead")
          steptol <- 1e-05
       }
         if(  gradtol < 0) 
       { 
warning("the value of  gradtol supplied is zero or negative the default value of  1e-05 was used instead")
          gradtol <- 1e-05
       } 
         if( print.level < 0) 
       {
 warning("the value of  print.level supplied is negative the default value of  0 was used instead")
      print.level <- 0
       }           
       list(fscale = fscale, typsize = typsize, stepmax = stepmax,  
            iterlim = iterlim, 
            ndigit = ndigit, steptol = steptol, gradtol = gradtol,  
            print.level =  print.level,
            check.analyticals = as.logical(check.analyticals)[1],
            hessian = as.logical(hessian)[1])
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
gamlssNonLinear <- function (
                y = NULL, 
       mu.formula = ~1, 
    sigma.formula = ~1, 
       nu.formula = ~1, 
      tau.formula = ~1,
           mu.fix = FALSE, 
        sigma.fix = FALSE, 
           nu.fix = FALSE, 
          tau.fix = FALSE, 
          all.fix = FALSE, 
         mu.start = NULL, 
      sigma.start = NULL, 
         nu.start = NULL, 
        tau.start = NULL, 
           family = NO(),
          weights = 1, 
            exact = FALSE, 
            delta = 1, 
             data = parent.frame(),
  #       typsize = abs(p0),
  #        stepmax = sqrt(p0 %*% p0),
          control = NonLinear.control(...),         
      llik.output = FALSE,
                 ... )
{
.gamlss.multin.list<-c("MULTIN", "MN3", "MN4", "MN5")
#----------------------------------------------------------------------------------------
rqres <- function (pfun = "pNO", 
                   type = c("Continuous", "Discrete", "Mixed"),
               censored = NULL,  
                   ymin = NULL, 
                 mass.p = NULL, 
                prob.mp = NULL,
                      y = y,
                         ... )
{ }
body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
# this function is needed to get the infromation for the parameter formula
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
getParameters <- function(what = "mu", 
                       formula = mu.formula, 
                           fix = mu.fix, 
                         start = mu.start,
                            np = npmu,
                      start.at = "1",
                      start.v  = 1)
 {
  if (inherits(formula, "formula")) # if formula
   {
        fo2 <- inteprFormula(formula, .envir = envir, .start = start.v, .name = envname)
       npt1 <- length(attr(fo2, "parameters"))
      if (is.character(attr(fo2, "model"))) 
        {
            if (length(attr(fo2, "model")) == 1) 
            {
             fo1 <- function(p) eval(parse(text= paste(paste(paste(paste(what,
                 "h",sep="."),"(p[",sep=""), start.at, sep=""),"]*rep(1,N))",sep="")))                                                                          
             attributes(fo1) <- attributes(fo2)
             fo2 <- NULL
            }
        }
      else 
        {
            if (np  != npt1) 
            {
                cat("\nParameters are ")
                cat(attr(fo2, "parameters"), "\n")
                stop(paste(start, "should have", npt1, "estimates"))
            }
            if (is.list(start)) 
            {
                if (!is.null(names(start))) 
                {
                  o <- match(attr(fo2, "parameters"), names(start))
              start <- unlist(start)[o]
                  if (sum(!is.na(o)) != length(start)) 
                    stop("invalid estimates for", what, " - probably wrong names")
                }
                else start <- unlist(start)
            }
        }
      if (!is.null(fo2)) 
        {
               fo1 <- function(p) #mu.h(fo2(p))
                   { eval(parse(text=paste(paste(what,"h",sep="."),"(fo2(p))",sep="")))} 
   attributes(fo1) <- attributes(fo2)
        }
   } # end if formula
 else if (is.function(formula))  
   {# if function 
        fo1 <- formula
   }
   # end if
  if (!is.null(fo1) && is.null(attr(fo1, "parameters"))) 
    {
 attributes(fo1) <- if (is.function(formula)) 
                     {
                       if (!inherits(formula, "formulafn")) 
                        {
                        attributes(fnenvir(formula, .envir = envir))
                        }
                       else attributes(formula)
                     }
                    else 
                     {
            attributes(fnenvir(fo1, .envir = envir))
                     }
    }
    
    nlp <- if (is.function(fo1)) length(attr(fo1, "parameters"))
           else if (is.null(fo1))    NULL
           else npt1 #
    # end if else
   if (!is.null(nlp) && nlp != np) 
        stop(paste(what , "should have", nlp, "initial estimates"))

   if (inherits(formula, "formula") || is.function(formula)) 
    {
        if (!fix) 
        {
            if (is.numeric(start) && length(start) !=np) 
                stop(paste(what, "start must be of size ",np))
            if (!is.numeric(start)) 
                start <- rep(0,np)
                 fnfo <- fo1
        }
        else 
        {
            if (!is.numeric(start)) 
                stop("Missing initial conditions for mu")
            else if (length(start) !=np) 
                stop(paste(what, "start must be of size ",np))
            else fnfo <- function(p) fo1(start)
        }
    }
   else if (!fix) 
    {
        fnfo <- function(p) {function(p) eval(parse(text= paste(paste(paste(paste(what,"h",sep="."),"(p[",sep=""), start.at, sep=""),"]*rep(1,N))",sep="")))  }
        #mu.h(rep(p[1], N))
        if (!is.numeric(start)) 
        {
         stop(paste(what, "start must be of numeric"))
        }
    }
   else 
    {
        if (length(start) == N) 
            fnfo <- function(p) start
        else fnfo <- function(p) rep(start[1], N)
       np <- 1
    }
 list(func = fnfo, np = np , start = start)
 }
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
#  getParameters function ends here
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
## there are here to be consistent with gamlss but they have to change 
## gamlss.rc.list<-c("EX.rc","Exponential.rc") # the right censoring distribution list 
## gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
#----------------------------------------------------------------------------------------
      nlcall <- sys.call()
##-----------------------------------------------------
## this again meyby has to change if gamlssNonLinear() is callled from gamlss()
      envir <- parent.frame()
    respenv <- FALSE  
    envname <-  deparse(substitute(envir))
##-----------------------------------------------------
      if(!missing(data)) 
        if(any(is.na(data)))stop("The data contains NA's, use data = na.omit(mydata)") 
# get the data (I am not sure if this is the best way)
    if (is.data.frame(data)) { attach(data); on.exit(detach(data))}
##-----------------------------------------------------
##     get the family
##-----------------------------------------------------
      family <- as.gamlss.family(family)         
      if (!family$y.valid(y))  stop( "response variable out of range")
       fname <- family$family[1]
        dfun <- paste("d",fname,sep="")
        pfun <- paste("p",fname,sep="")
        lpar <- length(family$parameters)   
        mu.h <- family$mu.linkinv
     sigma.h <- family$sigma.linkinv 
        nu.h <- family$nu.linkinv
       tau.h <- family$tau.linkinv
#------------------------------------------------------
    if (all.fix) 
     mu.fix <- sigma.fix <- nu.fix <- tau.fix <- FALSE
       npmu <- length(mu.start)
    npsigma <- length(sigma.start)
       npnu <- length(nu.start)
      nptau <- length(tau.start)
        mu1 <- sigma1 <- nu1 <- tau1 <- NULL
       type <- "unknown"
#------ y variable
##-----------------------------------------------------------------------------------------
## This part deals with the response variable 
     if (any(is.na(y)))     stop("NAs in y - use na.omit()")
         Y <- y       
         if(is.null(dim(Y)))                       # if y not matrix
         N <- length(Y) else N <- dim(Y)[1]   # calculate the dimension for y  
## extracting now the y and the binomial denominator in case we use BI or BB
    if(any(family$family%in%gamlss:::.gamlss.bi.list)) 
    { 
       if (NCOL(Y) == 1) 
            {
            y <- if (is.factor(Y))  Y != levels(Y)[1] else Y
            bd <- rep(1, N)
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            } 
       else if (NCOL(Y) == 2) 
            {
            if (any(abs(Y - round(Y)) > 0.001)) {
            warning("non-integer counts in a binomial GAMLSS!")
                                                }
            bd <- Y[,1] + Y[,2]
            y <-  Y[,1]
            } 
       else stop(paste("For the binomial family, Y must be", 
            "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
            "and col 2 is no. failures"))
     }
     # multinomial checking
     else if(any(family$family%in%.gamlss.multin.list))
            {
               y <- if(is.factor(Y))   unclass(Y)
                    else Y
            } 
     else if(is.Surv(Y))
          { 
           ## checking that the family is censored
           if (length(grep("censored",family$family[[2]]))==0) 
            stop(paste("the family in not a censored distribution, use cens()"))
           ## checking compatability of Surv object and censored distribution
           if (length(grep(attr(Y,"type"),family$family[[2]]))==0) 
            stop(paste("the Surv object and the censored distribution are not of the same type"))
           y <- Y    
           }     
     else {y <- Y }
     if (!family$y.valid(y))  stop( "response variable out of range")
##---------------------------------------------------------------------------------------
##checking the permissible y values      
   if (!family$y.valid(y)) # MS Thursday, June 20, 2002 at 16:30 
       stop( "response variable out of range")
        censor <- FALSE
    if (length(weights) == 1)    weights <- rep(weights, N)
##---------------------mu----------------------------------------------------------------
if ("mu"%in%names(family$parameters))
{ 
  all <- getParameters(what = "mu", formula = mu.formula, fix = mu.fix, 
                       start = mu.start,  np = npmu, start.at = "1", start.v=1)
  fnmu <- all$func
  npmu <- all$np
  mu.start <- all$start
    npl2 <- if (!sigma.fix)   npmu *(!mu.fix) + 1
            else 1
}
#-----------------------sigma------------------------------------------------------------
if ("sigma"%in%names(family$parameters)) 
 {
     all <- getParameters(what = "sigma", formula = sigma.formula, fix = sigma.fix, 
                         start = sigma.start,  np = npsigma, 
                         start.at = "npl2", start.v = npl2)
     fnsigma <- all$func
     npsigma <- all$np
 start.sigma <- all$start
        npl3 <- if (!nu.fix)  npmu *(!mu.fix) +npsigma * (!sigma.fix) + 1
               else 1
 }
#----------------------------nu----------------------------------------------------------
 if ("nu"%in%names(family$parameters))
 { 
       all <- getParameters(what = "nu", formula = nu.formula, fix = nu.fix, 
                         start = nu.start,  np = npnu, 
                         start.at = "npl3", start.v = npl3)
     fnnu <- all$func
     npnu <- all$np
 start.nu <- all$start
    npl4 <- if (!tau.fix) 
       npmu *(!mu.fix) +npsigma * (!sigma.fix) +npnu * (!nu.fix) + 1
    else 1
  }
#---------------------------------tau----------------------------------------------------
  if ("tau"%in%names(family$parameters))
  {  
      all <- getParameters(what = "tau", formula = tau.formula, fix = tau.fix, 
                         start = tau.start,  np = nptau, 
                         start.at = "npl4", start.v = npl4)
     fntau <- all$func
     nptau <- all$np
 start.tau <- all$start
   }
 #------------------------------------------------------------------------------------
 ## definition of the log-likelihood function
 #------------------------------------------------------------------------------------
 #------------------------------------------------------------------------------------
   logLikelihood <- function(p) 
       {      
       if(lpar==1) 
         {   
                fmu <- fnmu(p)      
            if (exact) 
            {
                    if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu )
                     }   
                     else    
                     {
                tamp1 <- call(pfun, q = y + delta/2, mu = fmu )
                tamp2 <- call(pfun, q = y - delta/2, mu = fmu )
                     }
               
                 tamp <- eval(tamp1)-eval(tamp2)
             llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))   
                             {call(dfun, x = y, bd=bd, mu = fmu )}  
                     else    {call(dfun, x = y, mu = fmu )}   
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }
       if(lpar==2)
         {
              fmu <- fnmu(p)
            fsigma <- fnsigma(p)      
             if (exact) 
            {
             if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma )
                     }   
                     else    
                     {
                 tamp1 <- call(pfun, q = y + delta/2, mu = fmu )
                 tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma)
                     }      
                 tamp <- eval(tamp1)-eval(tamp2)
             llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
                ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    
                              {call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma)}  
                        else  {call(dfun, x = y, mu = fmu , sigma = fsigma)}  
                 tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }     
       if(lpar==3)
         {
           fmu <- fnmu(p)
        fsigma <- fnsigma(p)
           fnu <- fnnu(p)   
             if (exact) 
            {
            if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu )
                     }   
                     else    
                     {
               tamp1 <- call(pfun, q = y + delta/2, mu = fmu , sigma = fsigma, nu = fnu)
               tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma, nu = fnu)
                     }      
              tamp <- eval(tamp1)-eval(tamp2)
          llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    
                              {call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma, nu = fnu)}  
                        else  {call(dfun, x = y,        mu = fmu , sigma = fsigma, nu = fnu)}  
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
           }
       if(lpar==4)
         {
           fmu <- fnmu(p)
        fsigma <- fnsigma(p)
           fnu <- fnnu(p)
          ftau <- fntau(p)
             if (exact) 
            {
             if(any(family$family%in%gamlss:::.gamlss.bi.list))
                     {
                tamp1 <- call(pfun, q = y + delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu, tau = ftau )
                tamp2 <- call(pfun, q = y - delta/2, bd=bd, mu = fmu, sigma = fsigma, nu= fnu, tau = ftau )
                     }   
                     else    
                     {
               tamp1 <- call(pfun, q = y + delta/2, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)
               tamp2 <- call(pfun, q = y - delta/2, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)
                     }   
         tamp <- eval(tamp1)-eval(tamp2)
     llikcomp <- -(2*log(tamp)) * weights
            }
            else 
            {
             ctamp <-if(any(family$family%in%gamlss:::.gamlss.bi.list))    
                              {call(dfun, x = y, bd=bd, mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)}  
                        else  {call(dfun, x = y,        mu = fmu , sigma = fsigma, nu = fnu, tau = ftau)}  
                tamp <-eval(ctamp)
               # llikcomp <- -(log(tamp) + log(delta)) * weights
                llikcomp <- -2*log(tamp)* weights
            }
         }     
        llik <- sum(llikcomp)
        if (llik.output) 
        {
            if (length(p) == 0) 
                cat("-LogLik: ", sum(llikcomp), "\n")
            else cat("-LogLik: ", sum(llikcomp), " ", p, "\n")
        }
               
    z <-  if (lpar==1) list(llik = llik, llikcomp = llikcomp, mu = fmu)
          else if (lpar==2) list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma)
          else if (lpar==3) list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma, nu = fnu)#
          else  list(llik = llik, llikcomp = llikcomp, mu = fmu, sigma = fsigma, nu = fnu, tau = ftau)
    z
       }
    #------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------
    # the function to be mimimized by nlm()
    #------------------------------------------------------------------------------------
    optFunction <- function(p) 
    {
        tamp <-logLikelihood(p)$llik
        if (llik.output) 
            cat("-LogLik: ", tamp, " (", p, ")", "\n")
        if (is.na(tamp)) 
            1e+20
        else tamp
    }
    #----------------------------------------------------------------------------------
    p0 <- c()
    if (!mu.fix) {
        p0 <- c(mu.start)
        names(p0) <- c(rep("mu.start", length(mu.start)))
    }
    if (!sigma.fix) {
        tamp <- names(p0)
        p0 <- c(p0, sigma.start)
        names(p0) <- c(tamp, rep("sigma.start", length(sigma.start)))
    }
    if (!nu.fix) {
        tamp <- names(p0)
        p0 <- c(p0, nu.start)
        names(p0) <- c(tamp, rep("nu.start", length(nu.start)))
    }
    if (!tau.fix) {
        tamp <- names(p0)
        p0 <- c(p0, tau.start)
        names(p0) <- c(tamp, rep("tau.start", length(tau.start)))
    }
    if (llik.output) 
       {
        cat("No. of parameters: ",npmu, "",npsigma, "",npnu, 
            "",nptau, "\n")
        if (!mu.fix || !sigma.fix || !nu.fix || !tau.fix) 
        {
            cat("Vector of initial conditions on IR^p:", "\n")
            print(p0)
        }
       }
    #browser()
   # llik0 <-logLikelihood(p = p0)
    np0 <- length(p0)
    #----basic nlm run
    # getting the control parameters
    typsize <- if(is.null(control$typsize)) abs(p0) else control$typsize
    if(any(typsize <= 0)) 
      {
      warning("the value of typsize supplied is zero or negative the default value of abs(p0) was used instead")
      typsize <- abs(p0)+0.001 #MSThursday, June 1, 2006 at 10:13
      }
    stepmax <- if(is.null(control$stepmax)) sqrt(p0 %*% p0) else control$stepmax
    if(stepmax < 0) 
      {
    warning("the value of stepmax supplied is zero or negative the default value of  sqrt(p0 %*% p0) was used instead")
         stepmax <-  sqrt(p0 %*% p0)
      }
    hessian <- control$hessian
    fscale  <- control$fscale
print.level <- control$print.level
     ndigit <- control$ndigit
    gradtol <- control$gradtol
    steptol <- control$steptol
    iterlim <- control$iterlim
    #------------------------------------------------------------------------------------
    # the actual non linear fitting here                                                |
    #------------------------------------------------------------------------------------
    if (np0 > 0) {
        p.opt <- nlm(optFunction, p = p0, hessian = hessian, fscale = fscale, 
            typsize = rep(1, length(p0)), print.level = print.level, 
            ndigit = ndigit, gradtol = gradtol, steptol = steptol, 
            iterlim = iterlim, stepmax = stepmax)
        z <-logLikelihood(p.opt$estimate)
    }
    else z <-logLikelihood(p0)
    #------------------------------------------------------------------------------------
 coefficients <- p.opt$estimate
# browser()
   #  gradient <- p.opt$gradient
   #     error <- p.opt$error
   #np <- if (lpar==1) npmu   
   #     else if (lpar==2) npmu +npsigma  
   #     else if (lpar==3)npmu +npsigma +npnu
   #     else npmu +npsigma +npnu +nptau
         nobs <- sum(as.numeric(weights))
           mu <- as.vector(z$mu)
        sigma <- as.vector(z$sigma)
           nu <- as.vector(z$nu)
          tau <- as.vector(z$tau)
#----------------------------------------------------------------------------------------
  if (np0 > 0) 
   {
 #   browser()
        cov <- diag(np0)
        if (hessian) {
            if (np0 == 1) 
                cov <- 2/p.opt$hessian
                #cov <- 1/(2*p.opt$hessian) # bug corrected my Tom Jagger
            else {
                a <- if (any(is.na(p.opt$hessian)) || any(abs(p.opt$hessian) == Inf)) 0
                     else qr(p.opt$hessian)$rank
                if (a == np0) 
                  cov <- solve(p.opt$hessian/2)
                else cov <- matrix(NA, ncol = np0, nrow = np0)
            }
        }
          se <- if (hessian) sqrt(diag(cov))
                else NA
        corr <- if (hessian) cov2cor(cov) # cov/(se %o% se) 
                else NA
   }
  else coefficients <- se <- cov <- corr  <-  NULL
#------------------ output --------------------------------------------------------------
    out <- list(family = family$family , 
            parameters = names(family$parameters), 
                  call = nlcall, 
                     y = y, 
               control = control, 
               weights = weights, 
            G.deviance = z$llik, 
     #      P.deviance = z$llik, 
                     N = nobs, 
                 rqres = family$rqres,  
                  type = family$type, 
                   aic = z$llik + 2*np0, 
                   sbc = z$llik + log(nobs)*np0,
           df.residual = nobs - np0,
                df.fit = np0,
             converged = p.opt$code, 
                  iter = p.opt$iter, 
      #            pen = 0,
             residuals = eval(family$rqres), 
                method = "JL()",
          coefficients = coefficients,
                    se = se,
                   cov = cov,
                  corr = corr
                    )
   if(any(family$family%in%gamlss:::.gamlss.bi.list))  out$bd <- bd                           
#========================================================================================
##---------------------------------------------------------------------------------------
## this function is used in for outputing the parameters
##---------------------------------------------------------------------------------------
##=======================================================================================
parameterOut <- function(what="mu")
 {
  out <- list() 
  if(family$parameter[[what]]==TRUE && eval(parse(text=paste(what,".fix",sep="")))==FALSE)
     {              
  out$fv <- eval(parse(text=paste("z", what ,sep="$")))
      out$lp <-  eval(parse(text=paste(paste(paste("family", what ,sep="$"),"linkfun(",sep="."),"out$fv)", sep="")))
    out$link <- eval(parse(text=paste(paste("family", what ,sep="$"),"link",sep=".")))
 out$formula <- eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"formula\")",sep=",")))
       names <- if(is.expression(eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"model\")",sep=",")))))
                 {
               eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"parameters\")",sep=",")))
                 }
                else  
                 { 
               eval(parse(text=paste(paste("attr(",paste("fn",what,sep=""),sep=""),"\"model\")",sep=",")))
                 }
out$coefficients <- coefficients[paste(what,"start", sep=".")==names(p0)]
names(out$coefficients)<-names
      out$se <- se[paste(what,"start", sep=".")==names(p0)]
names(out$se)<- names   
      out$df <- eval(parse(text=(paste("np", what, sep="")))) 
   #out$nl.df <- 0
   #out$pen  <- 0 
     }
   else
     { out$fix <- eval(parse(text=paste(what,".fix",sep="")))
        out$df <- 0
        out$fv <- eval(parse(text=paste("z", what ,sep="$")))
     }
#if (is.data.frame(data)) detach(data)
  out
  }
##----------------------------------------------------------------------------------------
##========================================================================================
##  Output for mu model: -----------------------------------------------------------------
 if ("mu"%in%names(family$parameters))   out <- c(out, mu = parameterOut(what="mu") )
 else            out$mu.df <- 0
##  Output for sigma model: --------------------------------------------------------------
 if ("sigma"%in%names(family$parameters))out <- c(out, sigma = parameterOut(what="sigma"))
 else            out$sigma.df <- 0
## Output for nu model: ------------------------------------------------------------------
 if ("nu"%in%names(family$parameters))   out <- c(out, nu = parameterOut(what="nu") )
 else              out$nu.df <- 0
## output for tau model ------------------------------------------------------------------
 if ("tau"%in%names(family$parameters))  out <- c(out, tau = parameterOut(what="tau") )
else               out$tau.df <- 0
    class(out) <- list("gamlssNonLinear", "gamlss")
#if (is.data.frame(data)) detach(data)
    out
} 
#-----------------------------------------------------------------------------------------  
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
"funMobJ" <- function (z, envir = parent.frame()) 
{
    if (!inherits(z, "formula")) 
        return(NULL)
    local <- fac <- cov <- ch <- NULL
    ch1 <- deparse(z[[length(z)]])
    for (j in 1:length(ch1)) ch <- paste(ch, ch1[j], collapse = " ", 
        sep = "\n")
    ch <- gsub("\n", " \n ", gsub("\\^", " ^ ", gsub("\\)", " )", 
        gsub("\\[", " [ ", ch))))
    ch <- gsub("/", " / ", gsub(",", " ,", gsub("\\(", " ( ", 
        gsub(":", " : ", ch))))
    ch <- paste(" ", gsub(" -", " - ", gsub(" \\+", " + ", ch)), 
        " ", sep = "")
    mem <- all.vars(z)
    fcn <- all.names(z)
    fcn <- fcn[!(fcn %in% mem)]
    if (length(mem) > 0) {
        tmp <- vector(length = length(mem))
        for (i in 1:length(mem)) tmp[i] <- exists(mem[i], envir = envir) && 
            is.function(eval(parse(text = mem[i]), envir = envir)) && 
            (length(grep(paste(mem[i], ","), ch)) > 0 || length(grep(paste(",", 
                mem[i]), ch)) > 0 || length(grep(paste("=", mem[i]), 
                ch)) > 0 || length(grep(paste("\\(", mem[i], 
                "\\)"), ch)) > 0) && mem[i] != "times" && mem[i] != 
            "response"
        fcn <- unique(c(fcn, mem[tmp]))
        mem <- mem[!tmp]
    }
    if (length(mem) > 0) {
        for (j in 1:length(mem)) {
            local <- c(local, length(grep(paste(mem[j], "<-"), 
                ch)) > 0)
            cov <- c(cov, exists(mem[j], envir = envir) && is.numeric(eval(parse(text = mem[j]), 
                envir = envir)))
            fac <- c(fac, exists(mem[j], envir = envir) && !cov[j] && 
                is.factor(eval(parse(text = mem[j]), envir = envir)))
        }
        zz <- list(formula = ch, objects = mem, functions = fcn, 
            parameters = !cov & !fac & !local, covariates = cov, 
            factors = fac, local = local)
    }
    else zz <- list(formula = ch, objects = mem, functions = fcn)
    class(zz) <- "funMobJ"
    zz
}
#----------------------------------------------------------------------------------------
functionEnvironment <- function (.z, ...) 
UseMethod("functionEnvironment")
#----------------------------------------------------------------------------------------
functionEnvironment.default <- function (.z, .envir = parent.frame(), .name = NULL, .expand = TRUE, 
                            .response = FALSE) 
{
    if (!is.function(.z)) 
        return(NULL)
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    if (!is.environment(.envir)) {
        if (is.null(.name)) 
            .name <- paste(deparse(substitute(.envir)))
        #if (inherits(.envir, "repeated")) 
        #    return(functionEnvironment.repeated(.z, .envir, .name = .name, 
        #        .expand, .response))
        #if (inherits(.envir, "tccov")) 
        #    return(functionEnvironment.tccov(.z, .envir, .name = .name, .expand))
        #if (inherits(.envir, "tvcov")) 
        #    return(functionEnvironment.tvcov(.z, .envir, .name = .name, .expand))
        if (inherits(.envir, "data.frame")) 
            return(functionEnvironment.data.frame(.z, .envir, .name = .name, 
                .expand))
    }
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    if (length(.mem) > 0) {
        for (.j in 1:length(.mem)) {
            .ex <- c(.ex, exists(.mem[.j], envir = .envir))
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j])))
            } else FALSE)
        }
        for (.j in 1:length(.mem)) {
            if (!.fcn[.j] && .ex[.j] && is.factor(eval(parse(text = .mem[.j]), 
                envir = .envir))) 
                stop(paste(.mem[.j], "is a factor variable"))
        }
        .un <- unique(.mem[!.ex])
        if (length(unique(.mem[.ex & !.fcn])) == 0 && length(.un) == 
            0) 
            warning("functionEnvironment.default: no variables found")
    }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex <- if (length(.fcn) > 0 && !is.null(.ex)) 
        .ex & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
functionEnvironment.data.frame <- function (.z, .envir = NULL, .name = NULL, .expand = TRUE) 
{
    if (!is.function(.z)) 
        return(NULL)
    .ndata <- if (is.null(.name)) 
        paste(deparse(substitute(.envir)))
    else .name
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex1 <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    .cn <- colnames(.envir)
    if (length(.mem) > 0) {
        .ex1 <- match(.mem, .cn)
        for (.j in 1:length(.mem)) {
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j]))) && 
                  is.na(.ex1[.j])
            } else FALSE)
        }
        .un <- unique(.mem[is.na(.ex1) & !.fcn])
        if (length(unique(.mem[!is.na(.ex1) & !.fcn])) == 0 && 
            length(.un) == 0) 
            warning("functionEnvironment.data.frame: no variables found")
    }
    for (.j in 1:length(.ch1)) {
        .ch1[.j] <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", 
            " [", .ch1[.j])))
        .ch1[.j] <- gsub("-", "- ", gsub("/", " / ", gsub(",", 
            " ,", .ch1[.j])))
        .ch1[.j] <- paste(" ", gsub("\\(", " ( ", .ch1[.j]), 
            " ", sep = "")
    }
    .ex1a <- .ex1[!is.na(.ex1)]
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) {
            if (is.factor(.envir[, .ex1a[.j]])) 
                stop(paste(colnames(.envir)[.ex1a[.j]], "is a factor variable"))
            for (.k in 1:length(.ch1)) .ch1[.k] <- gsub(paste(" ", 
                .cn[.ex1a[.j]], " ", sep = ""), paste(" ", .ndata, 
                "$", .cn[.ex1a[.j]], sep = ""), .ch1[.k])
        }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex1 <- if (!is.null(.ex1) && length(.fcn) > 0) 
        !is.na(.ex1) & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex1]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
vcov.gamlssNonLinear <- function (object, ...) 
{ 
    object$cov
}
#========================================================================================
#========================================================================================
# 
fnenvir <- function (.z, ...) 
UseMethod("fnenvir")
#----------------------------------------------------------------------------------------
fnenvir.default <- function (.z, .envir = parent.frame(), .name = NULL, .expand = TRUE, 
                            .response = FALSE) 
{
    if (!is.function(.z)) 
        return(NULL)
    if (is.name(.envir)) {
        if (is.null(.name)) 
            .name <- as.character(.envir)
        .envir <- eval(.envir)
    }
    if (!is.environment(.envir)) {
        if (is.null(.name)) 
            .name <- paste(deparse(substitute(.envir)))
       # if (inherits(.envir, "repeated")) 
       #     return(fnenvir.repeated(.z, .envir, .name = .name, 
       #         .expand, .response))
       # if (inherits(.envir, "tccov")) 
       #     return(fnenvir.tccov(.z, .envir, .name = .name, .expand))
       # if (inherits(.envir, "tvcov")) 
       #     return(fnenvir.tvcov(.z, .envir, .name = .name, .expand))
        if (inherits(.envir, "data.frame")) 
            return(fnenvir.data.frame(.z, .envir, .name = .name, 
                .expand))
    }
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    if (length(.mem) > 0) {
        for (.j in 1:length(.mem)) {
            .ex <- c(.ex, exists(.mem[.j], envir = .envir))
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j])))
            } else FALSE)
        }
        for (.j in 1:length(.mem)) {
            if (!.fcn[.j] && .ex[.j] && is.factor(eval(parse(text = .mem[.j]), 
                envir = .envir))) 
                stop(paste(.mem[.j], "is a factor variable"))
        }
        .un <- unique(.mem[!.ex])
        if (length(unique(.mem[.ex & !.fcn])) == 0 && length(.un) == 
            0) 
            warning("fnenvir.default: no variables found")
    }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex <- if (length(.fcn) > 0 && !is.null(.ex)) 
        .ex & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
fnenvir.data.frame <- function (.z, .envir = NULL, .name = NULL, .expand = TRUE) 
{
    if (!is.function(.z)) 
        return(NULL)
    .ndata <- if (is.null(.name)) 
        paste(deparse(substitute(.envir)))
    else .name
    .ch1 <- deparse(.z, width.cutoff = 500)
    .ch2 <- .ch1[1]
    .ch1 <- .ch1[-1]
    .mem2 <- strsplit(gsub("[(),]", " ", .ch2), " ")[[1]]
    if (length(.mem2) > 0) 
        .mem2 <- .mem2[.mem2 != ""]
    if (length(.mem2) > 1) 
        .mem2 <- .mem2[2:length(.mem2)]
    else .mem2 <- NULL
    .fcn <- .ex1 <- .ch <- NULL
    for (.j in 1:length(.ch1)) .ch <- paste(.ch, .ch1[.j], collapse = " ")
    .mem <- strsplit(gsub("(\\[(0|1|2|3|4|5|6|7|8|9|:|,)+\\])|([][+*/^():!<>%&|~,{}\"\\=-])|( [0-9]+)|(\\.[0-9]+)|(^[0-9]+)", 
        " ", .ch), " ")[[1]]
    if (length(.mem) > 0) 
        .mem <- .mem[.mem != ""]
    .cn <- colnames(.envir)
    if (length(.mem) > 0) {
        .ex1 <- match(.mem, .cn)
        for (.j in 1:length(.mem)) {
            .fcn <- c(.fcn, if (exists(.mem[.j])) {
                if (.mem[.j] == "function" || .mem[.j] == "if" || 
                  .mem[.j] == "else" || .mem[.j] == "for" || 
                  .mem[.j] == "while" || .mem[.j] == "repeat") TRUE else is.function(eval(parse(text = .mem[.j]))) && 
                  is.na(.ex1[.j])
            } else FALSE)
        }
        .un <- unique(.mem[is.na(.ex1) & !.fcn])
        if (length(unique(.mem[!is.na(.ex1) & !.fcn])) == 0 && 
            length(.un) == 0) 
            warning("fnenvir.data.frame: no variables found")
    }
    for (.j in 1:length(.ch1)) {
        .ch1[.j] <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", 
            " [", .ch1[.j])))
        .ch1[.j] <- gsub("-", "- ", gsub("/", " / ", gsub(",", 
            " ,", .ch1[.j])))
        .ch1[.j] <- paste(" ", gsub("\\(", " ( ", .ch1[.j]), 
            " ", sep = "")
    }
    .ex1a <- .ex1[!is.na(.ex1)]
    if (length(.ex1a) > 0) 
        for (.j in 1:length(.ex1a)) {
            if (is.factor(.envir[, .ex1a[.j]])) 
                stop(paste(colnames(.envir)[.ex1a[.j]], "is a factor variable"))
            for (.k in 1:length(.ch1)) .ch1[.k] <- gsub(paste(" ", 
                .cn[.ex1a[.j]], " ", sep = ""), paste(" ", .ndata, 
                "$", .cn[.ex1a[.j]], sep = ""), .ch1[.k])
        }
    .ch <- gsub("\\^", " ^ ", gsub("\\)", " )", gsub("\\[", " [", 
        .ch)))
    .ch <- gsub("-", "- ", gsub("/", " / ", gsub(",", " ,", gsub("\\+", 
        "+ ", .ch))))
    .ch <- paste(" ", gsub("\\(", " ( ", gsub(":", " : ", .ch)), 
        " ", sep = "")
    .ch2 <- strsplit(.ch, " ")[[1]]
    .un <- .un0 <- .un1 <- NULL
    if (length(.mem2) > 0) 
        for (.j in 1:length(.mem2)) {
            .ex1a <- NULL
            for (.k in 1:length(.ch2)) if (.mem2[.j] == .ch2[.k]) {
                if (.k < length(.ch2) && length(grep("^\\[", 
                  .ch2[.k + 1])) > 0) {
                  .ex1a <- c(.ex1a, paste(.ch2[.k], .ch2[.k + 
                    1], sep = ""))
                  .un1 <- c(.un1, .ch2[.k])
                }
                else .un0 <- c(.un0, .ch2[.k])
            }
            if (!is.null(.ex1a)) {
                .ex1a <- unique(.ex1a)
                .o <- gsub("(^[[:alnum:]]\\[)|(\\])", "", .ex1a)
                .un <- if (length(grep("[[:alpha:]]", .o)) > 
                  0) 
                  c(.un, .ex1a)
                else c(.un, .ex1a[order(as.numeric(.o))])
            }
        }
    if (length(.un0) > 0) {
        if (length(.un1) > 0) {
            .tmp <- NULL
            for (.k in 1:length(.un1)) if (length(grep(.un1[.k], 
                .un0)) > 0) 
                .tmp <- c(.tmp, grep(.un1[.k], .un0))
            .un <- c(.un, unique(if (!is.null(.tmp)) .un0[-.tmp] else .un0))
        }
        else .un <- c(.un, unique(.un0))
    }
    .fnb <- eval(parse(text = paste("function(", paste(.mem2, 
        collapse = ","), ")", paste("eval(attr(.fnb,\"model\"))"))))
    .ex1 <- if (!is.null(.ex1) && length(.fcn) > 0) 
        !is.na(.ex1) & !.fcn
    else NULL
    attributes(.fnb) <- list(model = parse(text = .ch1), parameters = .un, 
        covariates = unique(.mem[.ex1]), class = "formulafn")
    .obj <- ls(all.names = TRUE)
    rm(list = .obj[.obj != ".fnb"])
    rm(.obj)
    return(.fnb)
}
#----------------------------------------------------------------------------------------
