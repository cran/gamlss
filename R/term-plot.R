### function to plot additive terms in the fitting of a GAMLSS model
### it is based on the termplot() function in R
### but there were are several bugs corrected
### 1. if what parameters is constant an error is given
### 2. if interaction are present then a warning is given 
###    rather that crashing as is termplot()
### TO DO
###    i) to allow having plots in one screen the "pages" option  OK
###   ii) to allow plotting as shades the "scheme" option         OK
###  iii) to allow plotting class object   Partial OK it needs to go trought all smoothers
### Author: Mikis Stasinopoulos
### bugs: if the additve terms depends on more that one variable 
###       i.e. loess(x1,x2) produces rubish
term.plot <- function (object, 
                        what = c("mu","sigma","nu","tau"),  
                        data = NULL, 
                       envir = environment(formula(object)), 
               partial.resid = FALSE, 
                         rug = FALSE, 
                       terms = NULL, 
                          se = TRUE, 
                        ylim = c("common","free"),
                      scheme = c( "shaded", "lines"),
                       xlabs = NULL, 
                       ylabs = NULL, 
                        main = NULL, 
                       pages = 0, #  New 
                    col.term = "darkred",
                      col.se = "orange", 
                   col.shaded = "gray", 
                     col.res = "lightblue", 
                     col.rug = "gray",
                    lwd.term = 1.5,   
                      lty.se = 2, 
                      lwd.se = 1,
                     cex.res = 1, 
                     pch.res = par("pch"),
                         ask = interactive() && nb.fig < n.tms &&.Device != "postscript", #dev.interactive() && nb.fig < n.tms
           use.factor.levels = TRUE, 
                 surface.gam = FALSE, 
#                 gam.scheme = 3,
                             ...) 
{
#-------------------------------------------------------------------------------
# Local functions
# i) CheckSmoList() to chech whether they are smoothers in terms
#-------------------------------------------------------------------------------
# the results is something like 
#     c(0,0,1,2,2,2,3,4) if smootherrs where used in positions
#     c(0,0,1,1,0,0,1,1) 
CheckSmoList <- function(termList)
  {
    gamlss.sm.list1 <- c( "cs","scs", "ps", "pb", "cy", "pvc", "pbm",  "pbj",   
                         "mrf",   "mrfa", "sap",  "krig",   "lo", "random",
                         "re",  "fp", "pp", "nl","ri","ridge","fk", "la",     
                         "tr",  "ga",   "nn" )
    # ideally this should be done autonmatically
    gamlss.sm.list2 <- c( "cs(","scs(", "ps(", "pb(", "cy(", "pvc(", "pbm(",  "pbj(",   
                         "mrf(",   "mrfa(", "sap(",  "krig(",   "lo(", "random(",
                         "re(",  "fp(", "pp(", "nl(","ri(","ridge(","fk(", "la(",     
                         "tr(",  "ga(",   "nn(" )
    lgamsmol  <- length(gamlss.sm.list1)
         lsm  <- length(termList)
          res <- rep(0, lsm) 
          att <- rep(0, lsm)
    for (i in 1:length(gamlss.sm.list2))
    {
           LL <-   grepl(gamlss.sm.list2[i], termList, fixed=TRUE)
         res  <- res+LL
      att[LL] <- gamlss.sm.list1[i]
    }
    res <- cumsum(res)
    attr(res, "whichSmo") <- att
    res
  }
#----end of local function------CheckSmoList------------------------------------  
CheckSmoWithPlot <- function(termList)
{
  #gamlss.Smo.plot.list <- c( "tr", "ga")
  gamlss.Smo.plot.list1 <- c( "tr(", "ga(", "nn(")
  lgamsmol  <- length(gamlss.Smo.plot.list1)
  lsm  <- length(termList)
  res <- rep(0, lsm) 
  for (i in gamlss.Smo.plot.list1)
  {
    LL <-   grepl(i, termList,  fixed=TRUE) 
    res  <- res+LL
  }
res  
}

# more local functions----------------------------------------------------------
carrier <- function(term) 
  {
  if (length(term) > 1)   carrier(term[[2]])
  else eval(term, data, enclos = pf)
  }# with envir  Global
#-------------------------------------------------------------------------------
carrier.name <- function(term) 
    {
    if (length(term) > 1) carrier.name(term[[2]])
    else as.character(term)
    }
#------MORE LOCAL FUNCTIONS-----------------------------------------------------
se.lines <- function(x, iy, i, ff = 2)             
{
  tt <- ff * terms$se.fit[iy, i]
  lines(x, tms[iy, i] + tt, lty = lty.se, lwd = lwd.se, col = col.se)
  lines(x, tms[iy, i] - tt, lty = lty.se, lwd = lwd.se, col = col.se)
}
#-------------------------------------------------------------------------------
se.shaded <- function(x, iy, i, ff = 2) 
{
  tt <- ff * terms$se.fit[iy, i]
  xx <- c(x,rev(x))
  yy <- c(tms[iy, i] - tt, rev(tms[iy, i] + tt))
  polygon(xx, yy, col = col.shaded, border = col.shaded)
}
#-------------------------------------------------------------------------------
#--------END of local functions-------------------------------------------------
#-----end of local functions----------------------------------------------------  
#-------------------------------------------------------------------------------
#  begining of the proper function
#-------------------------------------------------------------------------------
## only for gamlss objects 
    if (!is.gamlss(object))  stop(paste("This is not an gamlss object", "\n", "")) 
           what <- match.arg(what)
           ylim <- match.arg(ylim)
         scheme <- match.arg(scheme)
## checking the parameter
    if (!what%in%object$par) stop(paste(what,"is not a parameter in the object","\n"))
## getting the specific term for plotting
    which.terms <- terms
## get all terms and attributes 
      par.terms <- object[[paste(what, "terms", sep=".")]]
      par.attr  <- attributes(par.terms)
          terms <- if (is.null(terms)) lpred(object, what = what, type = "terms", se.fit = se)
                  else lpred(object, what = what, type = "terms", se.fit = se, terms = terms)
# the number of terms
          n.tms <- ncol(tms <- as.matrix(if (se) terms$fit  else terms))
## if the parameters has only a constant fitted stop
     if (n.tms == 0)
          stop("The model for ", what, " has only the constant fitted") 
## model frame    
# ff<-paste(paste("object$", what, sep=""),".formula", sep="")
# get_all_vars(object$mu.formula, data=object$call$data)
             mf <- model.frame(object, what = what) 
#            mf <- model.frame.gamlssT(object, what = what) 
## this take care if data is used         
    if (is.null(data))  # get the data from gamlss data 
           data <- eval(object$call$data, envir)
    if (is.null(data))  # if still null get from model frame
           data <- mf
    if (NROW(tms) < NROW(data)) 
     {
       use.rows <- match(rownames(tms), rownames(data))
     }
    else use.rows <- NULL
            nmt <- colnames(tms)                    # the names of terms 
## whether there are interactions in the model (they are difficult to plot)
   Interactions <- par.attr$order > 1 
## if interaction nmt has to change -------------------------------------------
      if (any(Interactions))
       {
                    nmt <- nmt[!Interactions] # take out interactions
         if (!se) 
          { 
                  terms <- terms[,nmt, drop = FALSE]
          }
         else 
          { 
              terms$fit <- terms$fit[,nmt,  drop = FALSE]  
           terms$se.fit <- terms$se.fit[,nmt,  drop = FALSE] 
          }
                  n.tms <- ncol(tms <- as.matrix(if (se) terms$fit  else terms))
          # I am assuming that 'terms' will be used wisely here 
          warning("interactions have been taken out from the plots")
       }
#-------------------------------------------------------------------------------
             cn <- parse(text = nmt) # as expression
ifSpecialSmo  <- CheckSmoWithPlot(nmt) #????????????????????????????????????????
whichValueSmo <- CheckSmoList(nmt)
# if (!is.null(smooth)) # I do not need this but match.fun() is very interesting 
#         smooth <- match.fun(smooth) 
    if (is.null(ylabs)) 
          ylabs <- paste("Partial for", nmt) # get the labels 
    if (is.null(main)) 
           main <- ""
    else if (is.logical(main)) 
           main <- if (main)  deparse(object$call, 500)
                  else ""
                 else if (!is.character(main)) 
                   stop("`main' must be TRUE, FALSE, NULL or character (vector).")
          main <- rep(main, length = n.tms)
            pf <- envir  # the carrier has different envir (the Global)

    if (is.null(xlabs))#  get the x lables
       {
       xlabs <- unlist(lapply(cn, carrier.name))
       }
    if (partial.resid ) # || !is.null(smooth) 
       {
        pres <- residuals(object, what = what, type="partial")
        #pres <- pres # mayby not - mean(pres)	# Dilip suggest  centering
        if (!is.null(which.terms)) 
        pres <- pres[, which.terms, drop = FALSE]
        if (any(Interactions))
        pres <- pres[, nmt, drop = FALSE]
       }
        is.fac <- sapply(nmt, function(i) is.factor(mf[, i])) # whether factors


        nb.fig <- prod(par("mfcol")) # the number of figures
    if (ask) {
            op <- par(ask = TRUE)
             on.exit(par(op))
              }
         ylims <- ylim # default "common"
    if (identical(ylims, "common"))  # whether common limit in y
    {
         ylims <- if (!se) 
        range(tms, na.rm = TRUE)
      else range(tms + 1.05 * 2 * terms$se.fit, tms - 1.05 * 
                   2 * terms$se.fit, na.rm = TRUE)
      if (partial.resid) 
        ylims <- range(ylims, pres, na.rm = TRUE)
      if (rug) 
        ylims[1L] <- ylims[1L] - 0.07 * diff(ylims)
    } # finnish  limits for y
#-------------------------------------------------------------------------------
# pages 
     n.plots <- n.tms
    if (pages > n.plots)  pages <- n.plots # if pages bigger use no of plots 
    if (pages < 0)        pages <- 0       # set to zero
    if (pages != 0)                        # if pages are not zero 
       {
         ppp <- n.plots%/%pages            # no of plots per page
         if (n.plots%%pages != 0) 
          {
            while (ppp * (pages - 1) >= n.plots) pages <- pages - 1
          }
                      c <- r <- trunc(sqrt(ppp))
          if (c < 1)  r <- c <- 1
          if (c * r < ppp) c <- c + 1
          if (c * r < ppp) r <- r + 1
                 oldpar <- par(mfrow = c(r, c))
        }
   else {  # if pages are zero
        ppp <- 1
    oldpar <- par()
        }
   if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) || 
      pages > 1 && dev.interactive()) 
         ask <- TRUE
    else ask <- FALSE
   if (!is.null(terms)) { ask <- FALSE}
   if (ask) 
      {
     oask <- devAskNewPage(TRUE)
     on.exit(devAskNewPage(oask))
       }
#-------------------------------------------------------------------------------
    for (i in 1:n.tms) # START of the 1:n.tms loop
    {
      # if we need differnt y limit for each variable 
        if (identical(ylim, "free"))
          { 
               ylims <- range(tms[, i], na.rm = TRUE)
            if (se) 
               ylims <- range(ylims, tms[, i] + 1.05 * 2 * terms$se.fit[, 
                            i], tms[, i] - 1.05 * 2 * terms$se.fit[, i], 
                            na.rm = TRUE)
            if (partial.resid)  
              ylims <- range(ylims, pres[, i], na.rm = TRUE)
            if (rug) 
           ylims[1] <- ylims[1] - 0.07 * diff(ylims)
           }
        if (!ifSpecialSmo[i])#--------------------------------------------------
        { #   if is a normal smoother use this
        # if is a factor--------------------------------------------------------
          if (is.fac[i]) # if the term is a factor
          {
            ff <- mf[, nmt[i]]
            if (!is.null(object$na.action)) 
              ff <- naresid(object$na.action, ff)
              ll <- levels(ff)
           xlims <- range(seq(along = ll)) + c(-0.5, 0.5)
              xx <- as.numeric(ff)
            if (rug) 
            {
              xlims[1] <- xlims[1] - 0.07 * diff(xlims)
              xlims[2] <- xlims[2] + 0.03 * diff(xlims)
            }
            plot(1, 0, type = "n", xlab = xlabs[i], ylab = ylabs[i], 
                 xlim = xlims, ylim = ylims, main = main[i], xaxt = "n", 
                 ...)
            if (use.factor.levels) 
              axis(1, at = seq(along = ll), labels = ll, ...)
            else axis(1)
            for (j in seq(along = ll)) 
            {
              ww <- which(ff == ll[j])[c(1, 1)]
              jf <- j + c(-0.4, 0.4)  
              if (se) 
                if (identical(scheme, "lines"))
                {
                  se.lines(jf, iy = ww, i = i)
                  lines(jf, tms[ww, i], col = col.term, lwd = lwd.term, 
                        ...)  
                }
                if (identical(scheme, "shaded"))
                {
                se.shaded(jf, iy = ww, i = i)
                lines(jf, tms[ww, i], col = col.term, lwd = lwd.term, 
                      ...)  
               }
              else lines(jf, tms[ww, i], col = col.term, lwd = lwd.term, 
                         ...) 
            }
          } # end if factor  ---------------------------------------------------
          else 
          { # here is where changes had to be made at the moment every pass
            # cn is expression
              xx <- carrier(cn[[i]]) # ds Friday, October 9, 2009 at 13:13
            # why we need this ??? is it for random??
            if (is.factor(xx)) xx <- seq(along = levels(xx))
            if (!is.null(use.rows)) # in case some of the rows are not used 
               xx <- xx[use.rows]
            xlims <- range(xx, na.rm = TRUE)
            if (rug) 
              xlims[1] <- xlims[1] - 0.07 * diff(xlims)
               oo <- order(xx)
            if (identical(scheme, "lines"))
            {
              plot(xx[oo], tms[oo, i], type = "l", xlab = xlabs[i], 
                   ylab = ylabs[i], xlim = xlims, ylim = ylims, 
                   main = main[i], col = col.term, lwd = lwd.term, 
                   ...)
              if (se) 
                se.lines(xx[oo], iy = oo, i = i)
            }
            if (identical(scheme, "shaded"))
            {
              if (se) 
                plot(xx[oo], tms[oo, i], type = "n", xlab = xlabs[i], 
                     ylab = ylabs[i], xlim = xlims, ylim = ylims, 
                     main = main[i], col = "black", lwd = lwd.term, 
                     ...)          
              se.shaded(xx[oo], iy = oo, i = i)
              lines(xx[oo], tms[oo, i], col = col.term,lwd = lwd.term,...)          
            }  
           } # end continuous and normal smoothers
          if (partial.resid) 
          {
            points(xx, pres[, i], cex = cex.res, pch = pch.res, 
                   col = col.res)
          }
          if (rug) 
          {
                n <- length(xx)
            lines(rep.int(jitter(xx), rep.int(3, n)), rep.int(ylims[1] + 
                                                                c(0, 0.05, NA) * diff(ylims), n), col=col.rug)
            if (partial.resid) 
              lines(rep.int(xlims[1] + c(0, 0.05, NA) * diff(xlims), 
                            n), rep.int(pres[, i], rep.int(3, n)),  col=col.rug)
          }  
        } # end of all normal smoother 
        else
        { # the special smoothers who have a plotting function 
          if (attr(whichValueSmo, "whichSmo")[i]=="ga"&&surface.gam==TRUE)
           {
            plot(getSmo(object, what, which=whichValueSmo[i]), scheme=1)
           } 
          if (attr(whichValueSmo, "whichSmo")[i]=="nn")
          {
            plot(getSmo(object, what, which=whichValueSmo[i]), y.lab=expression(eta))
          } 
          else
          {
            plot(getSmo(object, what, which=whichValueSmo[i]))
          }  
          if (attr(whichValueSmo, "whichSmo")[i]=="tr") 
             text(getSmo(object, what, which=whichValueSmo[i]))
        }
    } # end of the  1:n.tms lop -----------------------------------------------
    if (pages > 0) par(oldpar)
    invisible(n.tms)
}

# ################################################################################
# # my attempt to modify model.frame.gamlss but 
# #  get_all_vars() do not get what I want
# # new MS Thursday, June 24, 2004 at 13:45
# model.frame.gamlssT <-function (formula, what = c("mu", "sigma", "nu", "tau"), ...) 
# {
#   object <- formula
#   dots <- list(...)
#   what <- match.arg(what) 
#   Call <- object$call
#   parform <- formula(object, what)
#   #parform <- object[[paste(what, "formula", sep=".")]]
#   data <- if (!is.null(Call$data)) eval(Call$data)
#   else environment(formula$terms)
#   Terms <- terms(parform)
#   #XXLEV<-object[[paste(what,"xlevels",sep=".")]]
#   mf <- try(model.frame(Terms, data, xlev = object[[paste(what,"xlevels",sep=".")]]), silent = TRUE)
#   if (any(class(mf)%in%"try-error")){ mf <-  get_all_vars(parform, data)}
#   mf
# }




  
