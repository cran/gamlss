#------------------------------------------------------------------------------------------
# Mikis Stasinopoulos
# 7 - 4 - 11
# This is a new version wihch can take any object with resid()
# or the actual residuals using the argument resid
##  worm plots based on the worm plot function in S-plus by Van Buuren and 
##  Fredriks (2001)  Worm plot: a simple diagnostic device for modelling 
## growth reference curves, Statist. Med.
##  original created by Mikis Stasinopoulos  and Bob Rigby, Friday, November 8, 2002.
##  last change Tuersday 
wp  <- function (object = NULL, 
		           xvar = NULL,
	               resid = NULL,
		        n.inter = 4,
		    xcut.points = NULL, 
		        overlap = 0,
		       xlim.all = 4,
		      xlim.worm = 3.5, 
		     show.given = TRUE,
		           line = TRUE,
		       ylim.all = 12*sqrt(1/length(resid)),
		      ylim.worm = 12*sqrt(n.inter/length(resid)),   
		            cex = 1, 
		            pch = 21,
		                  ...)
{
	##-----------------------------------------------------------------------------
	panel.fun<-function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
			cex = par("cex"), col.smooth = "red", span = 2/3, iter = 3, ...) 
	{
		qq <- as.data.frame(qqnorm(y, plot = FALSE))
		qq$y <- qq$y - qq$x    
		grid(nx=NA,ny=NA, lwd = 2) # grid only in y-direction 
		#     points(qq$x, qq$y, pch = 1, col = col, bg = bg, cex = 0.7)
		points(qq$x, qq$y, pch = pch, col = col, bg = bg, cex = cex)
		abline(0, 0, lty = 2, col = 1)
		abline(0, 100000, lty = 2, col = 1)
		yuplim <- 10*sqrt(1/length(qq$y))
		level <- .95    
		lz <- -xlim.worm
		#  lz <- -2.75
		hz <-  xlim.worm 
		#  hz <-  2.75 
		dz <- 0.25
		z <- seq(lz,hz,dz)
		p <- pnorm(z)   
		se <- (1/dnorm(z))*(sqrt(p*(1-p)/length(qq$y)))     
		low <- qnorm((1-level)/2)*se
		high <- -low      
		# if ( any(abs(qq$y)>ylim.worm) ) 
		{
			no.points <- length(qq$y)
			total.points<<-total.points+no.points
			no.mis<-sum(abs(qq$y)>ylim.worm)
			cat("number of missing points from plot=", no.mis," out of ",no.points,"\n") 
			if ( any(abs(qq$y)>ylim.worm) )  warning("Some points are missed out ", "\n","increase the y limits using ylim.worm" )
		}
		if ( any(abs(qq$x)>xlim.worm) ) 
		{
			warning("Some points are missed out ", "\n","increase the x limits using xlim.worm" )
		}
		# se <- (1/dnorm(z))*(sqrt(p*(1-p)*n.inter/length(xvar)))    
		# low <- qnorm((1-level)/2)*se
		# high <- -low
		lines(z, low,  lty=2,  lwd=0.01)
		lines(z, high, lty=2,  lwd=0.01)
		if(line == TRUE)
		{
			fit <- lm(qq$y ~ qq$x+I(qq$x^2)+I(qq$x^3)) #poly(qq$x,3)) 
			s <- spline(qq$x, fitted(fit))
			flags <- s$x > -2.5 & s$x < 2.5
			lines(list(x = s$x[flags], y = s$y[flags]), col="red", lwd=0.01)
			assign("coef1",coef(fit),envir = parent.frame(n=3))
			assign( "coef" , c(coef, coef1),envir = parent.frame(n=3))
		}
	}                
	##-----------------------------------------------------------------------------
	check.overlap <- function(interval)
	{
		if (!is.matrix(interval) ) {stop(paste("The interval specified is not a matrix."))}
		if (dim(interval)[2] !=2) {stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))}
		crows = dim(interval)[1]
		for (i in 1:(crows-1))
		{
			#if (interval[i,2] != interval[i+1,1]) {interval[i+1,1]=interval[i,2]}
			if (!(abs(interval[i,2]-interval[i+1,1])<0.0001)) {interval[i+1,1]=interval[i,2]}
		}
		return(interval)
	}
	##-----------------------------------------------------------------------------
# the problem here is that coplot uses  for given.values      MS Tuesday, September 21, 2004 
#  min <=  x <= 20                   min <=  x < 20             min <=  x <= 20-extra  
#   20 <=  x <= 30    instead   of    20 <=  x < 30  here uses   20 <=  x <= 30-extra  
#   30 <=  x <= max                   30 <=  x <= max            30 <=  x <= max+extra 
	get.intervals <- function (xvar, xcut.points ) 
	{
		if (!is.vector(xcut.points))  {stop(paste("The interval is not a vector."))}
		if ( any((xcut.points < min(xvar)) | any(xcut.points > max(xvar))))
		{stop(paste("The specified `xcut.points' are not within the range of the x: (", min(xvar),
							" , ", max(xvar), ")"))}
		extra <-(max(xvar)-min(xvar))/100000
		int <- c(min(xvar), xcut.points, (max(xvar)+2*extra))
		ii <- 1:(length(int)-1)
		r <- 2:length(int)
		x1 <- int[ii]
		xr <- int[r]-extra
		if (any(x1>xr)) {stop(paste("The interval is are not in a increasing order."))}
		cbind(x1,xr)
	}
	##-here is the main function ------------------------------------------------------------
	if (is.null(object)&&is.null(resid))  stop(paste("A fitted object with resid() method or the argument resid should be used ", "\n", ""))
	resid <- if (is.null(object)) resid else resid(object)  
	if(is.null(xvar)) 
	{
		qq   <- as.data.frame(qqnorm(resid, plot = FALSE))
		qq$y <- qq$y - qq$x
		level <- .95    
		lz <- -xlim.all
		hz <- xlim.all 
		dz <- 0.25
		z <- seq(lz,hz,dz)
		p <- pnorm(z)   
		se <- (1/dnorm(z))*(sqrt(p*(1-p)/length(qq$y)))    
		low <- qnorm((1-level)/2)*se
		high <- -low
		if ( any(abs(qq$y)>ylim.all) ) 
		{
			warning("Some points are missed out ", "\n","increase the y limits using ylim.all" )
		}
		if ( any(abs(qq$x)>xlim.all) ) 
		{
			warning("Some points are missed out ", "\n","increase the x limits using xlim.all" )
		}
		plot(qq$x,qq$y,ylab="Deviation", xlab="Unit normal quantile", 
				xlim=c(-xlim.all,xlim.all), ylim=c(-ylim.all , ylim.all), cex=cex, pch=pch ,   bg = "wheat", )
		grid(lty = "solid")
		abline(0, 0, lty = 2, col = 2)
		abline(0, 100000, lty = 2, col = 2)
		lines(z, low, lty=2)
		lines(z, high, lty=2)
		if(line == TRUE)
		{
			fit <- lm(qq$y ~ qq$x+I(qq$x^2)+I(qq$x^3)) #poly(qq$x,3)) 
			s <- spline(qq$x, fitted(fit))
			flags <- s$x > -3 & s$x < 3
			lines(list(x = s$x[flags], y = s$y[flags]), col="red", lwd=0.01)
		}
	}   
	else   
	{
		w  <- if (is.null(object)) rep(1, length(resid)) else  object$weights # geting the weights
		if (all(trunc(w)==w))  xvar <- rep(xvar, w) # if w=0 reduce if w=freq expand if all(w=1) same 
		if(is.null(xcut.points)) 
		{ # getting the intervals automatic
			given.in <- co.intervals(xvar, number=n.inter, overlap=overlap )
			if (overlap==0) given.in <- check.overlap(given.in) 
		}                 
		else
		{ # if xcut.points is set
			given.in <- get.intervals(xvar, xcut.points )
		}
		total.points <-0
		coef<-coef1<-NULL
		y <- resid   
		x <- resid  
		coplot(y~x | xvar , 
				given.values = given.in, 
				panel = panel.fun,
				ylim = c(-ylim.worm, ylim.worm), 
				xlim = c(-xlim.worm, xlim.worm),
				ylab = "Deviation", 
				xlab = "Unit normal quantile",
				show.given = show.given,
				bg = "wheat", 
				pch = pch,
				cex = cex, 
				bar.bg = c(num="light blue")) 
		if (overlap==0)
		{
			if  (total.points!=length(resid))
				warning("the total number of points in the plot is not equal \n to the number of observations in y \n")
		}
	}
	if(!is.null(xvar)){  mcoef <-matrix(coef,ncol=4,byrow=TRUE)        
		out<-list(classes = given.in, coef = mcoef)}
}
