# Authors Mikis Stasinopoulos and Bob Rigby with contribution from Elaine Borghie
# last change Monday, September 20, 2004 
Q.stats <- function(obj, 
                    xvar = NULL, 
                    xcut.points = NULL, 
                    n.inter = 10,
                    zvals = TRUE,
                    save = TRUE,
                    ...)
{
# function 1 -------------------------
  check.overlap <- function(interval)
  {
    if (!is.matrix(interval) ) {stop(paste("The interval specified is not a matrix."))}
    if (dim(interval)[2] !=2) 
     {
stop(paste("The interval specified is not a valid matrix.\nThe number of columns should be equal to 2."))
     }
    crows = dim(interval)[1]
    for (i in 1:(crows-1))
    {
        if (!(abs(interval[i,2]-interval[i+1,1])<0.0001)) {interval[i+1,1]=interval[i,2]}
    }
    return(interval)
  }
#------------------------------------------
# function 2 
  get.intervals <- function (xvar, xcut.points ) 
{
    if (!is.vector(xcut.points))  {stop(paste("The interval is not a vector."))}
    if ( any((xcut.points < min(xvar)) | any(xcut.points > max(xvar))))
        {stop(paste("The specified `xcut.points' are not within the range of the x: (", min(xvar),
        " , ", max(xvar), ")"))}
    int <- c(min(xvar), xcut.points, max(xvar))
     ii <- 1:(length(int)-1)
      r <- 2:length(int)
     x1 <- int[ii]
     xr <- int[r]
     if (any(x1>xr)) {stop(paste("The interval is are not in a increasing order."))}
    cbind(x1,xr)
}
#------------------------------------------
# function 3
Qtest <- function(x)
 {
       N <- length(x)
      mx <- mean(x)
      m2 <- sum((x-mx)^2)/N
      m3 <- sum((x-mx)^3)/N
      m4 <- sum((x-mx)^4)/N
  sqrtB1 <- m3/(m2^1.5)
    yval <- sqrtB1*sqrt((N+1)*(N+3)/(6*(N-2))) 
    b2b1 <- 3*(N*N+27*N-70)*(N+1)*(N+3)/((N-2)*(N+5)*(N+7)*(N+9)) 
     wsq <- -1 +sqrt(2*b2b1-2)
   delta <- 1/sqrt(0.5*log(wsq))
   alpha <- sqrt(2/(wsq-1))
  yalpha <- yval/alpha
     zb1 <- delta *log(yalpha+sqrt(yalpha^2+1)) 
     #------zbi
      b2 <- m4/(m2*m2)
  meanb2 <- 3*(N-1)/(N+1)
   varb2 <- 24*N*(N-2)*(N-3)/((N+1)*(N+1)*(N+3)*(N+5))
    xval <- (b2-meanb2)/sqrt(varb2)
   sb1b2 <- 6*(N*N-5*N+2)/((N+7)*(N+9)) 
   sb1b2 <- sb1b2*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)))
  isb1b2 <- 1/sb1b2
       A <- 6+8*isb1b2*(2*isb1b2+sqrt(1+4*isb1b2^2))
     zb2 <- (1-2/(9*A))-((1-(2/A))/(1+xval*sqrt(2/(A-4))))^(1/3)
     zb2 <- zb2/sqrt(2/(9*A))
      K2 <- (zb1*zb1)+(zb2*zb2)
      Q1 <- N*(mx^2)
      sdy <- sqrt(m2)
      Q2  <- ((sdy^(2/3)-(1-(2/(9*N-9))))^2)/(2/(9*N-9))
      Q3  <- zb1^2
      Q4  <- zb2^2
      
      z1  <- sqrt(N)*mx
      z2  <-  sqrt(Q2)*sign(sdy^(2/3)-(1-(2/(9*N-9))))
    list(Q1=Q1,  Q2=Q2, Q3=Q3, Q4=Q4, z1=z1, z2=z2, z3=zb1, z4=zb2, AgostinoK2=K2, N=N )      
 }
#---------------------------
# function 4
get.par.df <- function(obj)
      {
         mu.df <- obj$mu.df
      sigma.df <- if ("sigma"%in%obj$parameters) (obj$sigma.df+1)/2 else 0
         nu.df <- if ("nu"%in%obj$parameters & !(obj$family[1]%in%c("TF","PE"))) obj$nu.df else 0
        tau.df <- if ("tau"%in%obj$parameters & !(obj$family[1]%in%c("TF","PE"))) obj$tau.df 
                  else if ((obj$family[1]%in%c("TF","PE"))) obj$nu.df
                  else 0
        maxall <- max(c(mu.df,sigma.df,nu.df,tau.df)) 
      list(location=mu.df, scale=sigma.df, skewness=nu.df, kurtosis=tau.df, maximum=maxall )             
      }
# main function starts here
    overlap <- 0
    if (is.gamlss(obj))  
      {
      var <- resid(obj)
      df <- get.par.df(obj)
      } 
    else
      {
      stop("this is not a GAMLSS object /n")
      #   var <- obj  
      }              
    if(is.null(xvar)) 
      { xvar<-seq(1,length(obj$y),1)
        warning(paste("The xvar has been replace by an index 1:n", "\n", ""))  
      }                  
    if(is.null(xcut.points)) 
      { # getting the intervals automatic
      if (n.inter < df$maximum+1.9) 
         {
          n.inter <- ceiling(df$maximum+1.9)
          warning("the number of intervals have change to ", n.inter,"\n ")
         } 
      g.in <- co.intervals(xvar, number=n.inter, overlap=overlap)
      if (overlap==0) g.in <- check.overlap(g.in) 
      }                 
    else
      {
        if (length(xcut.points)+1 < df$maximum+1.9) 
         {
          stop("the number of cut points must be larger than ", ceiling(df$maximum+1.9),"\n ")
         } 
       # if xcut.points is set
      g.in <- get.intervals(xvar, xcut.points ) 
      }
     howmany <- dim(g.in)[1]
           X <- matrix(0, nrow = howmany, ncol = 6, 
                       dimnames=list(as.character(seq(1,howmany)),
                                     as.character(seq(1,6))))
           Z <- matrix(0, nrow = howmany, ncol = 6, 
                       dimnames=list(as.character(seq(1,howmany)),
                                     as.character(seq(1,6))))                         
     dimnames(Z)[[2]] <- as.character(c("Z1","Z2","Z3","Z4","AgostinoK2","N"))                              
     dimnames(X)[[2]] <- as.character(c("Q1","Q2","Q3","Q4","AgostinoK2","N"))
    dimnames(X)[[1]] <- paste(substr(as.character(g.in[1:howmany,1]),1,7),"to",substr(as.character(g.in[1:howmany,2]),1,7))
    dimnames(Z)[[1]] <- paste(substr(as.character(g.in[1:howmany,1]),1,7),"to",substr(as.character(g.in[1:howmany,2]),1,7))
       oxvar <- xvar[order(xvar)]
       oyvar <- var[order(xvar)] 
   for (i in 1:howmany)
       {  
          if(i==howmany) {             ##### Include points at the end of the last interval
                          yvar1 <- subset(oyvar, oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
                          xvar1 <- subset(oxvar, oxvar>=g.in[i,1]&oxvar<=g.in[i,2])
                          }
          else
          {
          yvar1 <- subset(oyvar, oxvar>=g.in[i,1]&oxvar<g.in[i,2])
          xvar1 <- subset(oxvar, oxvar>=g.in[i,1]&oxvar<g.in[i,2]) 
          }        
          nlist <- Qtest(yvar1)
           Z[i,] <- c(nlist$z1, nlist$z2, nlist$z3, nlist$z4, nlist$AgostinoK2, nlist$N)            
           X[i,] <- c(nlist$Q1, nlist$Q2, nlist$Q3, nlist$Q4, nlist$AgostinoK2, nlist$N)
       }
          Q1 <- sum(X[,1])
          Q2 <- sum(X[,2])
          Q3 <- sum(X[,3])
          Q4 <- sum(X[,4])
          K2 <- sum(X[,5])
          N  <-sum(X[,6])
         df1 <- howmany-df$location
         df2 <- howmany-df$scale
         df3 <- howmany-df$skewness
         df4 <- howmany-df$kurtosis
         df5 <- df3+df4
         pv1 <- pchisq(Q1, df=df1, ncp=0, lower.tail = FALSE)
         pv2 <- pchisq(Q2, df=df2, ncp=0, lower.tail = FALSE)
         pv3 <- pchisq(Q3, df=df3, ncp=0, lower.tail = FALSE)
         pv4 <- pchisq(Q4, df=df4, ncp=0, lower.tail = FALSE)
         pv5 <- pchisq(K2, df=df5, ncp=0, lower.tail = FALSE)
          nc <- matrix(c(Q1,Q2,Q3,Q4,K2,N),nrow=1, dimnames=list(as.character(seq(1,1)),
                                     as.character(seq(1,6))) )
          dimnames(nc)[[1]]<-as.character(c("TOTAL Q stats"))
         ndf <- matrix(c(df1,df2,df3,df4,df5,0),nrow=1, dimnames=list(as.character(seq(1,1)),
                                     as.character(seq(1,6))) )
          dimnames(ndf)[[1]]<-as.character(c("df for Q stats"))
         pval <- matrix(c(pv1,pv2,pv3,pv4,pv5,0),nrow=1, dimnames=list(as.character(seq(1,1)),
                                     as.character(seq(1,6))) )
          dimnames(pval)[[1]]<-as.character(c("p-val for Q stats"))
          X <- rbind( X , nc ,ndf, pval)
          Z <- rbind( Z , nc ,ndf, pval)      
if (save) {
           if (zvals) return(Z)  
           else return(X)
           }
}      
          
