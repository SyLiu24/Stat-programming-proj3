library(MASS)
x <- mcycle[,1]
y<- mcycle[,2]


pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  D <- diff(diag(k),differences=pord)
  
  n <- length(x)
  
  qrx <- qr(X)
  Q <- qr.Q(qrx)
  R <- qr.R(qrx)
  qty=qr.qty(qrx,y)[1:k]
  
  eRD <- tcrossprod(forwardsolve(t(R),t(D))) |> eigen()
  U <- eRD$vectors
  Lambda <- eRD$values
  
  if (length(logsp != 1)) {
    lambda <- seq(logsp[1],logsp[2],length=ngrid) |> exp()
  }
  else lambda <- exp(logsp)
  

  coef <- matrix(0,k,ngrid)
  fitted <- matrix(0,n,ngrid)
  edf <- rep(0,ngrid)
  for (i in 1:ngrid) {
    diag_lambda <- 1/(1+lambda[i]*Lambda)
    
    coef[,i] <- backsolve(R,U %*% (diag_lambda * (t(U) %*% qty)))
    edf[i] <- sum(diag_lambda)
    
  }
  fitted <- X %*% coef
  sig2 <- colSums((y-fitted)^2)/(n-edf)
  gcv <- sig2/(n-edf)
  print(gcv)
  lambda_i <- order(gcv)[1]
  
  temp <- backsolve(R,U %*% ((1/(1+lambda[lambda_i]*Lambda)) * (t(U)))) %*% solve(t(R))
  V <- sig2[lambda_i]*temp
  
  r2 <- 1-(n-1)*sig2[lambda_i]/sum((y-mean(y))^2)
  
  smoother <- list(data=cbind(x,y),ord=c(bord,pord),k=k,knots=knots,coef=coef[,lambda_i],
                   fitted=fitted[,lambda_i],sig2=sig2[lambda_i],edf=edf[lambda_i],
                   rsd=sqrt(sig2[lambda_i]),r2=r2,gcv=gcv[lambda_i],cov=V)
  class(smoother) <- "pspline"
  smoother
}

m <- pspline(x,y)
plot(x,y)
plot(x,m$fitted)

print.pspline <- function(m) {
  cat()
}

predict.pspline <- function(m,x,se=TRUE) {
  Xp <- splines::splineDesign(m$knots,x,ord=m$ord[1]+1,outer.ok=TRUE)
  fitted <- Xp %*% m$coef
  if (se) {
    return(list(fit=fitted,se=rowSums(Xp*(Xp%*%m$cov))^.5))
  }
  else return(fitted)
}

plot.pspline <- function(m) {
  x <- m$data[,1]
  y <- m$data[,2]
  print(min(x))
  xp <- seq(from=min(x),to=max(x),length=300)
  predict <- predict.pspline(m,xp)
  yp <- predict$fit
  se <- predict$se
  print(yp)
  
  ll=yp-qnorm(0.975)*se
  ul=yp+qnorm(0.975)*se
  
  plot(x,y)
  lines(xp,yp)
  lines(xp,ll,lty=2)
  lines(xp,ul,lty=2)
  
  # lo <- loess(y-m$fitted~m$fitted)
  plot(m$fitted,y-m$fitted)
  abline(0,0,lty=2)
  
  qqnorm(y-m$fitted)
  qqline(y-m$fitted)
  
  return(list(ll=ll,ul=ul,x=xp))
  
}
plot.pspline(m)