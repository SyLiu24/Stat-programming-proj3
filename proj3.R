# Shuying Liu
# s2436365


# Code to fit x, y data using pspline smoother
# with GCV smoothing parameter selection.

# The challenge is to smooth x, y data,
# i.e. estimate the smooth function f in model: y=f(x)+e.
# We will use a pspline fit to estimate f: Use basis expansion
# f=sum(beta_jb_j) where b_i are k evenly spaced Bspline basis functions, so the
# model will be turned into a linear model y=Xbeta+e, where X_ij=b_j(x_i),
# and beta is the vector of coefficients of basis functions.
# And we impose a penalty, for example order 2 difference
# sum(beta_{i-1}-2beta_i+beta_{i+1})=beta^TD^TDbeta
# to avoid over-fitting, then estimate the linear model by penalized least squares 
# min |y-Xbeta|^2+lambda*penalty
# where lambda is the smoothing parameter.
# So the estimated beta^hat=((X^TX+lambda*D^TD)^{-1}X^Ty, fitted value=Xbeta^hat.
# For the best fit smoother, the smoothing parameter lambda will be selected by 
# finding the value of λ to minimize the generalized cross validation criterion
# GCV = estimated residual variance (sig2) / (n - effective degrees of freedom),
# where edf=tr{(X^TX+lambda*D^TD)^{-1} X^TX}, sig2=|y-Xbeta|^2/(n-edf).


pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  # Find the best pspline fit smoother for given x, y data by
  # computing GCV for many trial values of lambda in a given range
  # then finding the trial value of lambda that minimizes GCV.
  
  # Input:
  #   x,y – the vectors of x, y data to smooth; k - the number of basis functions;
  #   logsp – the ends of the interval to search for lambda (log lambda scale);
  #   bord – the B-spline order to use;
  #   pord – the order of difference to use in the penalty;
  #   ngrid – the number of smoothing parameter values to try
  
  # Return a list of class pspline defining the best fit spline smoother
  
  # Set up the basis
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  # Matrix for the penalty beta^TD^TDbeta
  D <- diff(diag(k),differences=pord)
  
  n <- length(x)
  
  # Perform QR decomposition X=QR 
  # and then eigen decomposition U*Lambda*U^T=R^{-T}D^TDR^{-1}
  # for more efficient matrix computation
  # beta^hat=R^{-1}U(I + lambda*Lambda)^{−1}U^TQ^Ty
  # edf=tr{(I+lambda*Lambda)^{-1}}
  qrx <- qr(X)
  Q <- qr.Q(qrx)
  R <- qr.R(qrx)
  qty=qr.qty(qrx,y)[1:k]
  
  eRD <- tcrossprod(forwardsolve(t(R),t(D))) |> eigen()
  U <- eRD$vectors
  Lambda <- eRD$values
  
  # Generate lambda trial values in the given range
  if (length(logsp) != 1) {
    lambda <- seq(logsp[1],logsp[2],length=ngrid) |> exp()
  }
  else lambda <- exp(logsp)
  
  
  coef <- matrix(0,k,ngrid) # estimated beta
  edf <- rep(0,ngrid) # effective degrees of freedom
  
  # Compute estimated beta and edf for all trial values of lambda
  RU <- backsolve(R,U)
  UQy <- t(U) %*% qty
  for (i in 1:ngrid) {
    diag_lambda <- 1/(1+lambda[i]*Lambda)
    edf[i] <- sum(diag_lambda)
    coef[,i] <- RU %*% (diag_lambda * UQy)
  }
  # Compute model fits for all trail values of lambda
  fitted <- X %*% coef # fitted values
  sig2 <- colSums((y-fitted)^2)/(n-edf) # estimated residual variance
  gcv <- sig2/(n-edf) # generalized cross validation

  # Find the lambda minimizing GCV to get the best fit
  lambda_i <- which.min(gcv)
  diag_lambda_fit <- 1/(1+lambda[lambda_i]*Lambda)
  
  # Get details about the best fit
  sig2_fit <- sig2[lambda_i] # estimated residual variance
  V_fit <- sig2_fit*RU %*% (diag_lambda_fit * t(RU)) # covariance matrix of beta
  r2_fit <- 1-(n-1)*sig2_fit/sum((y-mean(y))^2) # r-squared
  
  smoother <- list(data=cbind(x,y),k=k,ord=c(bord,pord),
                   knots=knots,coef=coef[,lambda_i],
                   fitted=fitted[,lambda_i],sig2=sig2_fit,edf=edf[lambda_i],
                   r2=r2_fit,gcv=gcv[lambda_i],cov=V_fit)
  class(smoother) <- "pspline"
  
  return(smoother)
}


print.pspline <- function(m) {
  # Report some details of the model fit m
  # Return a list containing gcv, edf and r2 of m.
  
  l1 <- sprintf("Order %d p-spline with order %d penalty",
                m$ord[1],m$ord[2])
  l2 <- sprintf("Effective degrees of freedom: %.7g\tCoefficients: %d",
                m$edf,m$k)
  l3 <- sprintf("residual std dev: %.7g\tr-squared: %.7g\tGCV: %.7g",
                m$sig2^.5,m$r2,m$gcv)
  cat(l1,l2,l3,sep='\n')
  
  invisible(list(gcv=m$gcv,edf=m$edf,r2=m$r2))
}


predict.pspline <- function(m,x,se=TRUE) {
  #  Make predictions from the smooth fit m for new x values
  # (within the range of the original data)
  # If se=FALSE, function returns the predicted values;
  # If se=TRUE, function returns a list containing the predicted values 'fit'
  # and the corresponding standard errors 'se'
  
  # Generate new X model matrix
  Xp <- splines::splineDesign(m$knots,x,ord=m$ord[1]+1,outer.ok=TRUE)
  fitted <- Xp %*% m$coef #new fitted values

  if (se) {
    # The  required standard errors come from the
    # square roots of the leading diagonals 
    # of the predicted value covariance matrix X_pVX_p^T,
    # where V is  the covariance matrix for the coefficients
    return(list(fit=fitted,se=rowSums(Xp * (Xp %*% m$cov))^.5))
  }
  else return(fitted)
}


plot.pspline <- function(m) {
  # Produce 3 plots:
  #   First plot is the plot of original x,y data,
  #   overlaid with the estimated smooth function and the
  #   approximate 95% credible intervals for the smooth
  #   Second plot is the model residuals against fitted values
  #   Third plot is the qqplot of the residuals
  # Return a list containing the lower and upper confidence limits ll and ul
  # and the corresponding x values
  
  x <- m$data[,1]
  y <- m$data[,2]
  
  # Generate appropriate even spaced x test values 
  # and the corresponding predicted y to plot estimated smooth function
  xp <- seq(from=min(x),to=max(x),length=500)
  predict <- predict.pspline(m,xp)
  yp <- predict$fit
  se <- predict$se
  
  # Credit interval limits
  ll=yp-qnorm(0.975)*se
  ul=yp+qnorm(0.975)*se
  
  # First plot
  plot(x,y,main='Pspline Smoothing')
  lines(xp,yp)
  lines(xp,ll,lty=2)
  lines(xp,ul,lty=2)
  legend('bottomright',
         legend=c('Estimated smooth function','95% CI'),lty=c(1,2))
  
  # Second plot
  plot(m$fitted,y-m$fitted,xlab='Fitted values',
       ylab='Residuals',main='Residuals vs Fitted')
  abline(0,0,lty=2)
  
  # Third plot
  qqnorm(y-m$fitted)
  qqline(y-m$fitted)
  
  invisible(list(ll=ll,ul=ul,x=xp))
}


test_example <- function() {
  # Test with data from MASS mcycle
  
  library(MASS)
  x <- mcycle[,1]
  y <- mcycle[,2]
  
  m <- pspline(x,y)
  
  print.pspline(m)
  
  plot.pspline(m)
}