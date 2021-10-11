
library(GLMsData)
data("quilpie")
print(quilpie)

matqui <- as.matrix(quilpie)
xsoi <- matqui[,3]
xsoi <- as.numeric(xsoi)

quilpie$Phase <- factor( quilpie$Phase )
xphase <- with( quilpie, model.matrix( ~ Phase ) )


y <- matqui[,6]
y <- as.numeric(y)

soi.quilpie <- FitModelMle(y, xsoi)
LLH <- soi.quilpie$LLH
soi.aic <- -2 * LLH + 2 * length(soi.quilpie$coef)
m1.bic <- -2 * LLH + log(length(y)) * length(soi.quilpie$coef)

phase.quilpie <- FitModelMle(y, xphase, maxits = 15, add.constant=FALSE )
LLH <- phase.quilpie$LLH
phase.aic <- -2 * LLH + 2 * length(phase.quilpie$coef)
m2.bic <- -2 * LLH + log(length(y)) * length(phase.quilpie$coef)

phasesoi <- glm(y ~ xsoi+xphase, family = binomial(link = "logit"))
phasesoi.aic <- phasesoi$aic

aic<-list(soi.aic,
          phase.aic, 
          phasesoi.aic)


#Source: Dun & Smyth, 2017. Page 203-204

# Function for computing the information matrix:
MakeExpInf <- function(x, mu){
  # Args:
  # x: The matrix of explanatory variables
  # mu: The fitted values
  #
  # Returns:
  # The expected information matrix
  if ( length(mu) == 1 ) mu <- rep( mu, dim(x)[1])
  mu <- as.vector(mu)
  return( t(x) %*% diag( mu * (1 - mu) ) %*% x )
}
# Function for computing mu:
MakeMu <- function(x, beta){
  # Args:
  # x: The matrix of explanatory variables
  # beta: The linear model parameter estimates
  #
  # Returns:
  # The value of mu
  eta <- x %*% beta
  return( 1 / ( 1 + exp( -eta ) ) )
}
# Function for computing the score vector:
MakeScore <- function(x, y, beta){
  # Args:
  # x: The matrix of explanatory variables
  # y: The response variable
  # beta: The linear model parameter estimates
  #
  # Returns:
  # The score matrix
  mu <- MakeMu(x, beta)
  return( t(x) %*% (y - mu) )
}
FitModelMle <- function(y, x=NULL, maxits=8, add.constant=TRUE){
  # Args:
  # y: The response variable
  # x: The matrix of explanatory variables
  # maxits: The maximum number of iteration for the algorithm
  # add.constant: If TRUE, a constant is added to the x matrix
  # (All models must have a constant term.)
  #
  # Returns:
  # Information about the fitted glm
  if ( is.null(x)){ # If no x given, ensure constant appears
    allx <- cbind( Constant=rep( 1, length(y) ) )
  } else {
    allx <- x
    if( add.constant ){
      allx <- cbind( Constant=rep(1, length(y)), x)
    }
  }
  num.x.vars <- dim(allx)[2] - 1 # Take one, because of constant
  # Find initials: beta_0 = mean(y), and the other beta_j are zero
  beta <- c( mean(y), rep( 0, num.x.vars ) )
  
  # Set up
  beta.vec <- array( dim=c(maxits, length(beta) ) )
  beta.vec[1,] <- beta
  mu <- MakeMu( allx, beta )
  score.vec <- MakeScore(allx, y, beta)
  inf.mat <- MakeExpInf( allx, mu )
  # Now iterate to update
  for (i in (2:maxits)){
    beta <- beta + solve( inf.mat ) %*% score.vec
    beta.vec[i,] <- beta
    mu <- MakeMu( allx, beta )
    score.vec <- MakeScore(allx, y, beta)
    inf.mat <- MakeExpInf( allx, mu )
  }
  # Compute log-likelihood
  LLH <- sum( y*log(mu) + (1-y)*log(1-mu) )
  return( list(coef = beta.vec[maxits,], # MLE of parameter estimates
               coef.vec = beta.vec, # Estimates at each iteration
               LLH = LLH, # The maximum log-likelihood
               inf.mat = inf.mat, # The information matrix
               score.vec = score.vec, # The score vector
               mu = mu) ) # The fitted values
}