#' QR2 Decomposition
#'
#' @param A, an n x m matrix
#' @return a list with Q.transpose and R
#' @export
#' @examples
myQR <- function(A){

  ## Perform QR decomposition on the matrix A
  ## Input:
  ## A, an n x m matrix

  ########################
  ## FILL IN CODE BELOW ##
  ########################

  R = A
  n = dim(A)[1]
  m = dim(A)[2]

  Q = diag(rep(1,n))

  for(k in 1:(m-1))
  {
    X = matrix(0,n,1)
    X[k:n,1] = R[k:n,k]
    V = X
    V[k] = X[k] + sign(X[k,1])*norm(X, type = "F")
    S = norm(V, type = "F")
    u = V/S
    R = R - 2 * (u %*% (t(u) %*% R))
    Q = Q - 2 * (u %*% (t(u) %*% Q))
  }
  ## Function should output a list with Q.transpose and R
  ## Q is an orthogonal n x n matrix
  ## R is an upper triangular n x m matrix
  ## Q and R satisfy the equation: A = Q %*% R
  return(list("Q" = t(Q), "R" = R))

}


#' myLM Decomposition
#'
#' @param X is an n x p matrix of explanatory variables, Y is an n dimensional vector of responses
#' @return Function returns the beta_ls, the least squares and standard err
#' @export
#' @examples
myLinearModel <- function(X, Y){

  ## Perform the linear regression of Y on X
  ## Input:
  ## X is an n x p matrix of explanatory variables
  ## Y is an n dimensional vector of responses
  ## Do NOT simulate data in this function. n and p
  ## should be determined by X.
  ## Use myQR inside of this function

  ########################
  ## FILL IN CODE BELOW ##
  ########################
  n <- nrow(X)
  p <- ncol(X)

  Z = cbind(rep(1,n), X, Y)
 # A = t(Z) %*% Z
  D = myQR(Z)

  R1 = D$R[1:(p+1), 1:(p+1)]
  Y1 = D$R[1:(p+1), p+2]

  beta_ls = solve(R1) %*% Y1
  X = cbind(rep(1,n), X)
  sig = sum((Y-(X %*% beta_ls))^2)/(n-p-1)
  err = sqrt(diag(sig*solve(t(X) %*% X)))
  
  
  output = list(beta_ls = beta_ls, std_err = err)
  
  ## Function returns the 1 x (p + 1) vector beta_ls,
  ## std error
  return(output)

}





