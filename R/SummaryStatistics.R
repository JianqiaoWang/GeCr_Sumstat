#' @title Generate Summary Statistics from generated X and Y.
#' @description From given X and Y, generate summary statistics, t score in this case.
#' @param X Design matrix
#' @param Y Response vector
#' @param adj.C adjusted for Covariates
#' @export <Do you want users to be able to use this function? If not, it is a
#' @keywords
#' @seealso
#' @import <Load entire package as dependency>
#' @importFrom <Load specific functions from packages as dependencies>
#' @return t statistics returned
#' @aliases
#' @examples

SumStat<- function(Y, X, adj.C = NULL)
{

  Y <- as.matrix(Y);
  n <- length(Y);
  X <- as.matrix(X);

  if(!is.null(adj.C)){
  adj.C <- cbind(1,as.matrix(adj.C));
  U1 <- crossprod(adj.C,Y);
  U2 <- solve(crossprod(adj.C),U1);
  ytr <- Y-adj.C%*%U2;
  U3 <- crossprod(adj.C,X);
  U4 <- solve(crossprod(adj.C),U3);
  Xtr <- X-adj.C%*%U4;
  }else{

  ytr = scale(Y, center = T, scale = F)

  Xtr = scale(X, center = T, scale = F)

  }

  Xtr2 <- colSums(Xtr^2);

  b <- as.vector(crossprod(ytr, Xtr)/Xtr2 );

  ncol.X = ifelse(is.null(adj.C), 1, ncol(adj.C))

  sig <- (sum(ytr^2)-b^2*Xtr2)/(n-ncol.X-1);

  err <- sqrt(sig*(1/Xtr2));
  ## return z-scores
  return(b/err);
} # Calculate the t scores.



z.to.r <- function(z, n){

  return(z /sqrt(z^2 + n - 2))

}


r.to.z <- function(r, n){

  return(r * sqrt(n-2) /sqrt(1- r^2))

}

