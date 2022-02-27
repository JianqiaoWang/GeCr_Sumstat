

Tracel2 = function (Sigma, Omega)
{
  return(sum(diag((Sigma %*% Omega - diag(1, dim(Omega)[1]))^2)))
}


 FASTCLIME = function(X){

   library(fastclime)

   S = cor(X)

   P = nrow(S)

   lambda.min.ratio = 0.2

   lambda.max.tmp1 = min(max(S - diag(diag(S))), -min(S -
                                                        diag(diag(S))))

   lambda.max.tmp2 = max(max(S - diag(diag(S))), -min(S -
                                                        diag(diag(S))))
   if (lambda.max.tmp1 == 0){

     lambda.max = lambda.max.tmp2

   }else {

     lambda.max = lambda.max.tmp1

   }

   #lambda.max = 2 * lambda.max

   lambda.range <- exp(seq(log(lambda.min.ratio*lambda.max), log(2*lambda.max), length.out = 40))

   lambda.opt = CV(X = X, fold = 2, lambda.range = lambda.range)

   lambda <- lambda.range[lambda.opt]

   gc()

   print(lambda)

   # Finds the estimated path corrsponding to lambda=0.2

   Omega.final.1 <- fastclime(X, lambda.min = 0.9*lambda, nlambda = 150)

   Omega.final.2 <- fastclime.selector(Omega.final.1$lambdamtx, Omega.final.1$icovlist, lambda)

   return(Omega.final.2$icov)
 }

CV = function(X, fold = 2, lambda.range){

  N = nrow(X)

  N.train <- ceiling(N/fold)

  N.valid <- floor(N/fold)

  index <- sample(1:N, size = N.train, replace = FALSE)

  cv.index <- list(index, c(1:N)[-index])

  scalar = 1 - 1/N.train

  loss.re = matrix(0, nrow = fold, ncol = length(lambda.range))

  for (i in 1:fold) {

    X.train <- X[cv.index[[i]],]

    X.valid = X[-cv.index[[i]],]

    out1 = fastclime(X.train, lambda.min = min(lambda.range), nlambda = 150) # Estimate the solution path

    for( j in 1:length(lambda.range)){

      out2 = fastclime.selector(out1$lambdamtx, out1$icovlist, lambda.range[j])

      loss.re[i,j] <- Tracel2(Sigma = scalar*cov(X.valid), Omega = out2$icov)

    }

  }

  loss.mean = apply(loss.re, 2, mean)

  opt.idx = which.min(loss.mean)

  return(opt.idx)
}



Flare = function(X.design, fold = 5, nlambda = 5){

  library(flare)

  out1 = sugm(X.design, method = "clime", nlambda = nlambda)

  out1.select2 = sugm.select(out1, criterion = "cv", fold = fold, rep.num = 5)

  return(out1.select2$opt.icov)

}



Hardthresholding = function(X, fixed.lambda){

  S = cor(X)

  Temp <- diag(S)

  S[(abs(S) < fixed.lambda)] <- 0

  diag(S) <- Temp

  return(S)

}


inv.thresh = function(X){

 # library(CVTuningCov)

  n <- nrow(X);

  p <- ncol(X);

  k.grid <- sqrt(log(p)/n) * seq(0,20,by=0.05);

  CV.F.fit <- regular.CV(X, k.grid = k.grid, method = "HardThresholding",
             fold = 2, norm = "F", seed=10323)

  threvalue = CV.F.fit$k.grid[which.min(rowSums(CV.F.fit$CV.pre.error))]

  hard.Sigma<-Hardthresholding(X, fixed.lambda = threvalue);

  hard.sigma.decom <- eigen(hard.Sigma, symmetric = T)

  hard.sigma.decom$values <- pmax(hard.sigma.decom$values,1 - 2*sqrt(log(p)/N))

  invSigma <- hard.sigma.decom$vectors %*% diag(1/hard.sigma.decom$values) %*% t(hard.sigma.decom$vectors)

  return(invSigma)

}



thresh.loss <- function(mat1, mat2, threshold, method, norm) {

  cov.mat1 <- cov(mat1)

  cov.mat2 <- cov(mat2)

  if (method == "hard") {
    mat.diff <- hardthresholding(cov.mat1, threshold) -
      cov.mat2
  }

  if (norm %in% c("F", "f")) {
    loss <- F.norm2(mat.diff)
  }

  else if (norm %in% c("O", "o")) {

    loss <- O.norm2(mat.diff)

  }

  return(loss)
}


