#' @title Regressor Class
#' @description Generate design matrix with given autoregressive matrix
#' @param n sample size
#' @param p data dimension
#' @param rho autoregressive parameter
#' @export <Do you want users to be able to use this function? If not, it is a
#' @keywords
#' @seealso
#' @import mvtnorm
#' @importFrom <Load specific functions from packages as dependencies>
#' @return Object returned
#' @aliases
#' @examples
#'
#'

library(R6)

Regressor <- R6Class("Regressor",
                     list(P = NULL,
                          rho = NULL,
                          Cov = NULL,
                          initialize = function( P, rho, Cov = NULL
                          ){

                            self$P = P

                            self$rho = rho

                            self$Cov = Cov

                          },
                          normal = function(n){

                             if(is.null(self$Cov)){

                             self$Cov <- rho^abs(outer(1:P, 1:P, "-")) #correlation matrix

                             }

                             X <- rmvnorm(n, mean = rep(0,self$P), sigma = self$Cov)

                             return(X)
                           },
                          block.normal = function(n, Pblock, rho.vector, mc.cores){



                           # X <- matrix(NA, nrow = n, ncol = sum(Pblock))

                          #  for(i in 1:length(Pblock)){

                           #   sigma = rho.vector[i]^abs(outer(1:Pblock[i], 1:Pblock[i], "-"))


                          #    X[, ((i-1)*100 + 1):(100*i) ] =  rmvnorm(n, mean = rep(0,Pblock[i]), sigma = sigma)

                           # }

                            a <- Sys.time()

                            X  <- mclapply(1:length(Pblock), function(i){

                              sigma = rho.vector[i]^abs(outer(1:Pblock[i], 1:Pblock[i], "-"))

                              dat = rmvnorm(n, mean = rep(0,Pblock[i]), sigma = sigma)

                            }, mc.cores = mc.cores  ) #%>%

                            X <- do.call(cbind, X)

                            return(X)
                          },
                          block.geno = function(n, Pblock, rho.vector, mc.cores){

                            a <- Sys.time()

                            X  <- mclapply(1:length(Pblock), function(i){

                              rho = rho.vector[i]

                              m = Pblock[i]

                              n = n

                              p = 0.1

                              q = 1- p

                              normal.cop=normalCopula(rho,dim=m,dispstr="ar1")

                              u <- rCopula(n, normal.cop)

                              dat=array(qbinom(u,size=2,prob=p), dim=c(n, m))

                            }, mc.cores = mc.cores) #%>%

                            X <- do.call(cbind, X)

                            return(X)
                          },
                           SNP = function(n){

                             p= 0.2;

                             q= 1 - p; #binomial success and failure rates, could be other values;

                             normal.cop = normalCopula(self$rho, dim= self$P ,dispstr="ar1")

                             u.1 <- rCopula(n, normal.cop)

                             X <- array(qbinom(u.1,size=2,prob=p), dim=c(n, self$P))

                             return(X)
                           } )
)

