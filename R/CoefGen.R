#' @title Generate coefficients
#' @description Generate coefficients
#' @param p  dimension
#' @param MAG magnitude
#' @export <Do you want users to be able to use this function? If not, it is a
#' @keywords
#' @seealso
#' @import mvtnorm
#' @importFrom <Load specific functions from packages as dependencies>
#' @return Object returned
#' @aliases
#' @examples
#'

Coef <- setRefClass("Coefficient",
                    fields = list(p = "numeric",
                                  Beta = "vector", Gamma1 = "vector"),
                    methods = list(
                      sparse = function(s){
                        Beta <<- c(rep(1,s),rep(0,p-s))
                        #return(Beta)
                      },
                      dense = function(beta.abs){
                        if(beta.abs == 1){
                          Beta <<- abs(c(round(rnorm(P,0,1),4)))
                        }else{
                          Beta <<- (c(round(rnorm(P,0,1),4)))
                        }
                        #return(Beta)
                      },
                      Coefs = function(SIGMA, h2.X, h2.Z, GeCovaraince, Var.X = 1, Var.Z = 1){

                        Beta <<- as.vector(Beta)

                        #str(Beta)

                       # str(SIGMA)

                        Beta <<- Beta*sqrt(h2.X*Var.X)/sqrt(
                          as.numeric(t(Beta) %*% SIGMA %*% Beta)
                                                             ) %>% as.vector()

                        sigma2.eps.X = Var.X*(1 - h2.X)

                        U <- as.vector(SIGMA %*% Beta)

                        temp1 <- rnorm(P, 0, 1) %>% as.matrix()

                        temp2 <- rnorm(P, 0, 1) %>% as.matrix()

                        vpd1 <- crossprod(U, temp1)

                        vpd2 <- crossprod(U, temp2)

                        h1 = as.numeric(-vpd2/vpd1)

                        V1 <- temp1 * h1 + temp2

                        V1 = V1*sqrt(Var.X*h2.X)/as.numeric(sqrt(t(V1)%*%SIGMA%*%V1)) # normalize V1

                        #let Gamma = a1 * Beta + a2 * V1
                        # a1 ^2 = gc / h2.x
                        # a1^2 + a2^2 = h2.z/h2.x
                        #normalize gamma

                        a1 <- (sqrt(Var.X*Var.Z)*GeCovaraince/ (h2.X*Var.X) )
                        a2 <- sqrt((Var.Z*h2.Z)/(Var.X*h2.X) -  a1^2)
                        Gamma1 <<- as.vector( a1 * Beta + a2* V1)
                        Gamma1 <<- Gamma1*sqrt(Var.Z*h2.Z)/ as.numeric(sqrt(t(Gamma1)%*%SIGMA %*% Gamma1))

                        #return(list(Beta, Gamma1))
                      },
                      Coef.Prop = function(prop, SIGMA, h2.X, h2.Z, GeCovaraince, alpha = 0,  Var.X = 1, Var.Z = 1){

                        Beta.temp <- rep(0, P)

                        P.eff = floor(prop*P)

                        index <- sample(1:P, P.eff)

                        Beta.temp[index] = rnorm(P.eff)

                        Beta <<- Beta.temp

                        Beta <<- Beta*sqrt(h2.X*Var.X)/sqrt(
                          as.numeric(t(Beta) %*% SIGMA %*% Beta)
                        ) %>% as.vector()

                        U <- SIGMA %*% Beta

                        sigma2.eps.X = Var.X*(1 - h2.X)

                        U <- as.vector(SIGMA %*% Beta)

                        temp1 <- rep(0, P)

                        temp2 <- rep(0, P)

                        temp1[index] <- rnorm(P.eff)

                        temp2[index] <- rnorm(P.eff)

                        vpd1 <- crossprod(U, temp1)

                        vpd2 <- crossprod(U, temp2)

                        h1 = as.numeric(-vpd2/vpd1)

                        V1 <- temp1 * h1 + temp2

                        V1 = V1*sqrt(Var.X*h2.X)/as.numeric(sqrt(t(V1)%*%SIGMA%*%V1)) # normalize V1
                        #let Gamma = a1 * Beta + a2 * V1
                        # a1 ^2 = gc / h2.x
                        # a1^2 + a2^2 = h2.z/h2.x
                        #normalize gamma
                        a1 <- (sqrt(Var.X*Var.Z)*GeCovaraince/ (h2.X*Var.X) )
                        a2 <- sqrt((Var.Z*h2.Z)/(Var.X*h2.X) -  a1^2)
                        Gamma1 <<- as.vector( a1 * Beta + a2* V1)
                        Gamma1 <<- Gamma1*sqrt(Var.Z*h2.Z)/ as.numeric(sqrt(t(Gamma1)%*%SIGMA %*% Gamma1))

                      },
                      Coef.LDAK.Prop = function(prop, SIGMA, h2.X, h2.Z, GeCovaraince, alpha = 0.75, Var.X = 1, Var.Z = 1, MAF,LDW){

                        Beta.temp <-  rnorm(P, mean = rep(0, P), sd = (MAF *(1 - MAF))^(-alpha) * 1/LDW  )

                        P.eff = floor(prop*P)

                        index <- sample(1:P, P.eff)
                          if(P.eff != P){
                        Beta.temp[-index] = 0
                          }

                        Beta <<- Beta.temp

                        Beta <<- Beta*sqrt(h2.X*Var.X)/sqrt(
                          as.numeric(t(Beta) %*% SIGMA %*% Beta)
                        ) %>% as.vector()

                        U <- SIGMA %*% Beta

                        sigma2.eps.X = Var.X*(1 - h2.X)

                        U <- as.vector(SIGMA %*% Beta)

                        temp1 <- rep(0, P)

                        temp2 <- rep(0, P)

                        temp1[index] <- rnorm(P.eff)

                        temp2[index] <- rnorm(P.eff)

                        vpd1 <- crossprod(U, temp1)

                        vpd2 <- crossprod(U, temp2)

                        h1 = as.numeric(-vpd2/vpd1)

                        V1 <- temp1 * h1 + temp2

                        V1 = V1*sqrt(Var.X*h2.X)/as.numeric(sqrt(t(V1)%*%SIGMA%*%V1)) # normalize V1
                        #let Gamma = a1 * Beta + a2 * V1
                        # a1 ^2 = gc / h2.x
                        # a1^2 + a2^2 = h2.z/h2.x
                        #normalize gamma
                        a1 <- (sqrt(Var.X*Var.Z)*GeCovaraince/ (h2.X*Var.X) )
                        a2 <- sqrt((Var.Z*h2.Z)/(Var.X*h2.X) -  a1^2)
                        Gamma1 <<- as.vector( a1 * Beta + a2* V1)
                        Gamma1 <<- Gamma1*sqrt(Var.Z*h2.Z)/ as.numeric(sqrt(t(Gamma1)%*%SIGMA %*% Gamma1))

                      },
                      Orthogonal = function(x, SIGMA, GeCov, seed = 123){

                        h2.X = 1

                        Var.X = 1

                        h2.Z = 1

                        Var.Z = 1

                        U <- SIGMA %*% x

                        set.seed(seed)

                        temp1 <- rnorm(P, 0, 1)

                        temp2 <- rnorm(P, 0, 1)

                        vpd1 <- crossprod(U, temp1)

                        vpd2 <- crossprod(U, temp2)

                        h1 = -vpd2/vpd1

                        V1 <- h1 * temp1 + temp2

                        V1 = V1*sqrt(Var.X*h2.X)/sqrt(t(V1)%*%SIGMA%*%V1) # normalize V1

                        #let Gamma = a1 * Beta + a2 * V1
                        # a1 ^2 = gc / h2.x
                        # a1^2 + a2^2 = h2.z/h2.x
                        #normalize gamma

                        a1 <- (sqrt(Var.X*Var.Z)*GeCov/ (h2.X*Var.X) )
                        a2 <- sqrt((Var.Z*h2.Z)/(Var.X*h2.X) -  a1^2)
                        y <- a1 * x + a2* V1

                        return(y)

                        #return(list(Beta, Gamma1))
                      }

                    )
)

bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}
