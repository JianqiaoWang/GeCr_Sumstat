#' @title Genetic correlation class
#' @description This class is only defined for the estimation of genetic correlation and heritability.
#' @param input_parameter <Input Parameter Description>
#' @export <Do you want users to be able to use this function? If not, it is a
#' @keywords
#' @seealso
#' @import <Load entire package as dependency>
#' @importFrom <Load specific functions from packages as dependencies>
#' @return Object returned
#' @aliases
#' @examples

library(R6)

##############################################
#
# Calculating the esimated variance of heritability
#
##############################################

norm_vec <- function(x) sqrt(sum(x^2))


h2.est.var = function(n = NA, y, S, h2.est = NA, h2.beta.A.2 = NA, Fnorm = NA, Tr = NA, sumsta = NA, exact = T){


  stopifnot(!is.na(n),
            !is.na(h2.est),
            !is.na(h2.beta.A.2),
            !is.na(Fnorm),
            !is.na(sumsta)
            )

  h2.beta.A = h2.est

  h2.beta.A.2 = h2.beta.A.2

 # if(exact){ # potential feature

  # if(is.na(h2.beta.A.2)){
  #
  #   h2.beta.A.2 = h2.est
  # }
  #
  # h2.beta.A = h2.est
  #
  # }else{
  #
  #   h2.beta.A.2 = h2.est
  #
  #   h2.beta.A = h2.est
  #
  # }

  if(sumsta != TRUE){

    res <- var(y)^2*2/n*
      (

        #  1/n * sum(diag((S %*% A %*% S %*% A))) +
        1/n * Fnorm^2 +

          2* h2.beta.A.2  +

          2* h2.beta.A^2

      )

  }else{

    res <- 2/n*
      (

        #  1/n * sum(diag((S %*% A %*% S %*% A))) +
        1/n * Fnorm^2 +

          2* h2.beta.A.2  -

          h2.beta.A^2

      )
  }

  return(res)

}


GC.est.var = function(n1, n2, y.s = NULL, y2.s = NULL, h2.beta.est, h2.gamma.est,
                      h2.beta.A.2 = NA, h2.gamma.A.2 = NA, rho.est, Tr , Fnorm, sumsta, exact){


    if(is.na(h2.beta.A.2) & is.na(h2.gamma.A.2)){

    h2.beta.A.2 = max(h2.beta.est,0)

    h2.gamma.A.2 = max(h2.gamma.est,0)
    }

if(sumsta == 1) #print("variance is biased")

  Delta.2.X = 1 # Potential feature

  Delta.2.Z = 1

  rho.A = rho.est

  res <- 1/(n1*n2)*(Delta.2.X)*(Delta.2.Z)*

    # sum(diag((S %*% A %*% S %*% A))) +
    #  sum(diag((Sigma %*% A %*% Sigma %*% A))) +
    Fnorm^2 +

    1/n2 *(Delta.2.X)*(Delta.2.Z) * h2.beta.A.2  +

    1/n1 *(Delta.2.X)*(Delta.2.Z) * h2.gamma.A.2  +

    (1/n1 + 1/n2) * (Delta.2.X)*(Delta.2.Z)*(rho.A)^2

  return(res)
}

GeCov <- R6Class("GeCr", list(
  Beta_mg = NULL,
  Gamma_mg = NULL,
  N1 = NA,
  N2 = NA,
  P = NA,
  y.s = NULL,
  y2.s = NULL,
  sumsta = T,
  simplified = F,
  h2.beta.est = NA,
  h2.gamma.est = NA,
  h2.beta.A.2 = NA,
  h2.gamma.A.2 = NA,
  rho.est = NA,
  rho.est.A.2 = NA,
  GeCr.est = NA,
  h2.beta.est.se = NA,
  h2.gamma.est.se = NA,
  rho.est.se = NA,
  GeCr.est.se = NA,
  cov.1.se = NA,
  cov.2.se = NA,
  Tr = NA,
  Fnorm = NA,
  Type = NA,

  initialize = function(Beta_mg,  Gamma_mg,  N1, N2, y.s = NA, y2.s = NA, sumsta.only = 1){

    stopifnot(is.vector(as.vector(Beta_mg)), is.vector(as.vector(Gamma_mg)), length(Beta_mg) == length(Gamma_mg))

    stopifnot(is.numeric(N1), is.numeric(N2))

    self$Beta_mg = Beta_mg

    self$Gamma_mg = Gamma_mg

    self$N1 = N1

    self$N2 = N2

    self$P = length(Beta_mg)

   # self$Fnorm = sqrt(self$P)

    self$sumsta = sumsta.only

    self$y.s = y.s

    self$y2.s = y2.s

  },

  print = function(...){
    cat(" A class for genetic correlation:\n")
    cat("The heritability  of beta_mg is: ",  self$h2.beta.est, "\n", "standard error is ", self$h2.beta.est.se,"\n", sep = "")
    cat("The heritability  of gamma_mg is: ",  self$h2.gamma.est, "\n", "standard error is ", self$h2.gamma.est.se,"\n", sep = "")
    cat("The Genetic Correlation is: ",  self$GeCr.est, "\n", "standard error is ", self$GeCr.est.se, sep = "")
  }

))

GeCov$set( "public", "Estimation", function(Sigma = NULL,
                                            Omega = NULL,
                                            LDSCORE = NULL,
                                            Type = NULL){

  library(Matrix)

  stopifnot(!is.null(Omega))

  stopifnot(!is.null(Sigma))

  ###################################################################################
  #
  # Calculate the moments of LD matrix,  Trace( A Sigma) and Trace( A Sigma  A Sigma)
  #
  ##################################################################################

  if(Type == "proposed"){

    self$Tr = as.numeric(sum(Omega * Sigma))

    self$Fnorm = sqrt(as.numeric(self$P))

  }else if(Type == "Dicker"){

      self$Tr = as.numeric(self$P)

      self$Fnorm = sqrt(as.numeric(self$P))

  }else if(Type == "Identity"){

    # calculate the Trace ( Sigma ) and Trace ( Sigma^2 )

      self$Tr = as.numeric(sum(diag(Sigma)))

      self$Fnorm = ifelse((self$sumsta != 1 | is.null( LDSCORE )) ,
             sqrt(sum( Sigma^2 )) ,
             sqrt(as.numeric(sum( LDSCORE)) ) )

    }

#####################################################################
#
# heritability, genetic covariance,  and genetic correlation estimates
#
######################################################################

    self$Type = Type


    self$rho.est <- as.numeric( t(self$Beta_mg) %*% Omega %*% self$Gamma_mg)

  if(self$sumsta == 1){

    Q.beta =  as.numeric(

      (t(self$Beta_mg) %*% Omega %*% (self$Beta_mg) - self$Tr/ self$N1)
    )


    Q.gamma <- as.numeric(
      (t(self$Gamma_mg) %*% Omega %*% (self$Gamma_mg) - self$Tr/self$N2)
    )

  }else{

    Q.beta <-  as.numeric(
      (t(self$Beta_mg) %*% Omega %*% (self$Beta_mg) - self$Tr/self$N1^2 * norm_vec(self$y.s)^2)
    )

    Q.gamma <- as.numeric(
      (t(self$Gamma_mg) %*% Omega %*% (self$Gamma_mg) - self$Tr/self$N2^2* norm_vec(self$y2.s)^2)
    )
  }

  #Q.beta =   Q.beta * (Q.beta > 0 & Q.beta < 1) + 1 * (Q.beta >= 1) + 0*(Q.beta <= 0)

  #Q.gamma =  Q.gamma * (Q.gamma > 0 & Q.gamma < 1) + 1 * (Q.gamma >= 1) + 0*(Q.gamma <= 0)

  self$h2.beta.est = Q.beta

  self$h2.gamma.est = Q.gamma

  try({
  self$GeCr.est = self$rho.est/sqrt(self$h2.beta.est * self$h2.gamma.est)
  })


#########################################################
#
# Higher order estimates (to calculate the variance)
#
############################################################

  if(Type == "Identity"){

  self$h2.beta.A.2 = as.numeric(

    (t(self$Beta_mg) %*% Sigma %*% (self$Beta_mg) - self$Fnorm^2/ self$N1)
  )

  self$h2.gamma.A.2 = as.numeric(

    (t(self$Gamma_mg) %*% Sigma %*% (self$Gamma_mg) - self$Fnorm^2/ self$N2)
  )

  self$rho.est.A.2 <- as.numeric( t(self$Beta_mg) %*% Sigma %*% self$Gamma_mg)

  }else{


    self$h2.beta.A.2 = self$h2.beta.est

    self$h2.gamma.A.2 = self$h2.gamma.est

    self$rho.est.A.2  = self$rho.est

  }




}

)

GeCov$set( "public", "Var.Est", function(exact){


  self$h2.beta.est.se = sqrt(
    as.numeric(h2.est.var(n = self$N1, y = self$y.s, h2.est = self$h2.beta.est,
                          h2.beta.A.2 =  self$h2.beta.A.2,
                          Tr = self$Tr , Fnorm = self$Fnorm, sumsta = self$sumsta, exact = exact))
  )

  self$h2.gamma.est.se = sqrt(
    as.numeric(
      h2.est.var(n = self$N2, y = self$y2.s,
                 h2.est = self$h2.gamma.est,  h2.beta.A.2 =  self$h2.gamma.A.2,
                 Tr = self$Tr , Fnorm = self$Fnorm, sumsta = self$sumsta, exact = exact))
  )

  self$rho.est.se = sqrt(
    as.numeric(
      GC.est.var(n1 = self$N1, n2 = self$N2, h2.beta.est = self$h2.beta.est,
                 h2.gamma.est = self$h2.gamma.est, h2.beta.A.2 = self$h2.beta.A.2,
                 h2.gamma.A.2 = self$h2.gamma.A.2,
                 rho.est = self$rho.est,
                Tr = self$Tr , Fnorm = self$Fnorm, sumsta = self$sumsta, exact = exact))
  )

  n1 = self$N1

  n2 = self$N2

  phi.1 = 1/(n1*n2)*

    self$Fnorm^2 +

    1/n1 * self$h2.beta.est  +

    1/n2 * self$h2.gamma.est  +

    (1/n1 + 1/n2) * (self$rho.est)^2


  phi.6 = 1/(n1)*( self$Fnorm^2/n1 -

                     2*self$h2.beta.est -

                     2*self$h2.beta.est^2 )


  phi.7 =  1/(n2)*(self$Fnorm^2/n2 -

                      2*self$h2.gamma.est -

                      2*self$h2.gamma.est^2)

  Heritability.1 <- self$h2.beta.est

  Heritability.2 <- self$h2.gamma.est

  Genetic.Covariance <- self$rho.est

  self$cov.1.se <- 2/n1* (self$rho.est + Genetic.Covariance * Heritability.1)

  self$cov.2.se <- 2/n2* (self$rho.est + Genetic.Covariance * Heritability.2)

  self$GeCr.est.se =  sqrt(

    1/(Heritability.1 * Heritability.2) * phi.1 +
      1/2 * Genetic.Covariance^2/(Heritability.1^3 * Heritability.2)* phi.6 +
      1/2 * Genetic.Covariance^2/(Heritability.1 * Heritability.2^3)*  phi.7

  )

  if(self$Type == "Identity"){

    phi.1 = 1/(n1*n2)*

      self$Fnorm^2 +

      1/n1 * self$h2.beta.A.2  +

      1/n2 * self$h2.gamma.A.2  +

      (1/n1 + 1/n2) * (self$rho.est)^2


    phi.6 = 2/(n1)*( self$Fnorm^2/n1 +

                       2*self$h2.beta.A.2 +

                       2*self$h2.beta.est^2 )


    phi.7 =  2/(n2)*(self$Fnorm^2/n2 +

                       2*self$h2.gamma.A.2 +

                       2*self$h2.gamma.est^2)


    self$cov.1.se <- 2/n1* (self$rho.est.A.2 + Genetic.Covariance * Heritability.1)

    self$cov.2.se <- 2/n2* (self$rho.est.A.2 + Genetic.Covariance * Heritability.2)

    self$GeCr.est.se  =  sqrt(

      1/(Heritability.1 * Heritability.2) * phi.1 +
        1/4 * Genetic.Covariance^2/(Heritability.1^3 * Heritability.2)* phi.6 +
        1/4 * Genetic.Covariance^2/(Heritability.1 * Heritability.2^3)* phi.7 -
       Genetic.Covariance/(Heritability.1^2 * Heritability.2)* self$cov.1.se -
          Genetic.Covariance/(Heritability.1 * Heritability.2^2)* self$cov.2.se

    )

  }



}
)


GeCov$set( "public", "Var.Est.limited", function(N){


  self$h2.beta.est.se = sqrt(
    2*(self$h2.beta.est)^2/N
  )

  self$h2.gamma.est.se = sqrt(
    2*(self$h2.gamma.est)^2/N
  )

  self$rho.est.se = sqrt(
    as.numeric(
      (self$h2.gamma.est * self$h2.beta.est + self$self$rho.est^2)/N
  )
)

  Heritability.1 <- self$h2.beta.est

  Heritability.2 <- self$h2.gamma.est

  Genetic.Covariance <- self$rho.est

  cov.gamma.I = 2/N * self$rho.est* self$h2.gamma.est

  cov.beta.I = 2/N * self$rho.est* self$h2.beta.est

  cov.beta.gamma = 2/N * self$rho.est^2

  # self$GeCr.est.se =  sqrt(
  #
  #   1/(Heritability.1 * Heritability.2) *self$rho.est^2 +
  #     1/4 * Genetic.Covariance^2/(Heritability.1^3 * Heritability.2)* self$h2.beta.est.se^2 +
  #     1/4 * Genetic.Covariance^2/(Heritability.1 * Heritability.2^3)* self$h2.gamma.est.se^2 +
  #     2*(-1/2)*Genetic.Covariance/(Heritability.1^2 * Heritability.2)*cov.beta.I +
  #     2*(-1/2)*Genetic.Covariance/(Heritability.1 * Heritability.2^2)*cov.gamma.I +
  #     2*(1/4)*(Genetic.Covariance^2/Heritability.1^2 * Heritability.2^2)*cov.beta.gamma
  # )

  self$GeCr.est.se =  sqrt(

    1/N -2/N* self$GeCr.est^2 + self$GeCr.est^4/N
  )



}
)



GeCov$set( "public", "PVal", function(){

  Z.beta = self$h2.beta.est/self$h2.beta.est.se

  p.beta = 2*pnorm(abs(Z.beta), lower.tail = F)

  Z.gamma = self$h2.gamma.est/self$h2.gamma.est.se

  p.gamma = 2*pnorm(abs(Z.gamma),lower.tail = F)

  Z.rho = self$rho.est/self$rho.est.se

  p.rho = 2*pnorm(abs(Z.rho),lower.tail = F)

  Z.GeCr = self$GeCr.est/self$GeCr.est.se

  p.GeCr = 2*pnorm(abs(Z.GeCr), lower.tail = F)

  return(c(p.beta, p.gamma, p.rho, p.GeCr))

})

Theoretical.GECR.SE = function(N.1, N.2, P, h2.1, h2.2, R){

 theoretical.se =  sqrt(

   R^2/ (2* h2.1^2 ) * (P/N.1^2) +  R^2/ (2* h2.2^2 ) * (P/N.2^2) +  (1/(h2.1*h2.2))* (P/(N.2 * N.1 )) +
     (1/N.1) * (1 - R^2)/h2.1 +   (1/N.2) * (1 - R^2)/h2.2
  )

 return(theoretical.se)

}

