# simulation for normal covariates
rm(list = ls())
#library(foreach)
#library(doParallel)
library(foreach)
library(doParallel)
library(mvtnorm)
library(MASS)
library(psych)
library(copula)
library(parallel)
library(stringr)


##############################################
# function sections
##############################################
source('../R/LDScore.R')
source('../R/CoefGen.R')
source('../R/Genetic_Correlation_class_v2.R')
source('../R/GenX.R')
source('../R/SummaryStatistics.R')
source('../R/Precsion_Matrix_Estimation.R')

##############################################
# Information needed for simulating dataset
##############################################
nblock = 100
n = 1000
m = n
n1 = n # number of subject of study 1
n2 = m # number of subject of study 2
Pblock = rep(100, nblock) #number of SNP 100 within each block
P = sum(Pblock)
rho.vector = runif(nblock, min = 0, max = 0.9)
SIGMA = lapply(1:nblock, function(i) {
  rho.vector[i] ^ abs(outer(1:Pblock[i], 1:Pblock[i], "-"))
}) %>% bdiag_m

Omega0 = lapply(1:nblock, function(i) {
solve(rho.vector[i] ^ abs(outer(1:Pblock[i], 1:Pblock[i], "-")))
}) %>% bdiag_m

nsim = 300 #number of dataset simulated, if nsim=1, only 1 dataset is generated
GeCovaraince = 0.2 # genectic covariance
h2.X = 0.4 # heritability of trait 1
h2.Z = 0.4 # heritability of trait 2
R = GeCovaraince / (sqrt(h2.X) * sqrt(h2.Z)) # genetic correlation
alpha = 0.05
Var.Z = 1
Var.X = 1
type2 = F
sigma2.eps.X = Var.X * (1 - h2.X) #  variance of error term 1
sigma2.eps.Z = Var.Z * (1 - h2.Z) #  variance of error term 2

##################### generation of LD score  ##################################
if(F){
LDSCORE <- Matrix::rowSums(SIGMA^2)
ldscore.data = data.frame(cbind(1:length(SIGMA), 1,
                                paste0("rs",1:length(SIGMA)), LDSCORE))
colnames(ldscore.data) = c("BP", "CHR", "SNP", "LD")
fwrite(ldscore.data, file = "simu.l2.ldscore", quote = F, sep = "\t",
       col.names=T, na = "NA" )
fwrite(data.frame(P), file = "simu.l2.M_5_50", quote = F, sep = "\t",
       col.names=F, na = "NA" )
weight.ld.data = data.frame(cbind(1:length(SIGMA), 1,
                                  paste0("rs",1:length(SIGMA)), 1))
colnames(weight.ld.data) = c("BP", "CHR", "SNP", "LD")
fwrite(weight.ld.data, file = "weight.l2.ldscore", quote = F, sep = "\t",
       col.names=T, na = "NA" )
fwrite(data.frame(P), file = "weight.l2.M_5_50", quote = F, sep = "\t",
       col.names=F, na = "NA" )
}

######################## Single simulation ###################################
######### step.1 generation of (X, Y1) ans (Z, Y2)
coef.beta.gamma = Coef$new(p =P)
coef.beta.gamma$Coef.Prop(prop = prop, SIGMA, h2.X, h2.Z, GeCovaraince)
Beta = coef.beta.gamma$Beta
Gamma1 = coef.beta.gamma$Gamma1



Simulation = function(i, N) {

  set.seed(i)

  ######################
  #Generate data;
  ######################

  X = Design$block.normal(n, Pblock, rho.vector, mc.cores = 15)

  X.eps <- rnorm(n, 0, sd = sqrt(sigma2.eps.X))
  Y1 = X %*% Beta + X.eps
  beta_mg = as.vector(crossprod(X, Y1)/n)

  rm(X)

  Z = Design$block.normal(m, Pblock, rho.vector, mc.cores = 15)
  Z.eps <- rnorm(m, 0, sd = sqrt(sigma2.eps.Z))
  Y2 = Z %*% Gamma1 + Z.eps
  gamma_mg = as.vector(crossprod(Z, Y2)/m)


  rm(Z)

  gc()

  beta_t = beta_mg * sqrt(n)

  gamma_t = gamma_mg * sqrt(m)

  ################  Split external data into two parts and estimate covariance

  I = Matrix::Diagonal(P)

  GeCor.Analysis = GeCov$new(Beta_mg = beta_mg,
                             Gamma_mg = gamma_mg,
                             N1 = n,
                             N2 = m,
                             sumsta = F,
                             y.s = Y1,
                             y2.s = Y2
  )

  GeCor.Analysis$Estimation(Omega = Omega0, Sigma = SIGMA, Type = "proposed")

  GeCor.Analysis$Var.Est(TRUE)
  #
  Proposed = c(GeCor.Analysis$GeCr.est,
               GeCor.Analysis$GeCr.est.se,
               GeCor.Analysis$h2.beta.est,
               GeCor.Analysis$h2.beta.est.se,
               GeCor.Analysis$h2.gamma.est,
               GeCor.Analysis$h2.gamma.est.se,
               GeCor.Analysis$rho.est,
               GeCor.Analysis$rho.est.se,
               GeCor.Analysis$cov.1.se,
               GeCor.Analysis$cov.2.se,
               Theoretical.GECR.SE(N.1 = n1, N.2 = n2,
                                   P = P, h2.1 = h2.X,
                                   h2.2 = h2.Z, R = R))


  print(Proposed)

  GeCor.Analysis$Estimation(Omega = I, Sigma = SIGMA, Type = "Identity")

  GeCor.Analysis$Var.Est(TRUE)

  Identity = c(GeCor.Analysis$GeCr.est,
               GeCor.Analysis$GeCr.est.se,
               GeCor.Analysis$h2.beta.est,
               GeCor.Analysis$h2.beta.est.se,
               GeCor.Analysis$h2.gamma.est,
               GeCor.Analysis$h2.gamma.est.se,
               GeCor.Analysis$rho.est,
               GeCor.Analysis$rho.est.se)

  LDSC.CSTR =GR.LDSC.EST.SE(J = i, Z1 = beta_t, Z2 = gamma_t, N1 = n1, N2 = n2,
                            LDSCORE =  LDSCORE, CSTR = T,  weight=T)

  LDSC.NO.Weight =GR.LDSC.EST.SE(J = i, Z1 = beta_t, Z2 = gamma_t, N1 = n1, N2
                                 = n2, LDSCORE =  LDSCORE, CSTR = T,  weight=F)

  # LDSC.NO.CSTR =GR.LDSC.EST.SE(Z1 = beta_t, Z2 = gamma_t, #
  #N1 = n1, N2 = n2, LDSCORE =  LDSCORE, CSTR = F)

  print(Identity)

  LDSC.NO.CSTR = NA

  return(list(Proposed = Proposed, Identity = Identity, LDSC.CSTR = LDSC.CSTR, LDSC.NO.Weight = LDSC.NO.Weight))


}

list.all = list()



for(alternative1 in c(1)){


  for(rho in c(0.8)){
    for(GeCovaraince in c(0.2 )){

      X.norm = TRUE

      R = GeCovaraince/(sqrt(h2.X)*sqrt(h2.Z))

      z_alpha <- qnorm(1 - alpha/2)

      ##############################################
      #Create variables to hold the results;

      df = data.frame()

      Design = Regressor$new(P = P , rho = rho)

      coef.beta.gamma = Coef$new(p = P)

      ### Generate beta and gamma

      if(alternative1 == 1){

        coef.sigma = matrix(c(h2.X, GeCovaraince, GeCovaraince, h2.Z),2,2)/P

        coef.matrix = rmvnorm(n=P, mean= c(0,0), sigma= coef.sigma)

        Beta = coef.matrix[, 1]

        Gamma1 = coef.matrix[, 2]

      }

      if(alternative1  == 2){

        prop = 0.0001

        set.seed(123)

        coef.beta.gamma$Coef.Prop(prop = prop, SIGMA, h2.X, h2.Z, GeCovaraince)

        Beta = coef.beta.gamma$Beta

        Gamma1 = coef.beta.gamma$Gamma1

      }

      if(alternative1 == 3){


        edc = eigen(SIGMA, symmetric = T)

        index.1 = which(edc$values < 1 )[1]

        index.2 = 600

        index.3 = 800



        lambda.1 = edc$values[index.1]

        lambda.2 = edc$values[index.2]

        lambda.3 = edc$values[index.3]


        U.1 = edc$vectors[,index.1]

        U.2 = edc$vectors[,index.2]

        U.3 = edc$vectors[,index.3]


        fn <- function(a.1, a.2, a.3, b.1, b.2, b.3){

          beta.sigma.2.gamma <- a.1  * b.1 * lambda.1^2 + a.2  * b.2 * lambda.2^2 + a.3  * b.3 * lambda.3^2

          beta.sigma.gamma <- a.1  * b.1 * lambda.1 + a.2  * b.2 * lambda.2 +  a.3  * b.3 * lambda.3

          beta.sigma.beta <- a.1^2 * lambda.1 + a.2^2 * lambda.2 +  a.3^2 * lambda.3

          gamma.sigma.gamma <- b.1^2 * lambda.1 + b.2^2 * lambda.2 + b.3^2 * lambda.3

          return(c(gamma.sigma.gamma, beta.sigma.beta, beta.sigma.gamma, beta.sigma.2.gamma))

        }

        fn2 <- function(x) sum( (fn(x[1], x[2], x[3], x[4], x[5], x[6]) - c(h2.Z, h2.X, GeCovaraince,0))^2)

        fn3 <- function(x)  (fn(x[1], x[2], x[3], x[4], x[5], x[6]) - c(h2.Z, h2.X, GeCovaraince,0))

        opt.result =  nlm(fn2, c(0.5, 0.5, 0.5, 0.5,0.5,0.5))

        #  opt.result =  optim(c(0.5, 0.5, 0.5, 0.5,0.5,0.5), fn2)

        a.1 = opt.result$estimate[1]

        a.2 = opt.result$estimate[2]

        a.3 = opt.result$estimate[3]


        b.1 = opt.result$estimate[4]

        b.2 = opt.result$estimate[5]

        b.3 = opt.result$estimate[6]


        ### Generate orthogonal vectors U.1 and U.2 to the vector Gamma.

        Beta =   a.1 * U.1 + a.2 * U.2 + a.3 * U.3

        Gamma1 =  b.1 * U.1 + b.2 * U.2 + b.3 * U.3

        t(Beta) %*% SIGMA %*% SIGMA %*% Gamma1/ sqrt( t(Beta) %*% SIGMA %*% SIGMA %*% Beta * t(Gamma1) %*% SIGMA %*% SIGMA %*% Gamma1 )

        t(Beta) %*% SIGMA  %*% Gamma1/ sqrt( t(Beta) %*% SIGMA  %*% Beta * t(Gamma1) %*% SIGMA %*% Gamma1)

      }


      Flag.scale <- 1

      n1 = n

      n2 = m

      set.seed(123)

      temp = mclapply(1:nsim, Simulation, N = N, mc.cores = 3)

      name = paste0("Xnorm",X.norm,"rho",rho, "GeCv", GeCovaraince, "alternative1", alternative1)

      list.all[[name]] =  temp

    }

  }
}

save(list.all, file = paste0("Experiment2_v2.Rdata"))

####### Note :
##


#temp = (eval(parse(text = paste("result", 0, sep = ""))))

#H <- data.frame(matrix(t(temp), ncol = 2, byrow = T,
#                       dimnames = list(c(),c("est.mean","est.sd"))
#))

# H$name <- rep(c("GeCov","h2.X","h2.Z","R"), times = nsim)
#
# if(Flag.scale == 1){
#
# H$tp <- rep(c(GeCov,h2.X,h2.Z,R), times = nsim)
#
# }else{
#
#   H$tp <- rep(c( t(Beta)%*% SIGMA %*% Gamma1,Var.X * h2.X, Var.Z * h2.Z,R), times = nsim)
#
# }

# print(c(mean(H$est.mean[H$name == "h2.X"]),
#
#       mean(H$est.mean[H$name == "h2.Z"]),
#
#       mean(H$est.mean[H$name == "GeCov"]),
#
#       mean(H$est.mean[H$name == "R"],na.rm = T),
#
#       sd(H$est.mean[H$name == "h2.X"]),
#
#       sd(H$est.mean[H$name == "h2.Z"]),
#
#       sd(H$est.mean[H$name == "GeCov"]),
#
#       sd(H$est.mean[H$name == "R"],na.rm = T)))
#
# }
#
# print(H)
