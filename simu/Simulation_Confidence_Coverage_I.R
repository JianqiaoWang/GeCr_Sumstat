rm(list = ls())
library(foreach)
library(doParallel)
library(mvtnorm)
library(MASS)
library(flare)
library(MASS)
library(psych)
library(copula)
library(flare)

#setup parallel backend to use many processors
#library(pracma)
#library(glasso)

##############################################
# function sections
##############################################


#source("../PackageFunction/mat-choice-old.R")
#source("../PackageFunction/Thresholding-Sample-Covariance.R")
#source("../PackageFunction/FlareClime.R")
#source("../PackageFunction/LDScore.R")

##############################################
# Information needed for simulating dataset
##############################################
#args=(commandArgs(TRUE))

#X.norm = as.numeric(args[1])

#beta.abs = as.numeric(args[2])

#MatrixChoice = (args[3])


source('./R/CoefGen.R')
source('./R/Genetic_Correlation_class.R')
source('./R/GenX.R')
source('./R/SummaryStatistics.R')
source('./R/Precsion_Matrix_Estimation.R')


GeCovaraince = 0 #denominator
h2.X = 0.6
h2.Z = 0.5
R = GeCovaraince/(sqrt(h2.X)*sqrt(h2.Z))

P = 100 #number of SNP
n = 20000  #number of subject
m = n
nsim = 100 #number of dataset simulated, if nsim=1, only 1 dataset is generated
alpha = 0.05
Var.Z = 1
Var.X = 1

######
type2 = TRUE



z_alpha <- qnorm(1 - alpha/2)

##############################################
#Create variables to hold the results;
result <- vector()
result1 <- vector()
result2 <- vector()
result3 <- vector()
result4 <- vector()
result5 <- vector()
result6 <- vector()
result7 <- vector()
result8 <- vector()
result0 <- vector()



MaxVec <- vector()
MinVec <- vector()
Mean1 <- vector()
Mean2 <- vector()
SD <- vector()
SD2 <- vector()
TrueValue <- vector()

df = data.frame()

#for(rho in c(0, 0.4, 0.6, 0.8)){

for(rho in c(0.6)){

  SIGMA =  rho^abs(outer(1:P, 1:P, "-"))

  Design = Regressor$new(P = P , Cov = SIGMA)

  coef.beta.gamma = Coef$new(p = P)

  Beta = coef.beta.gamma$sparse(5)

  coef.beta.gamma$Coefs(SIGMA, h2.X, h2.Z, GeCovaraince)

  Beta = coef.beta.gamma$Beta

  Gamma1 = coef.beta.gamma$Gamma1

  sigma2.eps.X = Var.X*(1 - h2.X)

  sigma2.eps.Z = Var.Z*(1 - h2.Z)

}


########
# source("./PackageFunction/GeCr_old.R")

# source("../PackageFunction/GeCr_simplified.R")

Flag.scale <- 1

n1 = n

n2 = n

#result0 <- matrix(NA, nrow = nsim, ncol = 8)

#cores=detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
#registerDoParallel(cl)
#message('Number of cores detected ',cores[1]-1,'system core information', Sys.getenv('LSB_DJOB_NUMPROC'))


Simulation = function(i, N){

      set.seed(i)

      ######################
      #Generate data;
      ######################

      X = Design$normal(n)

      Z = Design$normal(m)

      X.eps <- rnorm(n,0,sd = sqrt(sigma2.eps.X))

      Z.eps <- rnorm(n,0,sd = sqrt(sigma2.eps.Z))

      Y1 = X %*% Beta + X.eps

      Y2 = Z %*% Gamma1 + Z.eps

      beta_t = SumStat(Y1, X)

      gamma_t = SumStat(Y2, Z)

      beta_mg <- beta_t /sqrt(beta_t^2 + n - 2)

      gamma_mg <- gamma_t /sqrt(gamma_t^2 + n - 2)

      ########### test ##############
      beta_mg <- as.vector(SIGMA %*% Beta)

      gamma_mg <- as.vector(SIGMA %*% Gamma1)

      ###########################

      ################## Simulate External data ############

      ExternalData <- Design$normal(N)

      ExternalData <- scale(ExternalData)

      N.1 <- ceiling(N/2)

      index <- sample(1:N, size = N.1, replace = FALSE)

      ExternalData.Omega <- ExternalData[index ,]

      ExternalData.Sigma <- ExternalData[-index ,]

      ################  Split external data into two parts and estimate covariance

      #Omega0 = solve(SIGMA)

      #Sigma0 = SIGMA

      # Omega0 = FASTCLIME(X = ExternalData)

      if(type2){

        Omega.temp = Flare(X = ExternalData.Omega)

        Sigma0 = cor(ExternalData.Sigma)

        Omega0 = Omega.temp + Omega.temp%*%(diag(P) - Sigma0 %*% Omega.temp)

        # diag(diag(cor(ExternalData.Omega) %*% Omega.temp)^{-1})%*%

      }else{

      Omega0 = Flare(X = ExternalData)

       Sigma0 = cor(ExternalData)
      }

      ############################# Test for the codes ##########################

      GeCor.Analysis = GeCov$new(Beta_mg = beta_mg,
                                 Gamma_mg = gamma_mg,
                                 N1 = n,
                                 N2 = m)

      GeCor.Analysis$Estimation(Omega = Omega0, Sigma = Sigma0, Type = "proposed")

      if(type2){

        GeCor.Analysis$Var.Est.limited(N/2)

      }else{

      GeCor.Analysis$Var.Est(TRUE)

      }

      Proposed = c(GeCor.Analysis$h2.beta.est,
                   GeCor.Analysis$h2.gamma.est,
                   GeCor.Analysis$rho.est,
                   GeCor.Analysis$GeCr.est,
                   GeCor.Analysis$h2.beta.est.se,
                   GeCor.Analysis$h2.gamma.est.se,
                   GeCor.Analysis$rho.est.se,
                   GeCor.Analysis$GeCr.est.se)


      GeCor.Analysis$Estimation(Omega = Omega0, Sigma = Sigma0, Type = "Dicker")

      if(type2){

        GeCor.Analysis$Var.Est.limited(N/2)

      }else{

        GeCor.Analysis$Var.Est(TRUE)

      }

      Dicker = c(GeCor.Analysis$rho.est,
                 GeCor.Analysis$rho.est.se)





      print(list(STAT = c(GeCor.Analysis$h2.beta.est,
                           GeCor.Analysis$h2.gamma.est,
                           GeCor.Analysis$rho.est,
                           GeCor.Analysis$GeCr.est,
                           GeCor.Analysis$h2.beta.est.se,
                           GeCor.Analysis$h2.gamma.est.se,
                           GeCor.Analysis$rho.est.se,
                           GeCor.Analysis$GeCr.est.se),
                  Pval = GeCor.Analysis$PVal()
      ))

      return(list(Proposed = Proposed, Dicker = Dicker
      ))

}

library(parallel)

for(N in c(200, 500, 1000, 1200)){

temp = mclapply(1:nsim, Simulation, N = N, mc.cores = 1)

save(temp, file = paste0("clime_N", N,".Rdata"))

}


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
