
#setwd("..")

alpha = 0.05

#R = GeCovaraince/(sqrt(h2.X)*sqrt(h2.Z))

R = 0.5

nsim = 200

rho = 0.8

X.norm = T

method = "clime"

#method = "invthreshold"

estimator = "Proposed" # or Dicker

#estimator = "Dicker"

extract_result = function(R, alpha, method, estimator){

z_alpha <- qnorm(1 - alpha/2)

#result.data = data.frame(c(NA,NA,NA,NA),col.names = c("N.200", "N.500", "N.1000", "N.4000"))

result.data = vector()

for(N in c(200, 500, 1000, 4000)){

 # rm(temp)

  load(file = paste0("clime_N", N,"Xnorm",X.norm,"rho",rho,"0729",".Rdata"))

  proposed = do.call(rbind, purrr::modify_depth(temp,1,estimator))

  proposed = as.data.frame(proposed)

  #print(proposed)

  proposed$lower <- proposed[,4] - proposed[,8] * z_alpha

  proposed$upper <- proposed[,4] + proposed[,8] * z_alpha

  count <- 0

  count = sum(((R >= proposed$lower) + (R <= proposed$upper) == 2),na.rm = T)

  result.data = rbind(result.data, c((count), count/nsim))

}

result.data = data.frame(result.data)

result.data$N = c("N.200", "N.500", "N.1000", "N.4000")

return(result.data)
}

vector1 = data.frame()


for(X.norm in c(TRUE, FALSE)){
  for(rho in c(0.4, 0.8)){
    for(R in c(0.5)){
      for(method in c("clime")){
        for(estimator in c("Proposed", "Dicker")){

          temp1 = extract_result(R, alpha, method, estimator)

          temp1$rho = rho

          temp1$R = R

          temp1$estimator = estimator

          temp1$X.norm = X.norm

          vector1 = rbind(vector1, temp1)

        }
      }
    }
  }
}

vector_temp = vector1[,-1]

simulation_result = tidyr::spread(vector_temp, N, value = c(X2))

name = c("rho","R", "X.norm", "estimator", "N.200",  "N.500", "N.1000" ,"N.4000")

simulation_result = simulation_result %>% select(name)
library(knitr)
kable(simulation_result, format = "latex")
