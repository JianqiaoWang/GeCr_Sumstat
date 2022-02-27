
Coefs = function(Beta, SIGMA, h2.X, h2.Z, GeCov){

Beta = Beta*sqrt(h2.X*Var.X)/sqrt(t(Beta)%*%SIGMA%*%Beta)
sigma2.eps.X = Var.X*(1 - h2.X)

U <- SIGMA %*% Beta

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
Gamma1 <- a1 * Beta + a2* V1
Gamma1=Gamma1*sqrt(Var.Z*h2.Z)/sqrt(t(Gamma1)%*%SIGMA%*%Gamma1)

return()

}
