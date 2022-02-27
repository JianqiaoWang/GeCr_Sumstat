library(data.table)
library(dplyr)

SHvr <- function(Z, N, LDSCORE,CSTR ){

  chi2 <-  Z^2

  if(CSTR == T){

    predictor = ((N/length(Z))) *LDSCORE

  lm.fit <- lm( (chi2 -1)  ~ predictor  + 0)

  h2 <- lm.fit$coefficients
  }else{

    predictor = ((N/length(Z))) *LDSCORE

    lm.fit <- lm( (chi2-1)  ~ predictor )

    print(lm.fit)

    h2 <- lm.fit$coefficients[2]

  }

  return(list(h2 = h2))
}


LDSC.EST <- function(J, Z1, N1, Z2, N2, LDSCORE, CSTR , Nc=0, W=NULL){




  Heritability.1 = SHvr(Z1, N1, LDSCORE, CSTR = CSTR)$h2

  Heritability.2 = SHvr(Z2, N2, LDSCORE, CSTR = CSTR)$h2

  Y = Z1*Z2

  if(CSTR == T){

  lm.fit <- lm(Y ~ LDSCORE  + 0)

  Genetic.Covariance <- lm.fit$coefficients[1]/(sqrt(N1 * N2)/length(Z1))

  }else{

    lm.fit <- lm(Y ~ LDSCORE)

    Genetic.Covariance <- lm.fit$coefficients[2]/(sqrt(N1 * N2)/length(Z1))

  }


  Genetic.Correlation = Genetic.Covariance/sqrt(Heritability.1*Heritability.2)

  return(Genetic.Correlation)

}


GR.LDSC.EST.SE <- function(J, Z1, N1, Z2, N2, LDSCORE,
                           Nc=0, weight=T, CSTR = T,
                           info = NA,
                           corenum = 1){


  if(length(Z1) != length(Z2) ){

    stop("Z scores should have the same length.")

  }

  if(!is.na(info)){

    print("update info")

    Z1.data = data.frame(SNP = info$SNP , Z = Z1, N =  N1,
                          A1 =  info$A1, A2 = info$A2, stringsAsFactors = F)

    colnames(Z1.data) = c("SNP",  "Z", "N", "A1", "A2")

    fwrite(Z1.data, file = paste0(J, "Z1.sumstats") , quote = F, sep = "\t", col.names=T, na = "NA" )

    Z2.data = data.frame(SNP = info$SNP , Z = Z2, N =  N2,
                          A1 = info$A1,  A2 =  info$A2, stringsAsFactors = F)

    colnames(Z2.data) = c("SNP", "Z", "N", "A1", "A2")

    fwrite(Z2.data, file = paste0(J, "Z2.sumstats"), quote = F, sep = "\t", col.names=T, na = "NA" )


  }else{

  Z1.data = as.data.frame(cbind(paste0("rs",1:length(Z1)), Z1, N1,  A1 = "A", A2 = "G"))

  colnames(Z1.data) = c("SNP",  "Z", "N", "A1", "A2")

  fwrite(Z1.data, file = paste0(J, "Z1.sumstats") , quote = F, sep = "\t", col.names=T, na = "NA" )

  Z2.data = as.data.frame(cbind(paste0("rs",1:length(Z2)), Z2, N2,  A1 = "A", A2 = "G"))

  colnames(Z2.data) = c("SNP", "Z", "N", "A1", "A2")

  fwrite(Z2.data, file = paste0(J, "Z2.sumstats"), quote = F, sep = "\t", col.names=T, na = "NA" )

  }

  if(weight == T){

  system(paste0("ldsc.py --rg ",paste0(J, "Z1.sumstats,"), paste0(J, "Z2.sumstats"),
                " --no-intercept --ref-ld simu --w-ld simu --out ",paste0(J, "simu")))
  }else{

    system(paste0("ldsc.py --rg ",paste0(J, "Z1.sumstats,"), paste0(J, "Z2.sumstats"),
                  " --no-intercept --ref-ld simu --w-ld weight --out ",paste0(J, "simu")))
  }

  H = fread(paste0(J, "simu.log"), fill = T)

  result.all = rep(NA, 8)

  result.gecr = unlist(H[(stringr::str_detect(unlist(H), "Genetic Correlation:"))]) %>%
    str_split(, pattern = " ") %>% unlist()

  result.heri = unlist(H[(stringr::str_detect(unlist(H), "Total Observed scale h2:"))]) %>%
    str_split(, pattern = " ") %>% unlist()

  result.gecov = unlist(H[(stringr::str_detect(unlist(H), "Total Observed scale gencov:"))]) %>%
    str_split(, pattern = " ") %>% unlist()

  if(!is.null(result.gecr)){

  result =  as.numeric(unlist(regmatches(result.gecr,gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",result.gecr, perl=TRUE))))

  result.all[1] = result[1]

  result.all[2] = result[2]

  }

  if(!is.null(result.heri)){

    result =  as.numeric(unlist(regmatches(result.heri,
                                           gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",
                                                    result.heri, perl=TRUE))))
    result.all[3] = result[2]

    result.all[4] = result[3]

    result.all[5] = result[5]

    result.all[6] = result[6]

  }

  if(!is.null(result.gecov)){

    result =  as.numeric(unlist(regmatches(result.gecov,
                                           gregexpr("-?\\ *[0-9]+\\.?[0-9]*(?:[Ee]\\ *-?\\ *[0-9]+)?",
                                                    result.gecov, perl=TRUE))))
    result.all[7] = result[1]

    result.all[8] = result[2]

  }


  #GR.LDSC.EST = LDSC.EST(Z1 = Z1, N1 = N1, Z2 = Z2, N2 = N2,
   #                      LDSCORE = LDSCORE, CSTR = CSTR, Nc=0, W=NULL)

  #Jackknife.Est = mclapply(1: length(Z1), function(i){

   # LDSCORE_i <- rowSums(SIGMA[-i,-i]^2)

    #LDSC.EST(Z1 = Z1[-i], N1 = N1, Z2 = Z2[-i], N2 = N2,
     #        LDSCORE = LDSCORE_i, CSTR = CSTR, Nc=0, W=NULL)

  #}, mc.cores = corenum)

 # GR.LDSC.SE = sd(unlist(Jackknife.Est))

  return(result.all)

}

