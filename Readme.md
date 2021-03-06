Estimating the genetic correlation with GWAS summary association
statistics
================
Jianqiao Wang
2/27/2022

## Input

To estimate the genetic correlation, the proposed method requires the
information of summary association statistics and external LD matrix

### summary statistics

The summary statistics file should have the following information:

  - SNP: the name of the genetic variants
  - A1: the effect allele.
  - n: number of samples used when testing the predictor.
  - Z: The z test staitistics of the genetic variants.

Note that the Z statistics could be calculated based on estimate of
effect size (or log odds) of a predictor, and its standard deviation, or
converted from the p values along with effects size.

### external LD matrix

The external LD matrix could be calculated from external genotype data.
The corresponding variant name and the effect allele should also be kept

### Before the input

  - The summary statistics for outcome 1
  - The summary statistics for outcome 2
  - external LD information

We note theat two summary statistics should come from differet studies
so that they are independent with each other. **Make sure that all files
have the same genetic variants and effect alleles**

## implement the method

``` r
source('./R/LDScore.R')
source('./R/Genetic_Correlation_class_v2.R')
source('./R/SummaryStatistics.R')
source('./R/Precsion_Matrix_Estimation.R')
source('./R/Helper.R')

# read the sumstat file 
sumstat.beta = read.table("./data/FG_female.sumstats", header = T) 
sumstat.gamma = read.table("./data/HG_male.sumstats", header = T)
n = mean(sumstat.beta$N, na.rm = T)
m = mean(sumstat.gamma$N, na.rm = T)

# convert to beta_mg and gamma_mg
beta_mg = Z.to.Mg(sumstat.beta$Z, n = n)
gamma_mg = Z.to.Mg(sumstat.gamma$Z, n = m)


# Omega0, SIGMA #  LD matrix calculated from the external genotype data
SIGMA = readRDS("./data/SIGMA.rds")
Omega0 = readRDS("./data/Omega0.rds")


# Inititate the estimation object
GeCor.Analysis = GeCov$new(Beta_mg = beta_mg,
                             Gamma_mg = gamma_mg,
                             N1 = n,
                             N2 = m)

# plug-in the external Omega and Sigma
# Estimation 
GeCor.Analysis$Estimation(Omega = Omega0, Sigma = SIGMA, Type = "proposed")
# variance
GeCor.Analysis$Var.Est(TRUE)
# print out the result 
print(GeCor.Analysis)
```

If we ignore the dependence, the estimator is

``` r
ldscore = read.table("./data/simu.l2.ldscore", header = T)

# ignore the dependence
GeCor.Analysis$Estimation(Omega = I, Sigma = SIGMA, LDSCORE = ldscore$LD, Type = "Identity")

GeCor.Analysis$Var.Est(TRUE)

print(GeCor.Analysis)
```
