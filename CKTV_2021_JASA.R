###################################################################################
# Extrapolating TE in Multiple-cutoff RDDs
# Authors: Matias D. Cattaneo, Luke Keele, Rocio Titiunik and Gonzalo Vazquez-Bare
# Replication file
###################################################################################

rm(list = ls())

#####################################################################
## Load packages
#####################################################################

#install.packages('rdrobust')
#install.packages('nprobust')
#install.packages('rdmulti')


library(rdrobust)
library(nprobust)
library(rdmulti)


#####################################################################
## Setup
#####################################################################

data = read.csv('CKTV_2021_JASA.csv')

Y = data$ingresa_u3       # Outcome variable
X = data$icfes_puesto     # Running variable
C = data$cutoff           # Cutoff variable
D = as.numeric(X>=C)      # Treatment variable
X.norm = X - C            # Normalized running variable

c0 = -850                 # Low-cutoff group
c1 = -571                 # High-cutoff group
cc = -650                 # Evaluation point

#####################################################################
## Pooled and cutoff-specific effects
#####################################################################

Results = matrix(NA,11,7) # Matrix to collect results for Table 1

pooled = rdmc(Y,X,C)

Results[1,] = c(pooled$Coefs[1],
                pooled$H[1,1],
                pooled$Nh[1,1],
                pooled$V[1]^(1/2),
                pooled$Pv[1],
                pooled$CI[1,1],
                pooled$CI[2,1])

Results[2,] = c(pooled$Coefs[2],
                pooled$H[1,2],
                pooled$Nh[1,2],
                pooled$V[2]^(1/2),
                pooled$Pv[2],
                pooled$CI[1,2],
                pooled$CI[2,2])

Results[3,] = c(pooled$Coefs[3],
                pooled$H[1,3],
                pooled$Nh[1,3],
                pooled$V[3]^(1/2),
                pooled$Pv[3],
                pooled$CI[1,3],
                pooled$CI[2,3])

Results[4,] = c(pooled$Coefs[4],
                pooled$H[1,4],
                pooled$Nh[1,4],
                pooled$V[4]^(1/2),
                pooled$Pv[4],
                pooled$CI[1,4],
                pooled$CI[2,4])

#####################################################################
## Extrapolation: c0 to cc
#####################################################################

# Low-cutoff group: estimates at c0

aux1 = lprobust(Y[C==c0 & D==0],X[C==c0 & D==0],eval=c0)
h00 = aux1$Estimate[2]
N00 = aux1$Estimate[4]
m00 = aux1$Estimate[5]
m00.bc = aux1$Estimate[6]
V00 = aux1$Estimate[8]^2

Results[8,] = c(m00,h00,N00,NA,NA,NA,NA)


# Low-cutoff group: estimates at cc

aux2 = lprobust(Y[C==c0 & D==1],X[C==c0 & D==1],eval=cc)
h01 = aux2$Estimate[2]
N01 = aux2$Estimate[4]
m01 = aux2$Estimate[5]
m01.bc = aux2$Estimate[6]
V01 = aux2$Estimate[8]^2

Results[5,] = c(m01,h01,N01,NA,NA,NA,NA)


# High-cutoff group: estimates at c0 and cc (with covariance)

aux3 = lprobust(Y[C==c1 & D==0],X[C==c1 & D==0],eval=c(c0,cc),covgrid=TRUE,bwselect='mse-dpi')
h10 = aux3$Estimate[1,2]
N10 = aux3$Estimate[1,4]
m10 = aux3$Estimate[1,5]
m10.bc = aux3$Estimate[1,6]
V10 = aux3$Estimate[1,8]^2
h11 = aux3$Estimate[2,2]
N11 = aux3$Estimate[2,4]
m11 = aux3$Estimate[2,5]
m11.bc = aux3$Estimate[2,6]
V11 = aux3$Estimate[2,8]^2
cov.rb = aux3$cov.rb[1,2]

Results[9,] = c(m10,h10,N10,NA,NA,NA,NA)
Results[6,] = c(m11,h11,N11,NA,NA,NA,NA)


# TE: extrapolation Estimates

B = m00 - m10                                # Bias at low cutoff
Dif = m01- m11                               # Naive difference at evaluation point
TE = Dif - B                                 # Extrapolated TE estimate 
V.rb = V00 + V01 + V10 + V11 - 2*cov.rb      # Extrapolated TE RBC variance

TE.bc = m01.bc - m11.bc - m00.bc + m10.bc    # Extrapolated TE RBC estimate 
Tstat.bc = abs(TE.bc/sqrt(V.rb))             # Extrapolated TE RBC t-stat
pv.rb = 2*(1-pnorm(Tstat.bc))                # Extrapolated TE RBC p-value
CI.l = TE.bc - qnorm(0.975)*sqrt(V.rb)       # Extrapolated TE CI (lower limit)
CI.r = TE.bc + qnorm(0.975)*sqrt(V.rb)       # Extrapolated TE CI (upper limit)

Results[7,] = c(m01-m11,NA,NA,
                sqrt(V01+V11),
                2*(1-pnorm(abs((m01.bc-m11.bc)/sqrt(V01+V11)))),
                m01.bc-m11.bc - qnorm(0.975)*sqrt(V01+V11),
                m01.bc-m11.bc + qnorm(0.975)*sqrt(V01+V11))

Results[10,] = c(m00-m10,NA,NA,
                 sqrt(V00+V10),
                 2*(1-pnorm(abs((m00.bc-m10.bc)/sqrt(V00+V10)))),
                 m00.bc-m10.bc - qnorm(0.975)*sqrt(V00+V10),
                 m00.bc-m10.bc + qnorm(0.975)*sqrt(V00+V10))

Results[11,] = c(TE,NA,NA,sqrt(V.rb),pv.rb,CI.l,CI.r)

round(Results[,-4],3)

#####################################################################
## Parallel trends tests
#####################################################################

## Global approach: regression over X<=c0

library(car)
library(sandwich)

Pgl = matrix(NA,ncol=3,nrow=8)

Xsub = X[X<=c0]
Ysub = Y[X<=c0]
Csub = C[X<=c0]

reg =lm(Ysub ~ Xsub + I(Xsub^2) + (Csub==c1) + (Csub==c1)*(Xsub + I(Xsub^2)))

Pgl[1:6,1] = reg$coef
Pgl[1:6,2] = sqrt(diag(vcovHC(reg)))
Pgl[1:6,3] = 2*(1-pnorm(abs(Pgl[1:6,1]/Pgl[1:6,2])))

testcoef = linearHypothesis(reg,c("Xsub:Csub == c1TRUE = 0","I(Xsub^2):Csub == c1TRUE = 0"), test="Chisq",white.adjust=TRUE)
pvalue = testcoef$'Pr(>Chisq)'
pvalue.gl = pvalue[2]
pvalue.gl

Pgl[7,1] = length(resid(reg))
Pgl[8,3] = pvalue.gl

## Local approach: derivatives at c0

Ploc = matrix(NA,nrow=3,ncol=5)

deriv0 = lprobust(Y[C==c0 & D==0],X[C==c0 & D==0],eval=c0,deriv=1,p=2)
tau0 = deriv0$Estimate[5]
tau0.bc = deriv0$Estimate[6]
se0 = deriv0$Estimate[8]
Ploc[1,1] = tau0
Ploc[1,2] = deriv0$Estimate[2]
Ploc[1,3] = 2*(1-pnorm(abs(tau0.bc)/se0))
Ploc[1,4] = tau0.bc - qnorm(0.975)*se0
Ploc[1,5] = tau0.bc + qnorm(0.975)*se0

deriv1 = lprobust(Y[C==c1 & D==0],X[C==c1 & D==0],eval=c0,deriv=1,p=2)
tau1 = deriv1$Estimate[5]
tau1.bc = deriv1$Estimate[6]
se1 = deriv1$Estimate[8]
Ploc[2,1] = tau1
Ploc[2,2] = deriv1$Estimate[2]
Ploc[2,3] = 2*(1-pnorm(abs(tau1.bc)/se1))
Ploc[2,4] = tau1.bc - qnorm(0.975)*se1
Ploc[2,5] = tau1.bc + qnorm(0.975)*se1

Ploc[3,1] = tau0-tau1
Ploc[3,3] = 2*(1-pnorm(abs(tau0.bc-tau1.bc)/sqrt(se0^2+se1^2)))
Ploc[3,4] = tau0.bc-tau1.bc - qnorm(0.975)*sqrt(se0^2+se1^2)
Ploc[3,5] = tau0.bc-tau1.bc + qnorm(0.975)*sqrt(se0^2+se1^2)

Pgl
Ploc

#####################################################################
#####################################################################
## Simulations
#####################################################################
#####################################################################

rm(list = ls())

#install.packages('rdrobust')
#install.packages('nprobust')
#install.packages('rdmulti')
#install.packages('doParallel')
#install.packages('foreach')

library(doParallel)
library(foreach)
library(nprobust)
library(rdmulti)
library(rdrobust)


#####################################################################
## Setup and DGP
#####################################################################

coefs = c(-1.408947e+01,-7.395146e-02,-1.372308e-04,-1.123956e-07,-3.444331e-11)

#N = 1000 
#N = 2000
N = 5000
n.simul = 10000

clow = -850
chigh = -571
ceval = -650

set.seed(121919)

### Sample sizes

pr.clow = 1/2
N.clow = floor(N*pr.clow)
N.chigh = N - floor(N*pr.clow)

### Parameters and average potential outcomes

tau = 0.19
delta = -0.14

mu.10 = sum(c(1,poly(clow,degree=4,raw=TRUE))*coefs)
mu.00 = mu.10 + delta

mu.11 = sum(c(1,poly(ceval,degree=4,raw=TRUE))*coefs)
mu.01 = sum(c(1,poly(ceval,degree=4,raw=TRUE))*coefs) + delta + tau


parms = c(mu.00,mu.01,mu.10,mu.11,tau)

#####################################################################
## Simulations
#####################################################################

registerDoParallel(cores=detectCores(all.tests=TRUE)-1)

comb <- function(x, ...) {
  lapply(seq_along(x),function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))) # function to combine results from foreach into a list
}

time = proc.time()

results = foreach(i=1:n.simul,.combine='comb',.packages='nprobust') %dopar% {
  
  X.simul = runif(N,-1000,-1)
  Px = poly(X.simul,degree=4,raw=TRUE)
  Px = cbind(rep(1,nrow(Px)),Px)
  EY0.simul = Px%*%coefs
  Y0.simul = EY0.simul + rnorm(N,sd=.3)
  Y1.simul = Y0.simul + tau
  C.simul = sample(c(rep(clow,N.clow),rep(chigh,N.chigh)),N.clow+N.chigh,replace=FALSE)
  D.simul = X.simul>=C.simul
  Y.simul = Y0.simul + tau*D.simul + delta*(C.simul==clow)
  
  aux1 = lprobust(Y.simul[C.simul==clow & D.simul==0],X.simul[C.simul==clow & D.simul==0],eval=clow)
  N00 = aux1$Estimate[4]
  m00 = aux1$Estimate[5]
  m00.bc = aux1$Estimate[6]
  V00.rb = aux1$Estimate[8]^2
  
  aux2 = lprobust(Y.simul[C.simul==clow & D.simul==1],X.simul[C.simul==clow & D.simul==1],eval=ceval)
  N01 = aux2$Estimate[4]
  m01 = aux2$Estimate[5]
  m01.bc = aux2$Estimate[6]
  V01.rb = aux2$Estimate[8]^2
  
  aux3 = lprobust(Y.simul[C.simul==chigh & D.simul==0],X.simul[C.simul==chigh & D.simul==0],
                  eval=c(clow,ceval),covgrid=TRUE,bwselect='mse-dpi')
  N10 = aux3$Estimate[1,4]
  m10 = aux3$Estimate[1,5]
  m10.bc = aux3$Estimate[1,6]
  V10.rb = aux3$Estimate[1,8]^2
  N11 = aux3$Estimate[2,4]
  m11 = aux3$Estimate[2,5]
  m11.bc = aux3$Estimate[2,6]
  V11.rb = aux3$Estimate[2,8]^2
  cov.rb = aux3$cov.rb[1,2]
  
  B = m00 - m10
  Dif = m01- m11
  TE = Dif - B
  V.rb = V00.rb + V01.rb + V10.rb + V11.rb - 2*cov.rb
  
  TE.bc = m01.bc - m11.bc - m00.bc + m10.bc
  Tstat.bc = abs(TE.bc/sqrt(V.rb))
  
  Beta.cl.vec = c(m00,m01,m10,m11,TE)
  Beta.bc.vec = c(m00.bc,m01.bc,m10.bc,m11.bc,TE.bc)
  Var.rb.vec = c(V00.rb,V01.rb,V10.rb,V11.rb,V.rb)
  NH.vec = c(N00,N01,N10,N11)
  
  list(Beta.cl.vec,
       Beta.bc.vec,
       Var.rb.vec,
       NH.vec)
}

proc.time() - time

stopImplicitCluster()

Beta.cl = matrix(unlist(results[[1]]),nrow=5)
Beta.bc = matrix(unlist(results[[2]]),nrow=5)
Var.rb = matrix(unlist(results[[3]]),nrow=5)
NH = matrix(unlist(results[[4]]),nrow=4)

CI.l = Beta.bc - qnorm(0.975)*sqrt(Var.rb)
CI.r = Beta.bc + qnorm(0.975)*sqrt(Var.rb)

R.simul = matrix(NA,nrow=5,ncol=5)
R.simul[1:4,1] = rowMeans(NH)
R.simul[,2] = rowMeans(Beta.cl-parms)
R.simul[,3] = apply(Beta.cl,1,sd)^2
R.simul[,4] = sqrt(rowMeans((Beta.cl-parms)^2))
R.simul[,5] = rowMeans((parms>=CI.l & parms<=CI.r))

colnames(R.simul) = c("Nh","Bias","Variance","RMSE","Coverage")

round(R.simul,4)

#####################################################################
#####################################################################
### SUPPLEMENTAL APPENDIX
#####################################################################
#####################################################################

rm(list = ls())
library(rdrobust)
library(nprobust)
library(rdmulti)

data = read.csv('CKTV_2021_JASA_SA.csv')

Y = data$ingresa_u3       
X = data$icfes_puesto     
C = data$cutoff           
D = as.numeric(X>=C)      
X.norm = X - C            

c0 = -850                 
c1 = -571                 
cc = -650                 

#####################################################################
## Pooled and cutoff-specific effects
#####################################################################

Results = matrix(NA,11,7) 

pooled = rdmc(Y,X,C)

Results[1,] = c(pooled$Coefs[1],
                pooled$H[1,1],
                pooled$Nh[1,1],
                pooled$V[1]^(1/2),
                pooled$Pv[1],
                pooled$CI[1,1],
                pooled$CI[2,1])

Results[2,] = c(pooled$Coefs[2],
                pooled$H[1,2],
                pooled$Nh[1,2],
                pooled$V[2]^(1/2),
                pooled$Pv[2],
                pooled$CI[1,2],
                pooled$CI[2,2])

Results[3,] = c(pooled$Coefs[3],
                pooled$H[1,3],
                pooled$Nh[1,3],
                pooled$V[3]^(1/2),
                pooled$Pv[3],
                pooled$CI[1,3],
                pooled$CI[2,3])

Results[4,] = c(pooled$Coefs[4],
                pooled$H[1,4],
                pooled$Nh[1,4],
                pooled$V[4]^(1/2),
                pooled$Pv[4],
                pooled$CI[1,4],
                pooled$CI[2,4])

#####################################################################
## Extrapolation: c0 to cc
#####################################################################

# Low-cutoff group: estimates at c0

aux1 = lprobust(Y[C==c0 & D==0],X[C==c0 & D==0],eval=c0)
h00 = aux1$Estimate[2]
N00 = aux1$Estimate[4]
m00 = aux1$Estimate[5]
m00.bc = aux1$Estimate[6]
V00 = aux1$Estimate[8]^2

Results[8,] = c(m00,h00,N00,NA,NA,NA,NA)


# Low-cutoff group: estimates at cc

aux2 = lprobust(Y[C==c0 & D==1],X[C==c0 & D==1],eval=cc)
h01 = aux2$Estimate[2]
N01 = aux2$Estimate[4]
m01 = aux2$Estimate[5]
m01.bc = aux2$Estimate[6]
V01 = aux2$Estimate[8]^2

Results[5,] = c(m01,h01,N01,NA,NA,NA,NA)


# High-cutoff group: estimates at c0 and cc (with covariance)

aux3 = lprobust(Y[C==c1 & D==0],X[C==c1 & D==0],eval=c(c0,cc),covgrid=TRUE,bwselect='mse-dpi')
h10 = aux3$Estimate[1,2]
N10 = aux3$Estimate[1,4]
m10 = aux3$Estimate[1,5]
m10.bc = aux3$Estimate[1,6]
V10 = aux3$Estimate[1,8]^2
h11 = aux3$Estimate[2,2]
N11 = aux3$Estimate[2,4]
m11 = aux3$Estimate[2,5]
m11.bc = aux3$Estimate[2,6]
V11 = aux3$Estimate[2,8]^2
cov.rb = aux3$cov.rb[1,2]

Results[9,] = c(m10,h10,N10,NA,NA,NA,NA)
Results[6,] = c(m11,h11,N11,NA,NA,NA,NA)


# TE: extrapolation Estimates

B = m00 - m10                               
Dif = m01- m11                              
TE = Dif - B                                
V.rb = V00 + V01 + V10 + V11 - 2*cov.rb     

TE.bc = m01.bc - m11.bc - m00.bc + m10.bc    
Tstat.bc = abs(TE.bc/sqrt(V.rb))             
pv.rb = 2*(1-pnorm(Tstat.bc))               
CI.l = TE.bc - qnorm(0.975)*sqrt(V.rb)       
CI.r = TE.bc + qnorm(0.975)*sqrt(V.rb)       

Results[7,] = c(m01-m11,NA,NA,
                sqrt(V01+V11),
                2*(1-pnorm(abs((m01.bc-m11.bc)/sqrt(V01+V11)))),
                m01.bc-m11.bc - qnorm(0.975)*sqrt(V01+V11),
                m01.bc-m11.bc + qnorm(0.975)*sqrt(V01+V11))

Results[10,] = c(m00-m10,NA,NA,
                 sqrt(V00+V10),
                 2*(1-pnorm(abs((m00.bc-m10.bc)/sqrt(V00+V10)))),
                 m00.bc-m10.bc - qnorm(0.975)*sqrt(V00+V10),
                 m00.bc-m10.bc + qnorm(0.975)*sqrt(V00+V10))

Results[11,] = c(TE,NA,NA,sqrt(V.rb),pv.rb,CI.l,CI.r)

round(Results[,-4],3)

#####################################################################
## Local randomization
#####################################################################

rm(list = ls())

data = read.csv('CKTV_2021_JASA.csv')
Y = data$ingresa_u3       # Outcome variable
X = data$icfes_puesto     # Running variable
C = data$cutoff           # Cutoff variable
D = as.numeric(X>=C)      # Treatment variable
X.norm = X - C            # Normalized running variable

c0 = -850                 # Low-cutoff group
c1 = -571                 # High-cutoff group
cc = -650                 # Evaluation point

findw = function(X,nobs){
  Xaux = sort(X)
  wlength.x = as.numeric(which.max(cumsum(table(Xaux))>=nobs))
  tmp = merge(unique(X),table(X),by=1)
  w = tmp[wlength.x,1]
  return(w)
}

beta = .01    # 1 - level of CI for Delta
lgrid = 20    # Length of grid for Delta

perms = 5000
set.seed(20190909)

Res.LR = matrix(NA,nrow=7,ncol=9)

nobs = 50                         # desired number of obs in window

#### Estimation of Delta (low cutoff)

## Windows and sample sizes

Xaux = X[C==c0 & D==0]
w.low0 = findw(abs(Xaux),nobs)

Xaux = X[C==c1 & D==0 & X<c0]
w.low1L = findw(abs(Xaux),nobs/2)

Xaux = X[C==c1 & D==0 & X>=c0]
w.low1R = findw(Xaux,nobs/2)

n0 = length(Y[C==c0 & D==0 & X>=-w.low0])
n1 = length(Y[C==c1 & D==0 & X>=-w.low1L & X<=w.low1R])

Res.LR[1,1] = -w.low0
Res.LR[1,2] = c0
Res.LR[2,1] = -w.low1L
Res.LR[2,2] = w.low1R
Res.LR[1,3] = n0
Res.LR[2,3] = n1
Res.LR[3,3] = n0+n1

Yaux.l = Y[(C==c0 & D==0 & X>=-w.low0) | (C==c1 & D==0 & X>=-w.low1L & X<=w.low1R)]
Xaux.l = X[(C==c0 & D==0 & X>=-w.low0) | (C==c1 & D==0 & X>=-w.low1L & X<=w.low1R)] - c0
Caux.l = C[(C==c0 & D==0 & X>=-w.low0) | (C==c1 & D==0 & X>=-w.low1L & X<=w.low1R)]
Daux.l = as.numeric(Caux.l==c0)

## Difference in means

Delta = mean(Yaux.l[Daux.l==1])-mean(Yaux.l[Daux.l==0])

Res.LR[1,4] = mean(Yaux.l[Daux.l==1])
Res.LR[2,4] = mean(Yaux.l[Daux.l==0])
Res.LR[3,4] = Delta

var.Delta = var(Yaux.l[Daux.l==1])/length(Yaux.l[Daux.l==1])+var(Yaux.l[Daux.l==0])/length(Yaux.l[Daux.l==0])
Delta.grid = matrix(seq(from=Delta-qnorm(1-beta/2)*sqrt(var.Delta),
                        to=Delta+qnorm(1-beta/2)*sqrt(var.Delta),length.out=lgrid),nrow=1)

## Linear adjustment 

reg.l = lm(Yaux.l ~ Daux.l*Xaux.l)
Delta.lin = reg.l$coef['Daux.l']

Res.LR[1,7] = reg.l$coef['(Intercept)'] + Delta.lin
Res.LR[2,7] = reg.l$coef['(Intercept)']
Res.LR[3,7] = Delta.lin

var.Delta.lin = vcov(reg.l,type='HC0')['Daux.l','Daux.l']
Delta.grid.lin = matrix(seq(from=Delta.lin-qnorm(1-beta/2)*sqrt(var.Delta.lin),
                            to=Delta.lin+qnorm(1-beta/2)*sqrt(var.Delta.lin),length.out=lgrid),nrow=1)

## Randomization inference on Delta

Yresid.l = reg.l$resid + reg.l$coef['(Intercept)'] + reg.l$coef['Daux.l']*Daux.l

deltavec = numeric(perms)
deltalinvec = numeric(perms)
for (k in 1:perms){
  dperm = sample(Daux.l,replace=FALSE)
  deltavec[k] = mean(Yaux.l[dperm==1])-mean(Yaux.l[dperm==0])
  deltalinvec[k]=mean(Yresid.l[dperm==1])-mean(Yresid.l[dperm==0])
}

Res.LR[3,5] = mean(abs(deltavec)>=abs(Delta))
Res.LR[3,8] = mean(abs(deltalinvec)>=abs(Delta.lin))

## Neyman inference on Delta

Res.LR[3,6] = 2*pnorm(-abs(Delta/sqrt(var.Delta)))
Res.LR[3,9] = 2*pnorm(-abs(Delta.lin/sqrt(var.Delta.lin)))

#### Estimation of extrapolated TE

## Windows and sample sizes

Xaux = X[C==c0 & D==1 & X<cc]
w.high0L = findw(abs(Xaux),nobs/2)

Xaux = X[C==c0 & D==1 & X>=cc]
w.high0R = findw(Xaux,nobs/2)

Xaux = X[C==c1 & D==0 & X<cc]
w.high1L = findw(abs(Xaux),nobs/2)

Xaux = X[C==c1 & D==0 & X>=cc]
w.high1R = findw(Xaux,nobs/2)

N0 = length(Y[C==c0 & D==1 & X>=-w.high0L & X<=w.high0R])
N1 = length(Y[C==c1 & D==0 & X>=-w.high1L & X<=w.high1R])

Res.LR[4,1] = -w.high0L
Res.LR[4,2] = w.high0R
Res.LR[4,3] = N0
Res.LR[5,1] = -w.high1L
Res.LR[5,2] = w.high1R
Res.LR[5,3] = N1
Res.LR[6,3] = N0+N1

Yaux.h = Y[(C==c0 & D==1 & X>=-w.high0L & X<=w.high0R) | (C==c1 & D==0 & X>=-w.high1L & X<=w.high1R)]
Xaux.h = X[(C==c0 & D==1 & X>=-w.high0L & X<=w.high0R) | (C==c1 & D==0 & X>=-w.high1L & X<=w.high1R)] - cc
Caux.h = C[(C==c0 & D==1 & X>=-w.high0L & X<=w.high0R) | (C==c1 & D==0 & X>=-w.high1L & X<=w.high1R)]
Daux.h = as.numeric(Caux.h==c0)

# Adjusted outcomes

Yadj = Yaux.h + Delta*(1-Daux.h)
Yadj.lin = Yaux.h + Delta.lin*(1-Daux.h)

Ymat = apply(Delta.grid,2,function(z) Yaux.h+z*(1-Daux.h))
Ymat.lin = apply(Delta.grid.lin,2,function(z) Yaux.h+z*(1-Daux.h))

## Difference in means

tauhat = mean(Yadj[Daux.h==1])-mean(Yadj[Daux.h==0])

Res.LR[4,4] = mean(Yadj[Daux.h==1])
Res.LR[5,4] = mean(Yadj[Daux.h==0])-Delta
Res.LR[6,4] = tauhat+Delta
Res.LR[7,4] = tauhat

tauhat.sup = apply(Ymat,2,function(x) mean(x[Daux.h==1])-mean(x[Daux.h==0]))

## Linear adjustment

reg.naive = lm(Yaux.h ~ Daux.h*Xaux.h)

reg.h = lm(Yadj.lin ~ Daux.h*Xaux.h)
tauhat.lin = reg.h$coef['Daux.h']

Res.LR[4,7] = reg.h$coef['(Intercept)'] + tauhat.lin
Res.LR[5,7] = reg.h$coef['(Intercept)'] - Delta.lin
Res.LR[6,7] = tauhat.lin + Delta.lin
Res.LR[7,7] = tauhat.lin

tauhat.sup.lin = apply(Ymat,2,function(x) lm(x~Xaux.h*Daux.h)$coef['Daux.h'])

Yresid.h = reg.naive$resid + reg.naive$coef['(Intercept)'] + reg.naive$coef['Daux.h']*Daux.h
Ymat.resid = apply(Ymat.lin,2,function(x) lm(x~Xaux.h*Daux.h)$resid + lm(x~Xaux.h*Daux.h)$coef['(Intercept)']+lm(x~Xaux.h*Daux.h)$coef['Daux.h']*Daux.h)

#### Randomization inference for TE

coefmat = matrix(NA,nrow=length(Delta.grid),ncol=perms)
coefmat.lin = matrix(NA,nrow=length(Delta.grid.lin),ncol=perms)
coefvec.naive = numeric(perms)
coefvec.naive.lin = numeric(perms)
for (k in 1:perms){
  
  dperm = sample(Daux.h,replace=FALSE)
  coefvec.naive[k] = mean(Yaux.h[dperm==1])-mean(Yaux.h[dperm==0])
  coefvec.naive.lin[k] = mean(Yresid.h[dperm==1])-mean(Yresid.h[dperm==0])
  coefmat[,k] = apply(Ymat,2,function(x) mean(x[dperm==1]) - mean(x[dperm==0]))
  coefmat.lin[,k] = apply(Ymat.resid,2,function(x) mean(x[dperm==1]) - mean(x[dperm==0]))
  
}


Res.LR[6,5] = mean(abs(coefvec.naive)>=abs(Res.LR[6,4]))
Res.LR[6,8] = mean(abs(coefvec.naive.lin)>=abs(Res.LR[6,7]))
Res.LR[7,5] = max(rowMeans(abs(coefmat)>=abs(tauhat.sup))) + beta
Res.LR[7,8] = max(rowMeans(abs(coefmat.lin)>=abs(tauhat.sup.lin))) + beta

# Neyman p-values

var.naive = var(Yaux.h[Daux.h==1])/length(Yaux.h[Daux.h==1])+var(Yaux.h[Daux.h==0])/length(Yaux.h[Daux.h==0])
var.naive.lin = vcov(reg.naive,type='HC0')['Daux.h','Daux.h']
var.tauhat = var(Yadj[Daux.h==1])/length(Yadj[Daux.h==1])+var(Yadj[Daux.h==0])/length(Yadj[Daux.h==0])
var.tauhat.lin = vcov(reg.h,type='HC0')['Daux.h','Daux.h']

Res.LR[6,6] = 2*pnorm(-abs(Res.LR[6,4]/sqrt(var.naive)))
Res.LR[6,9] = 2*pnorm(-abs(Res.LR[6,7]/sqrt(var.naive.lin)))
Res.LR[7,6] = 2*(pnorm(-abs(tauhat / sqrt(var.tauhat + var.Delta))))
Res.LR[7,9] = 2*pnorm(-abs(tauhat.lin / sqrt(var.tauhat.lin + var.Delta.lin)))

# Results

round(Res.LR,3)

