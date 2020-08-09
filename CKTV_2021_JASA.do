********************************************************************************
** Extrapolating TE in Multiple-cutoff RDDs
** Authors: Matias D. Cattaneo, Luke Keele, Rocio Titiunik and Gonzalo Vazquez-Bare
** Replication file
********************************************************************************

use CKTV_2021_JASA.dta, clear

gen double Y = ingresa_u3
gen double X = icfes_puesto
gen double C = cutoff
gen double D = X>=C
gen double Xnorm = X-C

global c0 = -850
global c1 = -571
global cc = -650

********************************************************************************
** Pooled and cutoff-specific effects
********************************************************************************

mat Results = J(11,7,.)

rdmc Y X, c(C)

mat Coefs = e(coefs)
mat H = e(H)
mat Nh = e(sampsis)
mat V = e(V)
mat P = e(pv_rb)
mat CI = e(CI_rb)

mat Results[1,1] = Coefs[1,1]
mat Results[1,2] = H[1,1]
mat Results[1,3] = Nh[1,1]
mat Results[1,4] = sqrt(V[1,1])
mat Results[1,5] = P[1,1]
mat Results[1,6] = CI[1,1]
mat Results[1,7] = CI[2,1]

mat Results[2,1] = Coefs[1,2]
mat Results[2,2] = H[1,2]
mat Results[2,3] = Nh[1,2]
mat Results[2,4] = sqrt(V[2,2])
mat Results[2,5] = P[1,2]
mat Results[2,6] = CI[1,2]
mat Results[2,7] = CI[2,2]

mat Results[3,1] = Coefs[1,3]
mat Results[3,2] = H[1,3]
mat Results[3,3] = Nh[1,3]
mat Results[3,4] = sqrt(V[3,3])
mat Results[3,5] = P[1,3]
mat Results[3,6] = CI[1,3]
mat Results[3,7] = CI[2,3]

mat Results[4,1] = Coefs[1,4]
mat Results[4,2] = H[1,4]
mat Results[4,3] = Nh[1,4]
mat Results[4,4] = sqrt(V[4,4])
mat Results[4,5] = P[1,4]
mat Results[4,6] = CI[1,4]
mat Results[4,7] = CI[2,4]

********************************************************************************
** Extrapolation: c0 to cc
********************************************************************************

* Low-cutoff group: estimates at c0

gen eval0 = $c0 in 1
lprobust Y X if C==$c0 & D==0, eval(eval0)

mat aux = e(Result)
gl h00 = aux[1,2]
gl N00 = aux[1,4]
gl m00 = aux[1,5]
gl m00_bc = aux[1,6]
gl V00 = aux[1,8]^2

mat Results[8,1] = $m00
mat Results[8,2] = $h00
mat Results[8,3] = $N00

* Low-cutoff group: estimates at cc

gen eval1 = $cc in 1
lprobust Y X if C==$c0 & D==1, eval(eval1)

mat aux = e(Result)
gl h01 = aux[1,2]
gl N01 = aux[1,4]
gl m01 = aux[1,5]
gl m01_bc = aux[1,6]
gl V01 = aux[1,8]^2

mat Results[5,1] = $m01
mat Results[5,2] = $h01
mat Results[5,3] = $N01

* High-cutoff group: estimates at c0 and cc (with covariance)

gen eval2 = $c0 in 1
replace eval2 = $cc in 2
lprobust Y X if C==$c1 & D==0, eval(eval2) covgrid bwselect("mse-dpi")

mat aux = e(Result)
mat Cov = e(cov_rb)
gl h10 = aux[1,2]
gl N10 = aux[1,4]
gl m10 = aux[1,5]
gl m10_bc = aux[1,6]
gl V10 = aux[1,8]^2
gl h11 = aux[2,2]
gl N11 = aux[2,4]
gl m11 = aux[2,5]
gl m11_bc = aux[2,6]
gl V11 = aux[2,8]^2
gl cov_rb = Cov[2,1]

mat Results[9,1] = $m10
mat Results[9,2] = $h10
mat Results[9,3] = $N10 
mat Results[6,1] = $m11
mat Results[6,2] = $h11
mat Results[6,3] = $N11

* TE: extrapolation estimates

gl B = $m00 - $m10 								 // Bias at low cutoff
gl Dif = $m01 - $m11							 // Naive difference at eval point
gl TE = $Dif - $B                                // Extrapolated TE estimate
gl V_rb = $V00 + $V01 + $V10 + $V11 - 2*$cov_rb  // Extrapolated TE RBC variance

gl TE_bc = $m01_bc - $m11_bc - $m00_bc + $m10_bc // Extrapolated TE RBC estimate
gl Tstat_bc = abs($TE_bc/sqrt($V_rb))            // Extrapolated TE RBC  t-stat
gl pv_rb = 2*(1-normal($Tstat_bc))               // Extrapolated TE RBC p-value
gl CI_l = $TE_bc - invnormal(0.975)*sqrt($V_rb)  // Extrapolated TE CI (left)
gl CI_r = $TE_bc + invnormal(0.975)*sqrt($V_rb)  // Extrapolated TE CI (right)

mat Results[7,1] = $m01 - $m11
mat Results[7,4] = sqrt($V01+$V11)
mat Results[7,5] = 2*(1-normal(abs($m01_bc-$m11_bc)/sqrt($V01+$V11)))
mat Results[7,6] = $m01_bc-$m11_bc - invnormal(0.975)*sqrt($V01+$V11)
mat Results[7,7] = $m01_bc-$m11_bc + invnormal(0.975)*sqrt($V01+$V11)

mat Results[10,1] = $m00 - $m10
mat Results[10,4] = sqrt($V00+$V10)
mat Results[10,5] = 2*(1-normal(abs($m00_bc-$m10_bc)/sqrt($V00+$V10)))
mat Results[10,6] = $m00_bc-$m10_bc - invnormal(0.975)*sqrt($V00+$V10)
mat Results[10,7] = $m00_bc-$m10_bc + invnormal(0.975)*sqrt($V00+$V10)

mat Results[11,1] = $TE
mat Results[11,4] = sqrt($V_rb)
mat Results[11,5] = $pv_rb
mat Results[11,6] = $CI_l
mat Results[11,7] = $CI_r

matlist Results

********************************************************************************
** Parallel trends tests
********************************************************************************

** Global approach: regression over X<=c0

gen C1 = C==$c1
reg Y i.C1##(c.X##c.X) if X<=$c0, vce(hc3)

** Local approach: derivatives at c0

mat Ploc = J(3,5,.)

lprobust Y X if C==$c0 & D==0, eval(eval0) deriv(1) p(2)
mat aux = e(Result)
gl tau0 = aux[1,5]
gl tau0_bc = aux[1,6]
gl se0 = aux[1,8]

mat Ploc[1,1] = $tau0
mat Ploc[1,2] = aux[1,2]
mat Ploc[1,3] = 2*(1-normal(abs($tau0_bc/$se0)))
mat Ploc[1,4] = $tau0_bc - invnormal(0.975)*$se0
mat Ploc[1,5] = $tau0_bc + invnormal(0.975)*$se0

lprobust Y X if C==$c1 & D==0, eval(eval0) deriv(1) p(2)
mat aux = e(Result)
gl tau1 = aux[1,5]
gl tau1_bc = aux[1,6]
gl se1 = aux[1,8]

mat Ploc[2,1] = $tau1
mat Ploc[2,2] = aux[1,2]
mat Ploc[2,3] = 2*(1-normal(abs($tau1_bc/$se1)))
mat Ploc[2,4] = $tau1_bc - invnormal(0.975)*$se1
mat Ploc[2,5] = $tau1_bc + invnormal(0.975)*$se1

mat Ploc[3,1] = $tau0-$tau1
mat Ploc[3,3] = 2*(1-normal(abs(($tau0_bc-$tau1_bc)/sqrt($se0^2+$se1^2))))
mat Ploc[3,4] = $tau0_bc-$tau1_bc - invnormal(0.975)*sqrt($se0^2+$se1^2)
mat Ploc[3,5] = $tau0_bc-$tau1_bc + invnormal(0.975)*sqrt($se0^2+$se1^2)

matlist Ploc

********************************************************************************
********************************************************************************
** Simulations
********************************************************************************
********************************************************************************

** See R code for replication

********************************************************************************
********************************************************************************
** Supplemental Appendix
********************************************************************************
********************************************************************************

use CKTV_2021_JASA_SA.dta, clear

gen double Y = ingresa_u3
gen double X = icfes_puesto
gen double C = cutoff
gen double D = X>=C
gen double Xnorm = X-C

global c0 = -850
global c1 = -571
global cc = -650

********************************************************************************
** Pooled and cutoff-specific effects
********************************************************************************

mat Results = J(11,7,.)

rdmc Y X, c(C)

mat Coefs = e(coefs)
mat H = e(H)
mat Nh = e(sampsis)
mat V = e(V)
mat P = e(pv_rb)
mat CI = e(CI_rb)

mat Results[1,1] = Coefs[1,1]
mat Results[1,2] = H[1,1]
mat Results[1,3] = Nh[1,1]
mat Results[1,4] = sqrt(V[1,1])
mat Results[1,5] = P[1,1]
mat Results[1,6] = CI[1,1]
mat Results[1,7] = CI[2,1]

mat Results[2,1] = Coefs[1,2]
mat Results[2,2] = H[1,2]
mat Results[2,3] = Nh[1,2]
mat Results[2,4] = sqrt(V[2,2])
mat Results[2,5] = P[1,2]
mat Results[2,6] = CI[1,2]
mat Results[2,7] = CI[2,2]

mat Results[3,1] = Coefs[1,3]
mat Results[3,2] = H[1,3]
mat Results[3,3] = Nh[1,3]
mat Results[3,4] = sqrt(V[3,3])
mat Results[3,5] = P[1,3]
mat Results[3,6] = CI[1,3]
mat Results[3,7] = CI[2,3]

mat Results[4,1] = Coefs[1,4]
mat Results[4,2] = H[1,4]
mat Results[4,3] = Nh[1,4]
mat Results[4,4] = sqrt(V[4,4])
mat Results[4,5] = P[1,4]
mat Results[4,6] = CI[1,4]
mat Results[4,7] = CI[2,4]

********************************************************************************
** Extrapolation: c0 to cc
********************************************************************************

* Low-cutoff group: estimates at c0

gen eval0 = $c0 in 1
lprobust Y X if C==$c0 & D==0, eval(eval0)

mat aux = e(Result)
gl h00 = aux[1,2]
gl N00 = aux[1,4]
gl m00 = aux[1,5]
gl m00_bc = aux[1,6]
gl V00 = aux[1,8]^2

mat Results[8,1] = $m00
mat Results[8,2] = $h00
mat Results[8,3] = $N00

* Low-cutoff group: estimates at cc

gen eval1 = $cc in 1
lprobust Y X if C==$c0 & D==1, eval(eval1)

mat aux = e(Result)
gl h01 = aux[1,2]
gl N01 = aux[1,4]
gl m01 = aux[1,5]
gl m01_bc = aux[1,6]
gl V01 = aux[1,8]^2

mat Results[5,1] = $m01
mat Results[5,2] = $h01
mat Results[5,3] = $N01

* High-cutoff group: estimates at c0 and cc (with covariance)

gen eval2 = $c0 in 1
replace eval2 = $cc in 2
lprobust Y X if C==$c1 & D==0, eval(eval2) covgrid bwselect("mse-dpi")

mat aux = e(Result)
mat Cov = e(cov_rb)
gl h10 = aux[1,2]
gl N10 = aux[1,4]
gl m10 = aux[1,5]
gl m10_bc = aux[1,6]
gl V10 = aux[1,8]^2
gl h11 = aux[2,2]
gl N11 = aux[2,4]
gl m11 = aux[2,5]
gl m11_bc = aux[2,6]
gl V11 = aux[2,8]^2
gl cov_rb = Cov[2,1]

mat Results[9,1] = $m10
mat Results[9,2] = $h10
mat Results[9,3] = $N10 
mat Results[6,1] = $m11
mat Results[6,2] = $h11
mat Results[6,3] = $N11

* TE: extrapolation estimates

gl B = $m00 - $m10 								 
gl Dif = $m01 - $m11							 
gl TE = $Dif - $B                                
gl V_rb = $V00 + $V01 + $V10 + $V11 - 2*$cov_rb                      

gl TE_bc = $m01_bc - $m11_bc - $m00_bc + $m10_bc 
gl Tstat_bc = abs($TE_bc/sqrt($V_rb))            
gl pv_rb = 2*(1-normal($Tstat_bc))               
gl CI_l = $TE_bc - invnormal(0.975)*sqrt($V_rb)  
gl CI_r = $TE_bc + invnormal(0.975)*sqrt($V_rb)  

mat Results[7,1] = $m01 - $m11
mat Results[7,4] = sqrt($V01+$V11)
mat Results[7,5] = 2*(1-normal(abs($m01_bc-$m11_bc)/sqrt($V01+$V11)))
mat Results[7,6] = $m01_bc-$m11_bc - invnormal(0.975)*sqrt($V01+$V11)
mat Results[7,7] = $m01_bc-$m11_bc + invnormal(0.975)*sqrt($V01+$V11)

mat Results[10,1] = $m00 - $m10
mat Results[10,4] = sqrt($V00+$V10)
mat Results[10,5] = 2*(1-normal(abs($m00_bc-$m10_bc)/sqrt($V00+$V10)))
mat Results[10,6] = $m00_bc-$m10_bc - invnormal(0.975)*sqrt($V00+$V10)
mat Results[10,7] = $m00_bc-$m10_bc + invnormal(0.975)*sqrt($V00+$V10)

mat Results[11,1] = $TE
mat Results[11,4] = sqrt($V_rb)
mat Results[11,5] = $pv_rb
mat Results[11,6] = $CI_l
mat Results[11,7] = $CI_r

matlist Results

********************************************************************************
** Local randomization
********************************************************************************

** See R code for replication

