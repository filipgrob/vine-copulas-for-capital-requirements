#libraries
rm(list = ls())
set.seed(129)
library(stats)
library(VineCopula)
library(skewt)
#global parameters
T=1000
d = 7
N=1000
dd = d*(d-1)/2
mat=matrix(c(7,4,5,1,2,3,6, 0,4,6,5,1,2,3, 0,0,6,5,1,2,3, 
             0,0,0,5,1,3,2, 0,0,0,0,1,3,2, 0,0,0,0,0,3,2,
             0,0,0,0,0,0,2),7,7)











#MONTE CARLO SIMULATION
#SCENARIO WITH FIXED BB1 COPULA

#initializations
k=runif(7, 0, 7)
gamma=runif(7, 1, 7)
par1_mat=array(0,dim=c(d,d))
par1_mat[lower.tri(par1_mat)]=k
par2_mat=array(0,dim=c(d,d))
par2_mat[lower.tri(par2_mat)]=gamma
fam_mat=matrix(0,nrow=d,ncol=d)
fam_mat[lower.tri(fam_mat,diag=FALSE)]=7
rvm_r=RVineMatrix(Matrix=mat,family=fam_mat,par=par1_mat,par2=par2_mat)
rvine=array(NA, dim = c(T, d, N))
rvine_est=array(NA, dim=c(T, d, N))
tau=rvm_r$tau
gen_tau_diff=array(NA, dim = c(d, d, N))
lower_tau_diff=array(NA, dim = c(d, d, N))
upper_tau_diff=array(NA, dim = c(d, d, N))

#simulation
for (i in 1:N){
  rvine[,,i]=RVineSim(T,rvm_r)
}
#check missing values
sum(is.na(rvine))

#estimation
for (i in 1:N){
  rvm_est=RVineStructureSelect(rvine[,,i],progress=FALSE) #RVM = RVINE MATRIX
  
  #general tau difference
  rvine_est[,,i]=RVineSim(T,rvm_est) #RVINE EST=SAMPLE FROM RVM
  gen_tau_diff[,,i]=abs(TauMatrix(rvine[,,i])-TauMatrix(rvine_est[,,i]))
  
  #lower tau difference
  rvine_lower=matrix(nrow = 0, ncol = d)
  rvine_est_lower=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]<=0.2 & rvine[,k,i]<=0.2
      rvine_lower=rbind(rvine_lower,rvine[flg,,i])
      flg_est=rvine_est[,j,i]<=0.2 & rvine_est[,k,i]<=0.2
      rvine_est_lower=rbind(rvine_est_lower,rvine_est[flg_est,,i])
    }
    if (dim(rvine_lower)[1]>1 & dim(rvine_est_lower)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      lower_tau_diff[,,i]=abs(TauMatrix(rvine_lower)-TauMatrix(rvine_est_lower))
    }
  }
  
  #upper tau difference
  rvine_upper=matrix(nrow = 0, ncol = d)
  rvine_est_upper=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]>0.8 & rvine[,k,i]>0.8
      rvine_upper=rbind(rvine_upper,rvine[flg,,i])
      flg_est=rvine_est[,j,i]>0.8 & rvine_est[,k,i]>0.8
      rvine_est_upper=rbind(rvine_est_upper,rvine_est[flg_est,,i])
    }
    if (dim(rvine_upper)[1]>1 & dim(rvine_est_upper)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      upper_tau_diff[,,i]=abs(TauMatrix(rvine_upper)-TauMatrix(rvine_est_upper))
    }
  }
}

#display error
mean_tau_diff=apply(gen_tau_diff, c(1, 2), mean)
print(paste("mean generic tau difference 1",mean(mean_tau_diff[lower.tri(mean_tau_diff)])))
#T=500, N=10: 0.01361244393549
#T=500, N=100: 0.0111957247828991
#T=500, N=1000: 0.0115575364061456
#T=1000, N=1000: 0.00931771161637828

mean_lower_tau_diff=apply(lower_tau_diff, c(1, 2), mean)
print(paste("mean lower tau difference",mean(lower_tau_diff[lower.tri(lower_tau_diff)])))
#T=500, N=10: 0.0450569702002633
#T=500, N=100: 0.0393313219607121
#T=500, N=1000: 0.040562476569338
#T=1000, N=1000: 0.0329022128240668

mean_upper_tau_diff=apply(upper_tau_diff, c(1, 2), mean)
print(paste("mean upper tau difference",mean(upper_tau_diff[lower.tri(upper_tau_diff)])))
#T=500, N=10: 0.0721471873130148
#T=500, N=100: 0.0567947864616883
#T=500, N=1000: 0.0608024698527986
#T=1000, N=1000: 0.055513015028374















#SCENARIO WITH BB1 COPULA RANDOM ROTATION

#initializations
par1_mat=array(0,dim=c(d,d))
par2_mat=array(0,dim=c(d,d))
fam_mat=matrix(0,nrow=d,ncol=d)
random_fam=sample(c(7,17,27,37), 7, replace = TRUE)
fam_mat[lower.tri(fam_mat,diag=FALSE)]=random_fam
par1_mat=array(0,dim=c(d,d))
par2_mat=array(0,dim=c(d,d))
for (i in 1:d){
  for (j in 1:d){
    if (fam_mat[i,j]==7 || fam_mat[i,j]==17){
      par1_mat[i,j]=runif(1, 0, 7)
      par2_mat[i,j]=runif(1, 1, 7)
    }
    else if (fam_mat[i,j]==27 || fam_mat[i,j]==37){
      par1_mat[i,j]=runif(1, -7, 0)
      par2_mat[i,j]=runif(1, -7, -1)
    }
  }
}
rvm_r=RVineMatrix(Matrix=mat,family=fam_mat,par=par1_mat,par2=par2_mat)
rvine=array(NA, dim = c(T, d, N))
rvine_est=array(NA, dim=c(T, d, N))
tau=rvm_r$tau
gen_tau_diff=array(NA, dim = c(d, d, N))
lower_tau_diff=array(NA, dim = c(d, d, N))
upper_tau_diff=array(NA, dim = c(d, d, N))

#simulation
for (i in 1:N){
  rvine[,,i]=RVineSim(T,rvm_r)
}
#check missing values
sum(is.na(rvine))

#estimation
for (i in 1:N){
  rvm_est=RVineStructureSelect(rvine[,,i],progress=FALSE) #RVM = RVINE MATRIX
  
  #general tau difference
  rvine_est[,,i]=RVineSim(T,rvm_est) #RVINE EST=SAMPLE FROM RVM
  gen_tau_diff[,,i]=abs(TauMatrix(rvine[,,i])-TauMatrix(rvine_est[,,i]))
  
  #lower tau difference
  rvine_lower=matrix(nrow = 0, ncol = d)
  rvine_est_lower=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]<=0.2 & rvine[,k,i]<=0.2
      rvine_lower=rbind(rvine_lower,rvine[flg,,i])
      flg_est=rvine_est[,j,i]<=0.2 & rvine_est[,k,i]<=0.2
      rvine_est_lower=rbind(rvine_est_lower,rvine_est[flg_est,,i])
    }
    if (dim(rvine_lower)[1]>1 & dim(rvine_est_lower)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      lower_tau_diff[,,i]=abs(TauMatrix(rvine_lower)-TauMatrix(rvine_est_lower))
    }
  }
  
  #upper tau difference
  rvine_upper=matrix(nrow = 0, ncol = d)
  rvine_est_upper=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]>0.8 & rvine[,k,i]>0.8
      rvine_upper=rbind(rvine_upper,rvine[flg,,i])
      flg_est=rvine_est[,j,i]>0.8 & rvine_est[,k,i]>0.8
      rvine_est_upper=rbind(rvine_est_upper,rvine_est[flg_est,,i])
    }
    if (dim(rvine_upper)[1]>1 & dim(rvine_est_upper)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      upper_tau_diff[,,i]=abs(TauMatrix(rvine_upper)-TauMatrix(rvine_est_upper))
    }
  }
}

#display error
mean_tau_diff=apply(gen_tau_diff, c(1, 2), mean)
print(paste("mean generic tau difference 1",mean(mean_tau_diff[lower.tri(mean_tau_diff)])))
#T=500, N=10: 0.0104550434201737
#T=500, N=100: 0.0121030327321309
#T=500, N=1000: 0.01346073967774
#T=1000, N=1000: 0.0116576354653324

mean_lower_tau_diff=apply(lower_tau_diff, c(1, 2), mean)
print(paste("mean lower tau difference",mean(lower_tau_diff[lower.tri(lower_tau_diff)])))
#T=500, N=10: 0.0214083181761
#T=500, N=100: 0.020322939800837
#T=500, N=1000: 0.028608124109949
#T=1000, N=1000: 0.0301880199408099

mean_upper_tau_diff=apply(upper_tau_diff, c(1, 2), mean)
print(paste("mean upper tau difference",mean(upper_tau_diff[lower.tri(upper_tau_diff)])))
#T=500, N=10: 0.0575692439084721
#T=500, N=100: 0.0314271738330079
#T=500, N=1000: 0.0230840521430407
#T=1000, N=1000: 0.0500424886405396

















#SCENARIO WITH MIXED COPULAS MIXED TAUS TAKEN FROM CZADO

#initializations
par1_mat=array(0,dim=c(d,d))
par2_mat=array(0,dim=c(d,d))
fam_mat=matrix(0,nrow=d,ncol=d)
fam_mat=matrix(c(0,1,5,1,4,5,14, 0,0,1,5,14,1,4, 0,0,0,1,4,5,14,
                 0,0,0,0,14,1,4, 0,0,0,0,0,2,2, 0,0,0,0,0,0,2,
                 0,0,0,0,0,0,0),7,7)
tau=matrix(c(0,0.05,0.10,0.15,0.20,0.25,0.50, 0,0,0.10,0.15,0.20,0.30,0.55,
             0,0,0,0.15,0.20,0.35,0.60, 0,0,0,0,0.20,0.40,0.65,
             0,0,0,0,0,0.45,0.70, 0,0,0,0,0,0,0.70,
             0,0,0,0,0,0,0),7,7)
deg_free_t=3
for (i in 1:7){
  for (j in 1:7){
    if (i>j){
      par1_mat[i,j]=BiCopTau2Par(family=fam_mat[i,j],tau=tau[i,j])
      if (fam_mat[i,j]==2){
        par2_mat[i,j]=deg_free_t
        deg_free_t=deg_free_t+1
      }
    }
  }
}
rvm_r=RVineMatrix(Matrix=mat,family=fam_mat,par=par1_mat,par2=par2_mat)
rvine=array(NA, dim = c(T, d, N))
rvine_est=array(NA, dim=c(T, d, N))
gen_tau_diff=array(NA, dim = c(d, d, N))
lower_tau_diff=array(NA, dim = c(d, d, N))
upper_tau_diff=array(NA, dim = c(d, d, N))

#simulation and estimation
for (i in 1:N){
  rvine[,,i]=RVineSim(T,rvm_r)
  rvm_est=RVineStructureSelect(rvine[,,i],progress=FALSE) #RVM = RVINE MATRIX
  
  #general tau difference:
  #simulate data with true model, compute empirical tau from simulated data (A),
  #estimate the model, simulate data from the estimated model, 
  #compute empirical tau on simulated data from estimated model (B), 
  #compute A-B
  rvine_est[,,i]=RVineSim(T,rvm_est) #RVINE EST=SAMPLE FROM RVM
  gen_tau_diff[,,i]=abs(TauMatrix(rvine[,,i])-TauMatrix(rvine_est[,,i]))
  
  #lower tau difference
  rvine_lower=matrix(nrow = 0, ncol = d)
  rvine_est_lower=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]<=0.2 & rvine[,k,i]<=0.2
      rvine_lower=rbind(rvine_lower,rvine[flg,,i])
      flg_est=rvine_est[,j,i]<=0.2 & rvine_est[,k,i]<=0.2
      rvine_est_lower=rbind(rvine_est_lower,rvine_est[flg_est,,i])
    }
    if (dim(rvine_lower)[1]>1 & dim(rvine_est_lower)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      lower_tau_diff[,,i]=abs(TauMatrix(rvine_lower)-TauMatrix(rvine_est_lower))
    }
  }
  
  #upper tau difference
  rvine_upper=matrix(nrow = 0, ncol = d)
  rvine_est_upper=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]>0.8 & rvine[,k,i]>0.8
      rvine_upper=rbind(rvine_upper,rvine[flg,,i])
      flg_est=rvine_est[,j,i]>0.8 & rvine_est[,k,i]>0.8
      rvine_est_upper=rbind(rvine_est_upper,rvine_est[flg_est,,i])
    }
    if (dim(rvine_upper)[1]>1 & dim(rvine_est_upper)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      upper_tau_diff[,,i]=abs(TauMatrix(rvine_upper)-TauMatrix(rvine_est_upper))
    }
  }
}

#display error
mean_tau_diff=apply(gen_tau_diff, c(1, 2), mean)
print(paste("mean generic tau difference 1",mean(mean_tau_diff[lower.tri(mean_tau_diff)])))
#T=500, N=10: 0.0240433629163088
#T=500, N=100: 0.0217010478099055
#T=500, N=1000: 0.0202402946846073
#T=1000, N=1000: 0.0151380412793746

mean_lower_tau_diff=apply(lower_tau_diff, c(1, 2), mean)
print(paste("mean lower tau difference",mean(lower_tau_diff[lower.tri(lower_tau_diff)])))
#T=500, N=10: 0.0643451872350793
#T=500, N=100: 0.0574120432934252
#T=500, N=1000: 0.0584803049953999
#T=1000, N=1000: 0.0460839739232319

mean_upper_tau_diff=apply(upper_tau_diff, c(1, 2), mean)
print(paste("mean upper tau difference",mean(upper_tau_diff[lower.tri(upper_tau_diff)])))
#T=500, N=10: 0.0623787921797875
#T=500, N=100: 0.0632510552615448
#T=500, N=1000: 0.0668592147271203
#T=1000, N=1000: 0.0509377528813044

#study of anomalous cases
#desktop_path<- file.path(Sys.getenv("USERPROFILE"), "Desktop", "rvine_save.csv")
#rvine_save<- as.data.frame(as.table(rvine[,,c(292,514,533,736,745,989)]))
#write.csv(rvine_save, file =desktop_path, row.names = FALSE)
#read_df <- read.csv("rvine_save.csv")
#rvine_na<- array(read_df$Freq, dim = c(500, 7, 6))
#rvm_est of rvine[,,745] has got a NA in position (3,2) of tau matrix 
#(same for rvine[,,533]), this implies that mean_tau_diff2 has got also NA,
#in particular in position (2,1) (3,2) (4,2) (5,2) (6,1))


















#SCENARIO WITH EVERY POSSIBLE FAMILY RANDOMLY

#initializations
par1_mat=array(0,dim=c(d,d))
par2_mat=array(0,dim=c(d,d))
rvine=array(NA, dim = c(T, d, N))
rvine_est=array(NA, dim=c(T, d, N))
fam_mat=matrix(0,nrow=d,ncol=d)
gen_tau_diff=array(NA, dim = c(d, d, N))
lower_tau_diff=array(NA, dim = c(d, d, N))
upper_tau_diff=array(NA, dim = c(d, d, N))

#choose copula family and parameters
copula_families=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 
                  14, 16, 17, 18, 19, 20, 23, 24, 26, 27, 
                  28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 
                  104, 114, 124, 134, 204, 214, 224, 234)
random_fam=sample(copula_families, dd, replace=TRUE)
fam_mat[lower.tri(fam_mat,diag=FALSE)]=random_fam
for (i in 1:d) {
  for (j in 1:d) {
    if (fam_mat[i, j] %in% c(1, 2)) {
      par1_mat[i, j] = runif(1, -1, 1)
      if (fam_mat[i, j] == 2) par2_mat[i, j] = runif(1, 2, 1e3)
    } else if (fam_mat[i, j] %in% c(3, 13)) {
      par1_mat[i, j] = runif(1, 0, 28)
    } else if (fam_mat[i, j] %in% c(4, 14)) {
      par1_mat[i, j] = runif(1, 1, 17)
    } else if (fam_mat[i, j] == 5) {
      prova=runif(1, -35, 35)
      while (prova==0){
        prova=runif(1, -35, 35)
      }
      par1_mat[i, j] = prova
    } else if (fam_mat[i, j] %in% c(6, 16)) {
      par1_mat[i, j] = runif(1, 1, 30)
    } else if (fam_mat[i, j] %in% c(7, 17)) {
      par1_mat[i, j] = runif(1, 0, 7)
      par2_mat[i, j] = runif(1, 1, 7)
    } else if (fam_mat[i, j] %in% c(8, 18)) {
      par1_mat[i, j] = runif(1, 1, 6)
      par2_mat[i, j] = runif(1, 1, 8)
    } else if (fam_mat[i, j] %in% c(9, 19)) {
      par1_mat[i, j] = runif(1, 1, 6)
      par2_mat[i, j] = runif(1, 0, 75)
    } else if (fam_mat[i, j] %in% c(10, 20)) {
      par1_mat[i, j] = runif(1, 1, 8)
      par2_mat[i, j] = runif(1, 1e-4, 1)
    } else if (fam_mat[i, j] %in% c(23, 33)) {
      par1_mat[i, j] = runif(1, -28, 0)
    } else if (fam_mat[i, j] %in% c(24, 34)) {
      par1_mat[i, j] = runif(1, -17, -1)
    } else if (fam_mat[i, j] %in% c(26, 36)) {
      par1_mat[i, j] = runif(1, -30, -1)
    } else if (fam_mat[i, j] %in% c(27, 37)) {
      par1_mat[i, j] = runif(1, -7, 0)
      par2_mat[i, j] = runif(1, -7, -1)
    } else if (fam_mat[i, j] %in% c(28, 38)) {
      par1_mat[i, j] = runif(1, -6, -1)
      par2_mat[i, j] = runif(1, -8, -1)
    } else if (fam_mat[i, j] %in% c(29, 39)) {
      par1_mat[i, j] = runif(1, -6, -1)
      par2_mat[i, j] = runif(1, -75, 0)
    } else if (fam_mat[i, j] %in% c(30, 40)) {
      par1_mat[i, j] = runif(1, -8, -1)
      par2_mat[i, j] = runif(1, -1, -1e-4)
    } else if (fam_mat[i, j] %in% c(41, 51)) {
      par1_mat[i, j] = runif(1, 0, 1e10)
    } else if (fam_mat[i, j] %in% c(61, 71)) {
      par1_mat[i, j] = runif(1, -1e3, 0)
    } else if (fam_mat[i, j] == 42) {
      b <- runif(1, -1, 1)
      limA <- (b - 3 - sqrt(9 + 6 * b - 3 * b^2))/2
      while (limA>1) {
        b <- runif(1, -1, 1)
        limA <- (b - 3 - sqrt(9 + 6 * b - 3 * b^2))/2
      }
      par1_mat[i, j] = runif(1,limA,1)
      par2_mat[i, j] = b
    } else if (fam_mat[i, j] %in% c(104, 114, 204, 214)) {
      par1_mat[i, j] = runif(1, 1, 1e3)
      par2_mat[i, j] = runif(1, 0, 1)
    } else if (fam_mat[i, j] %in% c(124, 134, 224, 234)) {
      par1_mat[i, j] = runif(1, -1e3, -1)
      par2_mat[i, j] = runif(1, 0, 1)
    }
  }
}
rvm_r=RVineMatrix(Matrix=mat,family=fam_mat,par=par1_mat,par2=par2_mat)
tau=rvm_r$tau 

#simulation
for (i in 1:N){
  rvine[,,i]=RVineSim(T,rvm_r)
}
#check missing values
sum(is.na(rvine))

#estimation
for (i in 1:N){
  rvm_est=RVineStructureSelect(rvine[,,i],progress=FALSE) #RVM = RVINE MATRIX
  
  #general tau difference
  rvine_est[,,i]=RVineSim(T,rvm_est) #RVINE_EST = SAMPLE FROM RVM
  gen_tau_diff[,,i]=abs(TauMatrix(rvine[,,i])-TauMatrix(rvine_est[,,i]))
  
  #lower tau difference
  rvine_lower=matrix(nrow = 0, ncol = d)
  rvine_est_lower=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]<=0.2 & rvine[,k,i]<=0.2
      rvine_lower=rbind(rvine_lower,rvine[flg,,i])
      flg_est=rvine_est[,j,i]<=0.2 & rvine_est[,k,i]<=0.2
      rvine_est_lower=rbind(rvine_est_lower,rvine_est[flg_est,,i])
    }
    if (dim(rvine_lower)[1]>1 & dim(rvine_est_lower)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      lower_tau_diff[,,i]=abs(TauMatrix(rvine_lower)-TauMatrix(rvine_est_lower))
    }
  }
  
  #upper tau difference
  rvine_upper=matrix(nrow = 0, ncol = d)
  rvine_est_upper=matrix(nrow = 0, ncol = d)
  for (j in 1:(d-1)){
    for (k in (j+1):d){
      flg=rvine[,j,i]>0.8 & rvine[,k,i]>0.8
      rvine_upper=rbind(rvine_upper,rvine[flg,,i])
      flg_est=rvine_est[,j,i]>0.8 & rvine_est[,k,i]>0.8
      rvine_est_upper=rbind(rvine_est_upper,rvine_est[flg_est,,i])
    }
    if (dim(rvine_upper)[1]>1 & dim(rvine_est_upper)[1]>1){
      rvine_lower=unique(rvine_lower)
      rvine_est_lower=unique(rvine_est_lower)
      upper_tau_diff[,,i]=abs(TauMatrix(rvine_upper)-TauMatrix(rvine_est_upper))
    }
  }
}

#display error
mean_tau_diff=apply(gen_tau_diff, c(1, 2), mean)
print(paste("mean generic tau difference",mean(mean_tau_diff[lower.tri(mean_tau_diff)])))
#T=500, N=10: 0.0238954098673538
#T=500, N=100: 0.0806560380616617
#T=500, N=1000: 0.0357822972386742
#T=1000, N=1000: 0.0174836106895621

mean_lower_tau_diff=apply(lower_tau_diff, c(1, 2), mean)
print(paste("mean lower tau difference",mean(lower_tau_diff[lower.tri(lower_tau_diff)])))
#T=500, N=10: 0.0341544492487925
#T=500, N=100: 0.0756173075316744
#T=500, N=1000: 0.0433710712192429
#T=1000, N=1000: 0.0297517949215619

mean_upper_tau_diff=apply(upper_tau_diff, c(1, 2), mean)
print(paste("mean upper tau difference",mean(upper_tau_diff[lower.tri(upper_tau_diff)])))
#T=500, N=10: 0.0413958695116491
#T=500, N=100: 0.0795217912376646
#T=500, N=1000: 0.0571722129818986
#T=1000, N=1000: 0.115660439599446
