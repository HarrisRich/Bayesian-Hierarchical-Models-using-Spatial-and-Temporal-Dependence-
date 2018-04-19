
library(MASS)
#install.packages("forcats")
library(forcats)


data3=read.csv("neighbor_matrix16.csv",header = 1)
stable3=read.csv("WorkingDataSetToRun16.csv",header = 1)
table3=stable3
B3=data3[,2:2898]
B3=B3[-which(table3$siteID=="P08"),-which(table3$siteID=="P08")]
table3=table3[- which(table3$siteID=="P08"),]
lion.present=vector(length = length(table3$lionFemale.present))
for(i in 1:length(table3$lionFemale.present)){
  lion.present[i]=max(table3$lionFemale.present[i],table3$lionMale.present[i])
}

table3=cbind(table3,lion.present)
A3=B3


n=dim(table3)[1]

TID16=seq(1,length(unique(table3$date)),1)
datesort16 = sort(unique(table3$date))

CID16 = seq(1,length(unique(table3$siteID)),1)
camsort16 = sort(unique(table3$siteID))

C16 = rep(0,n)
for(i in 1:n){
  Ind0 = which(camsort16==table3$siteID[i])
  C16[i]=CID16[Ind0]
}

TP16 = rep(0,n)
for(i in 1:n){
  Ind0 = which(datesort16==table3$date[i])
  TP16[i]=TID16[Ind0]
}

At16 = matrix(0,n,n) 
for(i in 1:n){
  Ind1 = which(C16==C16[i] & TP16==(TP16[i]-1))
  if(length(Ind1)>0){
    At16[i,Ind1] = 1
  }
}

b3=as.numeric(table3$date)
c=(1:(length(unique(table3$date))))
for(i in 1:length(unique(table3$date))){
  b3[b3==(sort(unique(table3$date[i])))]=c[i]
}
#we will start with the psi hats and mu
tree3=vector(length = length(table3$T50))
for(i in 1:length(table3$T50)){
  if(table3$T50[i]>-1 & table3$T50[i]<11){tree3[i]=0}  
  if(table3$T50[i]>10 & table3$T50[i]<21){tree3[i]=1}
  if(table3$T50[i]>20){tree3[i]=2}  
}

X16=matrix(c(rep(1,dim(A3)[1]),table3$ndvi,(table3$ndvi)^2,sqrt(table3$amRivDist),
             table3$wet,table3$TM50,table3$fire,tree3),nrow = dim(A3)[1])
#X16=matrix(c(rep(1,dim(A3)[1]),table3$ndvi,(table3$ndvi)^2,sqrt(table3$amRivDist),
#           table3$TM50,tree3),nrow = dim(A3)[1])
A3=unname(as.matrix(A3))
LLst16 = function(theta0) {
  eta=theta0[1]
  delta=theta0[2]
  beta=theta0[3:length(theta0)]
  mu = (exp(X16%*%beta))/(1+exp(X16%*%beta))
  Ind1 = which(TP16==1)
  Ind2 = which(TP16>1)
  -1*(t(Z)%*%(X16%*%beta+eta*A3%*%(Z-mu)+delta*At16%*%(Z-mu)) -
        t(rep(1,sum(TP16==1)))%*%log(1+exp(X16[Ind1,]%*%beta+eta*A3[Ind1,]%*%(Z-mu))) - 
        t(rep(1,sum(TP16>1)))%*%log(1+exp(X16[Ind2,]%*%beta+eta*A3[Ind2,]%*%(Z-mu)+delta*At16[Ind2,]%*%(Z-mu))))
}
psi16 = function(fit,Z) {
  eta=fit$par[1]
  delta=fit$par[2]
  beta=fit$par[3:length(fit$par)]
  mu = (exp(X16%*%beta))/(1+exp(X16%*%beta))
  Ind1 = which(TP16==1)
  Ind2 = which(TP16>1)
  (X16%*%beta+eta*A3%*%(Z-mu)+delta*At16%*%(Z-mu)) -c((log(1+exp(X16[Ind1,]%*%beta+eta*A3[Ind1,]%*%(Z-mu)))),
                                                      (log(1+exp(X16[Ind2,]%*%beta+eta*A3[Ind2,]%*%(Z-mu)+delta*At16[Ind2,]%*%(Z-mu)))))
}

psihat16Z=matrix(ncol = 3, nrow =  (length = dim(A3)[1]))
Z=table3$zebra.present
fit = optim(c(1,1,1,1,1,1,1,1,1,1),LLst16, method = "L-BFGS-B", hessian = TRUE) 
betahatZ16 = vector(length = dim(X16)[2])
betahatZ16 = fit$par[3:10]

psihat16Z[,1]=psi16(fit,table3$zebra.present)
psihat16Z[,2]=b3
psihat16Z[,3]=table3$siteID

save(psihat16Z,betahatZ16,file = "Beta16AndPsi16ForZ.Rda")