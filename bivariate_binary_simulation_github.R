##################################################################################################################################################################################################################################################
##### Simulation code for the joint model with principal stratification model and outcome regression model with bivariate binary outcomes and two study covariates in cluster randomized trials subject to missingness#####
##################################################################################################################################################################################################################################################
require("LaplacesDemon")
require("mvtnorm")
require("pbv")
require("dplyr")
require("truncnorm")
### function to sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{ 
  
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

###################################
## 1. Data generating process  ####
###################################

# N is the total number of individuals
# s is the number of clusters. 
gendata<-function(N=1800, s=60,cv=0.3){
 
  m=N/s   #mean cluster size
  tau2<-matrix(c(0.1,0.5*sqrt(0.01),0.5*sqrt(0.01),0.1),ncol=2,nrow=2)#cluster variance-covariance matrix
  
  ##cluster ids and two covariates
  a0=1/cv^2
  b0=1/m/cv^2
  cl_size=round(rgamma(n=s,shape=a0,rate=b0))
  cl_size[which(cl_size<1)]<-1 #set the lower bound to 1
  cl<-rep(1:s,cl_size) #cluster index
  size=rep(cl_size,cl_size)
  x1<-rnorm(N,0,1) #continuous 
  x2<-runif(N,0,3)  #uniform
  X<-as.matrix(data.frame(rep(1,N),x1,x2,size))   #covariate matrix
  N=length(cl)
  
  ##strata model parameters
  beta=c(-0.3,0.5,-0.7) #first Probit layer effect (strata proportion with 0.15 never survivors, 0.20 protected and 0.65 always survivors)
  gamma=c(-1.5,-0.6,0.4) #second Probit layer effect (strata proportion with 0.15 never survivors, 0.20 protected and 0.65 always survivors)
  #beta=c(-0.1,0.5,-0.7) #first Probit layer effect (strata proportion with 0.20 never survivors, 0.25 protected and 0.55 always survivors)
  #gamma=c(-1.3,-0.6,0.4) #second Probit layer effect (strata proportion with 0.20 never survivors, 0.25 protected and 0.55 always survivors)
  
  phi2=0.1#cluster variance
  chi<-rnorm(s,0,sqrt(phi2))  #cluster effect
  chin<-rep(chi,cl_size) #cluster effect at individual level
  
  #generate latent variable for G group status
  #principal strata membership (0: always non-survivor,1: protected, 2: always survivor)
  p<-pnorm(X[,1:3]%*%beta+chin) 
  G00 <- rbinom(N, 1, p)
  G01<-rep(0,N)
  tgt<-which(G00==0)
  p<-pnorm(X[tgt,1:3]%*%gamma+chin[tgt])
  G01[tgt]<-rbinom(length(tgt), 1, p)
  G11<-rep(1,N)-G01-G00
  G=1*G01+2*G11
  
  #outcome model parameters
  d<-rep(sample(rep(c(0,1),each=s/2),replace=F),cl_size)
  
  #outcome model parameters
  alpha111a<-c(-3.4,-0.1,-0.1,0.1)
  alpha101a<-c(-2.3,-0.3,-0.6,0.1)
  alpha110a<-c(-3.4,0.1,0.1,0.1)
  
  alpha111b<-c(-3.3,-0.2,-0.1,0.1)
  alpha101b<-c(-2.4,-0.4,-0.5,0.1)
  alpha110b<-c(-3.3,0.2,0.1,0.1)
  eta<-rmvnorm(s,c(0,0),tau2) #cluster effect
  etan<-cbind(rep(eta[,1],cl_size),rep(eta[,2],cl_size)) #cluster effect at individual level
  
  #binary outcomes
  # truncation by death
  Y1<-rep(NA,N) #outcome (NA are death)
  Y1[intersect(which(G11==1),which(d==1))]<-rbinom(length(intersect(which(G11==1),which(d==1))), 1, pnorm(X[intersect(which(G11==1),which(d==1)),]%*%alpha111a+etan[intersect(which(G11==1),which(d==1)),1]))
  Y1[intersect(which(G01==1),which(d==1))]<-rbinom(length(intersect(which(G01==1),which(d==1))), 1, pnorm(X[intersect(which(G01==1),which(d==1)),]%*%alpha101a+etan[intersect(which(G01==1),which(d==1)),1]))
  Y1[intersect(which(G11==1),which(d==0))]<-rbinom(length(intersect(which(G11==1),which(d==0))), 1, pnorm(X[intersect(which(G11==1),which(d==0)),]%*%alpha110a+etan[intersect(which(G11==1),which(d==0)),1]))
  
  Y2<-rep(NA,N) #outcome (NA are death)
  Y2[intersect(which(G11==1),which(d==1))]<-rbinom(length(intersect(which(G11==1),which(d==1))), 1, pnorm(X[intersect(which(G11==1),which(d==1)),]%*%alpha111b+etan[intersect(which(G11==1),which(d==1)),2]))
  Y2[intersect(which(G01==1),which(d==1))]<-rbinom(length(intersect(which(G01==1),which(d==1))), 1, pnorm(X[intersect(which(G01==1),which(d==1)),]%*%alpha101b+etan[intersect(which(G01==1),which(d==1)),2]))
  Y2[intersect(which(G11==1),which(d==0))]<-rbinom(length(intersect(which(G11==1),which(d==0))), 1, pnorm(X[intersect(which(G11==1),which(d==0)),]%*%alpha110b+etan[intersect(which(G11==1),which(d==0)),2]))
  
  s<-ifelse(is.na(Y1),0,1) # survival status
  print(table(Y1,G))
  print(table(Y2,G))
  #true OR and RR
  P1a<-pnorm(X[which(G==2),1:4]%*%alpha111a)
  P0a<-pnorm(X[which(G==2),1:4]%*%alpha110a)
  ORa_i<-(mean(P1a)/(1-mean(P1a)))/(mean(P0a)/(1-mean(P0a))) #outcome 1 SIACE
  RRa_i<-mean(P1a)/mean(P0a) #outcome 1 SIACE
  
  P1b<-pnorm(X[which(G==2),1:4]%*%alpha111b)
  P0b<-pnorm(X[which(G==2),1:4]%*%alpha110b)
  ORb_i<-(mean(P1b)/(1-mean(P1b)))/(mean(P0b)/(1-mean(P0b)))  #outcome 2 SIACE
  RRb_i<-mean(P1b)/mean(P0b) #outcome 2 SIACE
  
  X1<-data.frame(cbind(P1a,P0a,P1b,P0b,cl[which(G==2)]))
  ORa_c<-mean(by((X1$X1/(1-X1$X1))/(X1$X2/(1-X1$X2)),X1$X5, mean))  #outcome 1 SCACE
  RRa_c<-mean(by(X1$X1/X1$X2,X1$X5, mean)) #outcome 1 SCACE
  ORb_c<-mean(by((X1$X3/(1-X1$X3))/(X1$X4/(1-X1$X4)),X1$X5, mean)) #outcome 2 SCACE
  RRb_c<-mean(by(X1$X3/X1$X4,X1$X5, mean)) #outcome 2 SCACE
  
  #missingness model
  ##missing both survival status and outcomes
  z1=c(1.3,0.5,-0.1)
  z3=c(2.1,0.5,-0.3)
  pr1=pnorm(X[,1:3]%*%z1)
  k<-rbinom(length((s)),1,pr1)
  s<-ifelse(k==1,s,NA)
  Y1[which(is.na(s))]<-NA
  Y2[which(is.na(s))]<-NA

  # only missing outcome
  mb<-which(is.na(Y1)&is.na(s))
  pr3=pnorm(X[-mb,1:3]%*%z3)
  k2<-rbinom((length(s)-length(mb)),1,pr3)
  Y1[-mb]<-ifelse(k2==1,Y1[-mb],NA)
  Y2[-mb]<-ifelse(k2==1,Y2[-mb],NA)
  #For G, the true status of 0, 1 or 2 will not be known, and death in treatment group and survival in control group have known G status   
  return(data.frame(cbind(Y1,Y2,x1,x2,d,cl,s,size,G,ORa_i,RRa_i,ORb_i,RRb_i,ORa_c,RRa_c,ORb_c,RRb_c)))
}


###################################
## 2. MCMC model estimation  ######
###################################
#df: the simulated data 
#S: length of chain

mvsampler<-function(df,S=5000,dau1=0.005){
  # extract data information
  N=nrow(df) # number of individuals
  cl=df$cl # cluster ID
  Nst<-length(which(!is.na(df$Y)))
  s<-length(unique(df$cl)) # number of clusters
  cs<-as.numeric(table(cl)) # size of cluster
  mbar=N/s
  p=3 #number of covariates(stra model)
  X<-as.matrix(cbind(rep(1,N),df$x1,df$x2,df$size)) #covariate matrix
  surv=df$s # survival status
  Y1=df$Y1 # outcome 1
  Y2=df$Y2 # outcome 2
  D<-df$d #treatment assignment
  
  
  ##starting values
  #outcome model
  
  alpha111a<-c(-3.4,-0.1,-0.1,0.1)
  alpha101a<-c(-2.3,-0.3,-0.6,0.1)
  alpha110a<-c(-3.4,0.1,0.1,0.1)
  
  alpha111b<-c(-3.3,-0.2,-0.1,0.1)
  alpha101b<-c(-2.4,-0.4,-0.5,0.1)
  alpha110b<-c(-3.3,0.2,0.1,0.1)

  #strata model
  beta=c(-0.3,0.5,-0.7) #first Probit layer effect (strata proportion with 0.15 never survivors, 0.20 protected and 0.65 always survivors)
  gamma=c(-1.5,-0.6,0.4) #second Probit layer effect (strata proportion with 0.15 never survivors, 0.20 protected and 0.65 always survivors)
  #beta=c(-0.1,0.5,-0.7) #first Probit layer effect (strata proportion with 0.20 never survivors, 0.25 protected and 0.55 always survivors)
  #gamma=c(-1.3,-0.6,0.4) #second Probit layer effect (strata proportion with 0.20 never survivors, 0.25 protected and 0.55 always survivors)
  
  Tau=matrix(c(0.1,0.5*sqrt(0.01),0.5*sqrt(0.01),0.1),ncol=2,nrow=2)
  Sigma<-diag(2)
  # random effect at individual level
  eta<-rmvnorm(s,c(0,0),Tau)
  etan1<-rep(eta[,1],cs)
  etan2<-rep(eta[,2],cs)
  
  #initial value for G 
  G=df$G
  
  ##weakly informative priors
  #outcome model regression coefficient
  Sigma111<-Sigma101<-Sigma110<-diag(2*(p+1))*1000
  a=0.001  #prior variance terms
  b=0.001
  cc=0.001
  dd=0.001
  
  #strata model
  #regression coefficients
  beta0=gamma0=rep(0,p)
  Lambda=Gamma=diag(p)*1000
  
  g=0.001  #prior variance terms
  h=0.001
  # random effects
  phi=0.1
  chi1<-rnorm(s,0,sqrt(phi))
  chin1<-rep(chi1,cs)
  
  
  
  tgts1=which(is.na(Y1)&surv==1&D==1) # survive but missing outcomes under treated
  tgts2=which(is.na(Y1)&surv==1&D==0) # survive but missing outcomes under control
  tgtss=which(is.na(Y1)&is.na(surv)) # missing both outcomes and survival status
  miss=c(which(is.na(Y1)&is.na(surv)),which(is.na(Y1)&surv==1)) # all individual with misingness
  nonmiss=c(1:N)[-miss] # individuals without missingness
  cs1<-as.numeric(table(cl[nonmiss]))
  Nst<-length(intersect(which(!is.na(Y1)),nonmiss))
  
  #missing value imputation for individuals who survive with missing outcomes
  Y1[intersect(tgts1,which(G==2))]<-rbinom(length(intersect(tgts1,which(G==2))),1, pnorm(X[intersect(tgts1,which(G==2)),]%*%alpha111a+etan1[intersect(tgts1,which(G==2))]))
  Y1[intersect(tgts2,which(G==2))]<-rbinom(length(intersect(tgts2,which(G==2))),1, pnorm(X[intersect(tgts2,which(G==2)),]%*%alpha110a+etan1[intersect(tgts2,which(G==2))]))
  Y1[intersect(tgts1,which(G==1))]<-rbinom(length(intersect(tgts1,which(G==1))),1, pnorm(X[intersect(tgts1,which(G==1)),]%*%alpha101a+etan1[intersect(tgts1,which(G==1))]))
  
  Y2[intersect(tgts1,which(G==2))]<-rbinom(length(intersect(tgts1,which(G==2))),1, pnorm(X[intersect(tgts1,which(G==2)),]%*%alpha111b+etan2[intersect(tgts1,which(G==2))]))
  Y2[intersect(tgts2,which(G==2))]<-rbinom(length(intersect(tgts2,which(G==2))),1, pnorm(X[intersect(tgts2,which(G==2)),]%*%alpha110b+etan2[intersect(tgts2,which(G==2))]))
  Y2[intersect(tgts1,which(G==1))]<-rbinom(length(intersect(tgts1,which(G==1))),1, pnorm(X[intersect(tgts1,which(G==1)),]%*%alpha101b+etan2[intersect(tgts1,which(G==1))]))
  
  ## impute S and Y for both missing
  surv[tgtss]<-ifelse(G[tgtss]==2,1,ifelse(G[tgtss]==1&D[tgtss]==1,1,0))
  tgtss1=which(surv[tgtss]==1&D[tgtss]==1)
  tgtss2=which(surv[tgtss]==1&D[tgtss]==0)
  
  Y1[intersect(tgtss[tgtss1],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111a+etan1[intersect(tgtss[tgtss1],which(G==2))]))
  Y1[intersect(tgtss[tgtss2],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss2],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110a+etan1[intersect(tgtss[tgtss2],which(G==2))]))
  Y1[intersect(tgtss[tgtss1],which(G==1))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==1))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101a+etan1[intersect(tgtss[tgtss1],which(G==1))]))
  
  Y2[intersect(tgtss[tgtss1],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111b+etan2[intersect(tgtss[tgtss1],which(G==2))]))
  Y2[intersect(tgtss[tgtss2],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss2],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110b+etan2[intersect(tgtss[tgtss2],which(G==2))]))
  Y2[intersect(tgtss[tgtss1],which(G==1))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==1))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101b+etan2[intersect(tgtss[tgtss1],which(G==1))]))
  
  #set Z/W using the known strata
  Z<-rep(NA,N)
  W<-rep(NA,N)
  Z[G==0]<-rtruncnorm(length(which(G==0)),a=0,mean=(X[G==0,1:3]%*%c(beta)+chin1[G==0]),sd=1)
  Z[G>0]<-rtruncnorm(length(which(G>0)),b=0,mean=(X[G>0,1:3]%*%c(beta)+chin1[G>0]),sd=1)
  W[G==1]<-rtruncnorm(length(which(G==1)),a=0,mean=(X[G==1,1:3]%*%c(gamma)+chin1[G==1]),sd=1)
  W[G==2]<-rtruncnorm(length(which(G==2)),b=0,mean=(X[G==2,1:3]%*%c(gamma)+chin1[G==2]),sd=1)
  W[G==0]<-NA #set back to 0 is very crucial
  # 
  
  #set U111, U110 U101 using information
  U111a<-U110a<-U101a<-rep(NA,N)
  U111b<-U110b<-U101b<-rep(NA,N)
  
  tgt<-intersect(which(Y1==1), intersect(which(G==2),which(D==1)))
  U111a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha111a)+etan1[tgt],sd=1)
  tgt<-intersect(which(Y1==0), intersect(which(G==2),which(D==1)))
  U111a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha111a)+etan1[tgt],sd=1)
  tgt<-intersect(which(Y1==1), intersect(which(G==1),which(D==1)))
  U101a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha101a)+etan1[tgt],sd=1)
  tgt<-intersect(which(Y1==0), intersect(which(G==1),which(D==1)))
  U101a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha101a)+etan1[tgt],sd=1)
  tgt<-intersect(which(Y1==1), intersect(which(G==2),which(D==0)))
  U110a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha110a)+etan1[tgt],sd=1)
  tgt<-intersect(which(Y1==0), intersect(which(G==2),which(D==0)))
  U110a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha110a)+etan1[tgt],sd=1)
  
  tgt<-intersect(which(Y2==1), intersect(which(G==2),which(D==1)))
  U111b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha111b)+etan2[tgt],sd=1)
  tgt<-intersect(which(Y2==0), intersect(which(G==2),which(D==1)))
  U111b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha111b)+etan2[tgt],sd=1)
  tgt<-intersect(which(Y2==1), intersect(which(G==1),which(D==1)))
  U101b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha101b)+etan2[tgt],sd=1)
  tgt<-intersect(which(Y2==0), intersect(which(G==1),which(D==1)))
  U101b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha101b)+etan2[tgt],sd=1)
  tgt<-intersect(which(Y2==1), intersect(which(G==2),which(D==0)))
  U110b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha110b)+etan2[tgt],sd=1)
  tgt<-intersect(which(Y2==0), intersect(which(G==2),which(D==0)))
  U110b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha110b)+etan2[tgt],sd=1)
  
  
  ##store chains for parameters/estimates with length S
  ALPHA111a<-ALPHA101a<-ALPHA110a<-matrix(NA,p+1,S)
  ALPHA111b<-ALPHA101b<-ALPHA110b<-matrix(NA,p+1,S)
  TAU<-matrix(NA,4,S)
  OR<-RR<-matrix(NA,4,S) #sample TE by estimating the potential outcomes of G=2.
  
  BETA<-GAMMA<-matrix(NA,p,S)
  PHI<-rep(0,S)
  GV<-matrix(NA,3,S)
  
  #store G trajectory
  GT<-matrix(NA, nrow(df), S)
  rownames(GT)<-df$GP
  
  for(i in 1:S){
    
    #update alpha
    alpha111<-c(alpha111a,alpha111b)
    tgt<-intersect(which(D==1),which(G==2)) 
    InvSigma<-solve(Sigma)
    a=InvSigma[1,1]
    b=c=InvSigma[1,2]
    d=InvSigma[2,2]
    
    res111<-cbind(U111a[tgt],U111b[tgt])-cbind(etan1[tgt],etan2[tgt])
    ssqx111<-crossprod(X[tgt,])
    V111<-solve(kronecker(InvSigma,ssqx111)+solve(Sigma111))
    ssqxy111<-crossprod(X[tgt,],res111)
    M111<-V111%*%(rbind(a*ssqxy111[,1,drop=F]+b*ssqxy111[,2,drop=F],
                        c*ssqxy111[,1,drop=F]+d*ssqxy111[,2,drop=F]))
    alpha111<-as.numeric(rmvnorm(1,M111,V111))
    alpha111a<-alpha111[1:4]
    alpha111b<-alpha111[5:8]
    ALPHA111a[,i]<-alpha111a
    ALPHA111b[,i]<-alpha111b
    
    #update alpha101
    alpha101<-c(alpha101a,alpha101b)
    tgt<-intersect(which(D==1),which(G==1)) 
    
    res101<-cbind(U101a[tgt],U101b[tgt])-cbind(etan1[tgt],etan2[tgt])
    ssqx101<-crossprod(X[tgt,])
    V101<-solve(kronecker(InvSigma,ssqx101)+solve(Sigma101))
    ssqxy101<-crossprod(X[tgt,],res101)
    M101<-V101%*%(rbind(a*ssqxy101[,1,drop=F]+b*ssqxy101[,2,drop=F],
                        c*ssqxy101[,1,drop=F]+d*ssqxy101[,2,drop=F]))
    alpha101<-as.numeric(rmvnorm(1,M101,V101))
    alpha101a<-alpha101[1:4]
    alpha101b<-alpha101[5:8]
    ALPHA101a[,i]<-alpha101a
    ALPHA101b[,i]<-alpha101b
    
    #update alpha110
    alpha110<-c(alpha110a,alpha110b)
    tgt<-intersect(which(D==0),which(G==2)) 
    res110<-cbind(U110a[tgt],U110b[tgt])-cbind(etan1[tgt],etan2[tgt])
    ssqx110<-crossprod(X[tgt,])
    V110<-solve(kronecker(InvSigma,ssqx110)+solve(Sigma110))
    ssqxy110<-crossprod(X[tgt,],res110)
    M110<-V110%*%(rbind(a*ssqxy110[,1,drop=F]+b*ssqxy110[,2,drop=F],
                        c*ssqxy110[,1,drop=F]+d*ssqxy110[,2,drop=F]))
    alpha110<-as.numeric(rmvnorm(1,M110,V110))
    alpha110a<-alpha110[1:4]
    alpha110b<-alpha110[5:8]
    ALPHA110a[,i]<-alpha110a
    ALPHA110b[,i]<-alpha110b
    
    #update eta part
    for(j in 1:s){
      if (D[which(cl==j)[1]]==0){
        tgt3<-intersect(which(G==2),intersect(which(D==0),which(df$cl==j)))
        V_eta=solve(length(tgt3)*InvSigma + solve(Tau))
        res_eta=matrix(cbind(U110a[tgt3],U110b[tgt3])-cbind(X[tgt3,]%*%alpha110a, X[tgt3,]%*%alpha110b),length(tgt3),2)
        M_eta=V_eta%*%(InvSigma %*% colSums(res_eta))
        eta[j,]<-as.numeric(rmvnorm(1,M_eta,V_eta))
        
      }else{
        tgt1<-intersect(which(G==2),intersect(which(D==1),which(df$cl==j)))
        tgt2<-intersect(which(G==1),intersect(which(D==1),which(df$cl==j)))
        
        V_eta=solve((length(tgt1)+length(tgt2))*InvSigma + solve(Tau))
        res_eta1=matrix(cbind(U111a[tgt1],U111b[tgt1])-cbind(X[tgt1,]%*%alpha111a, X[tgt1,]%*%alpha111b),length(tgt1),2)
        res_eta2=matrix(cbind(U101a[tgt2],U101b[tgt2])-cbind(X[tgt2,]%*%alpha101a, X[tgt2,]%*%alpha101b),length(tgt2),2)
        M_eta=V_eta%*%(InvSigma %*% (colSums(res_eta1) + colSums(res_eta2)))
        eta[j,]<-as.numeric(rmvnorm(1,M_eta,V_eta))
      }
    }
    etan1<-rep(eta[,1],cs)
    etan2<-rep(eta[,2],cs)
    
    #update tau2
    Tau<-solve(rwish(1,2+s,solve(matrix(c(0.1,0.05,0.05,0.1),2,2)+crossprod(eta))))
    TAU[,i]<-as.numeric(Tau)
    
    #update beta/gamma
    v=solve(crossprod(X[,1:3])+solve(Lambda))
    m=v%*%(crossprod(X[,1:3],Z-chin1)+solve(Lambda)%*%beta0)
    beta=rmvnorm(1,m,v)
    BETA[,i]<-beta
    
    tgt<-which(G>0)
    v=solve(crossprod(X[tgt,1:3])+solve(Gamma))
    m=v%*%(crossprod(X[tgt,1:3],W[tgt]-chin1[tgt])+solve(Gamma)%*%gamma0)
    gamma=rmvnorm(1,m,v)
    GAMMA[,i]<-gamma
    
    #update phi
    phi<-1/rgamma(1,g+s/2,h+0.5*sum(chi1^2))
    PHI[i]<-phi
    
    #update chi
    v=apply(as.data.frame(cs+1/phi),1,solve)
    ddd=as.data.frame(cbind(Z-(X[,1:3])%*%c(beta),cl))
    m=v*tapply(ddd$V1,ddd$cl,sum)
    
    chi1<-rnorm(s,m,sqrt(v))
    chin1<-rep(chi1,cs)
    
    #update membership
    tgt<-intersect(which(is.na(Y1)),which(D==0))
    A=pnorm(X[tgt,1:3]%*%t(beta)+chin1[tgt]) #probability always death/(always death+protected+always survival)
    B=pnorm(X[tgt,1:3]%*%t(gamma)+chin1[tgt]) #probability protected/(protected+always survival)
    PA<-(1-A)*B/(A+(1-A)*B)
    GP<-rbinom(length(tgt),1,PA)
    G[tgt]<-GP
    
    tgt1a<-intersect(which(Y1==1),which(D==1))
    tgt0a<-intersect(which(Y1==0),which(D==1))
    tgt1b<-intersect(which(Y2==1),which(D==1))
    tgt0b<-intersect(which(Y2==0),which(D==1))
    list=sort(c(tgt1a,tgt0a))
    
    P111<-rep(NA,length(c(list)))  
    p1<-rep(NA,N)  
    p1[tgt1a]<- X[tgt1a,]%*%(alpha111a)+etan1[tgt1a]
    p1[tgt0a]<- -(X[tgt0a,]%*%(alpha111a)+etan1[tgt0a])
    p1<-p1[!is.na(p1)]
    p2<-rep(NA,N)  
    p2[tgt1b]<- (X[tgt1b,]%*%(alpha111b)+etan2[tgt1b])
    p2[tgt0b]<- -(X[tgt0b,]%*%(alpha111b)+etan2[tgt0b])
    p2<-p2[!is.na(p2)]
    P111<-pbvnorm(x=p1,y=p2,rho=rep(0,length=length(P111)))
    
    P101<-rep(NA,length(c(tgt1a,tgt0a)))
    p1<-rep(NA,N)  
    p1[tgt1a]<- (X[tgt1a,]%*%(alpha101a)+etan1[tgt1a])
    p1[tgt0a]<- -1*(X[tgt0a,]%*%(alpha101a)+etan1[tgt0a])
    p1<-p1[!is.na(p1)]
    p2<-rep(NA,N)  
    p2[tgt1b]<- (X[tgt1b,]%*%(alpha101b)+etan2[tgt1b])
    p2[tgt0b]<- -(X[tgt0b,]%*%(alpha101b)+etan2[tgt0b])
    p2<-p2[!is.na(p2)]
    P101<-pbvnorm(x=p1,y=p2,rho=rep(0,length=length(P101)))
    
    A=(1-pnorm(X[list,1:3]%*%t(beta)+chin1[tgt]))*pnorm(X[list,1:3]%*%t(gamma)+chin1[list])* P101 #likelihood of protected
    B=(1-pnorm(X[list,1:3]%*%t(beta)+chin1[tgt]))*(1-pnorm(X[list,1:3]%*%t(gamma)+chin1[list]))*P111#likelihood of always survivor
    
    #missing value imputation for individuals who survive with missing outcomes
    Y1[intersect(tgts1,which(G==2))]<-rbinom(length(intersect(tgts1,which(G==2))),1, pnorm(X[intersect(tgts1,which(G==2)),]%*%alpha111a+etan1[intersect(tgts1,which(G==2))]))
    Y1[intersect(tgts2,which(G==2))]<-rbinom(length(intersect(tgts2,which(G==2))),1, pnorm(X[intersect(tgts2,which(G==2)),]%*%alpha110a+etan1[intersect(tgts2,which(G==2))]))
    Y1[intersect(tgts1,which(G==1))]<-rbinom(length(intersect(tgts1,which(G==1))),1, pnorm(X[intersect(tgts1,which(G==1)),]%*%alpha101a+etan1[intersect(tgts1,which(G==1))]))
    
    Y2[intersect(tgts1,which(G==2))]<-rbinom(length(intersect(tgts1,which(G==2))),1, pnorm(X[intersect(tgts1,which(G==2)),]%*%alpha111b+etan2[intersect(tgts1,which(G==2))]))
    Y2[intersect(tgts2,which(G==2))]<-rbinom(length(intersect(tgts2,which(G==2))),1, pnorm(X[intersect(tgts2,which(G==2)),]%*%alpha110b+etan2[intersect(tgts2,which(G==2))]))
    Y2[intersect(tgts1,which(G==1))]<-rbinom(length(intersect(tgts1,which(G==1))),1, pnorm(X[intersect(tgts1,which(G==1)),]%*%alpha101b+etan2[intersect(tgts1,which(G==1))]))
    
    #missing value imputation for individuals with both survival status and outcomes missing
    surv[tgtss]<-ifelse(G[tgtss]==2,1,ifelse(G[tgtss]==1&D[tgtss]==1,1,0))
    tgtss1=which(surv[tgtss]==1&D[tgtss]==1)
    tgtss2=which(surv[tgtss]==1&D[tgtss]==0)
    
    Y1[intersect(tgtss[tgtss1],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111a+etan1[intersect(tgtss[tgtss1],which(G==2))]))
    Y1[intersect(tgtss[tgtss2],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss2],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110a+etan1[intersect(tgtss[tgtss2],which(G==2))]))
    Y1[intersect(tgtss[tgtss1],which(G==1))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==1))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101a+etan1[intersect(tgtss[tgtss1],which(G==1))]))
    
    Y2[intersect(tgtss[tgtss1],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111b+etan2[intersect(tgtss[tgtss1],which(G==2))]))
    Y2[intersect(tgtss[tgtss2],which(G==2))]<-rbinom(length(intersect(tgtss[tgtss2],which(G==2))),1, pnorm(X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110b+etan2[intersect(tgtss[tgtss2],which(G==2))]))
    Y2[intersect(tgtss[tgtss1],which(G==1))]<-rbinom(length(intersect(tgtss[tgtss1],which(G==1))),1, pnorm(X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101b+etan2[intersect(tgtss[tgtss1],which(G==1))]))
    
    Y1[intersect(miss,which(G==0))]<-NA
    Y1[intersect(miss,which(G==1&D==0))]<-NA
    Y2[intersect(miss,which(G==0))]<-NA
    Y2[intersect(miss,which(G==1&D==0))]<-NA
    
    #Latent variable updates
    #set Z/W using the known strata
    Z[G==0]<-rtruncnorm(length(which(G==0)),a=0,mean=(X[G==0,1:3]%*%c(beta)+chin1[G==0]),sd=1)
    Z[G>0]<-rtruncnorm(length(which(G>0)),b=0,mean=(X[G>0,1:3]%*%c(beta)+chin1[G>0]),sd=1)
    W[G==1]<-rtruncnorm(length(which(G==1)),a=0,mean=(X[G==1,1:3]%*%c(gamma)+chin1[G==1]),sd=1)
    W[G==2]<-rtruncnorm(length(which(G==2)),b=0,mean=(X[G==2,1:3]%*%c(gamma)+chin1[G==2]),sd=1)
    W[G==0]<-NA #set back to 0 is very crucial
    
    #set U111, U110 U101 using information
    U111a<-U110a<-U101a<-rep(NA,N)
    U111b<-U110b<-U101b<-rep(NA,N)
    
    tgt<-intersect(which(Y1==1), intersect(which(G==2),which(D==1)))
    U111a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha111a)+etan1[tgt],sd=1)
    tgt<-intersect(which(Y1==0), intersect(which(G==2),which(D==1)))
    U111a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha111a)+etan1[tgt],sd=1)
    tgt<-intersect(which(Y1==1), intersect(which(G==1),which(D==1)))
    U101a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha101a)+etan1[tgt],sd=1)
    tgt<-intersect(which(Y1==0), intersect(which(G==1),which(D==1)))
    U101a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha101a)+etan1[tgt],sd=1)
    tgt<-intersect(which(Y1==1), intersect(which(G==2),which(D==0)))
    U110a[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha110a)+etan1[tgt],sd=1)
    tgt<-intersect(which(Y1==0), intersect(which(G==2),which(D==0)))
    U110a[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha110a)+etan1[tgt],sd=1)
    
    tgt<-intersect(which(Y2==1), intersect(which(G==2),which(D==1)))
    U111b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha111b)+etan2[tgt],sd=1)
    tgt<-intersect(which(Y2==0), intersect(which(G==2),which(D==1)))
    U111b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha111b)+etan2[tgt],sd=1)
    tgt<-intersect(which(Y2==1), intersect(which(G==1),which(D==1)))
    U101b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha101b)+etan2[tgt],sd=1)
    tgt<-intersect(which(Y2==0), intersect(which(G==1),which(D==1)))
    U101b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha101b)+etan2[tgt],sd=1)
    tgt<-intersect(which(Y2==1), intersect(which(G==2),which(D==0)))
    U110b[tgt]<-rtruncnorm(length(tgt),a=0,mean=X[tgt,]%*%(alpha110b)+etan2[tgt],sd=1)
    tgt<-intersect(which(Y2==0), intersect(which(G==2),which(D==0)))
    U110b[tgt]<-rtruncnorm(length(tgt),b=0,mean=X[tgt,]%*%(alpha110b)+etan2[tgt],sd=1)
    
    #potential outcomes
    P1a<-pnorm(X[which(G==2),1:4]%*%alpha111a)
    P0a<-pnorm(X[which(G==2),1:4]%*%alpha110a)
    OR[1,i]<-(mean(P1a)/(1-mean(P1a)))/(mean(P0a)/(1-mean(P0a)))
    RR[1,i]<-mean(P1a)/mean(P0a)
    
    P1b<-pnorm(X[which(G==2),1:4]%*%alpha111b)
    P0b<-pnorm(X[which(G==2),1:4]%*%alpha110b)
    OR[2,i]<-(mean(P1b)/(1-mean(P1b)))/(mean(P0b)/(1-mean(P0b)))
    RR[2,i]<-mean(P1b)/mean(P0b)
    
    X1<-data.frame(cbind(P1a,P0a,P1b,P0b,cl[which(G==2)]))
    OR[3,i]<-mean(by((X1$X1/(1-X1$X1))/(X1$X2/(1-X1$X2)),X1$X5, mean))
    RR[3,i]<-mean(by(X1$X1/X1$X2,X1$X5, mean))
    OR[4,i]<-mean(by((X1$X3/(1-X1$X3))/(X1$X4/(1-X1$X4)),X1$X5, mean))
    RR[4,i]<-mean(by(X1$X3/X1$X4,X1$X5, mean))
    
    #store G percentage
    GV[,i]<-c(length(which(G==0)),length(which(G==1)),length(which(G==2)))/N
    GT[,i]<-G
    if (i%%500==0) print(i)
  }   
  return(list(ALPHA111a=ALPHA111a,ALPHA101a=ALPHA101a,ALPHA110a=ALPHA110a,ALPHA111b=ALPHA111b,ALPHA101b=ALPHA101b,ALPHA110b=ALPHA110b,
              BETA=BETA,GAMMA=GAMMA,PHI=PHI,TAU=TAU,RR=RR,OR=OR,GV=GV,GT=GT))
}


###################################
## 3. Simulation Studies ##########
###################################
sim<-function(S=10000,nsim=50,N=6000, s=60,seed=1234567){
  p=3
  ALPHA111a<-ALPHA101a<-ALPHA110a<- ALPHA111b<-ALPHA101b<-ALPHA110b<-array(NA,c(nsim,4,S))
  BETA<-GAMMA<-GV<-array(NA,c(nsim,3,S))
  TAU<-array(NA,c(nsim,4,S))
  OR<-RR<-array(NA,c(nsim,4,S))
  true_OR<-array(NA,c(nsim,4))
  true_RR<-array(NA,c(nsim,4))
  true_GV<-array(NA,c(nsim,3))
  PHI<-matrix(NA,nrow=nsim,ncol=S)
  
  for(i in 1:nsim){
    print(paste("nsim=",i))
    df<-gendata(N=N,s=s)
    N=nrow(df)
    true_OR[i,]<-c(df$ORa_i[1],df$ORb_i[1],df$ORa_c[1],df$ORb_c[1])
    true_RR[i,]<-c(df$RRa_i[1],df$RRb_i[1],df$RRa_c[1],df$RRb_c[1])
    true_GV[i,]<-c(sum(df$G==0)/N,
                   sum(df$G==1)/N,
                   sum(df$G==2)/N)
    
    otp<-try(mvsampler(df=df,S=S,dau1=0.005))  #avoid potential cases with extreme data/no intial value; not observed in simulations.
    if (class(otp) =="try-error"){
      next
    }
    ALPHA111a[i,,]<-otp$ALPHA111a
    ALPHA101a[i,,]<-otp$ALPHA101a
    ALPHA110a[i,,]<-otp$ALPHA110a
    ALPHA111b[i,,]<-otp$ALPHA111b
    ALPHA101b[i,,]<-otp$ALPHA101b
    ALPHA110b[i,,]<-otp$ALPHA110b
    BETA[i,,]<-otp$BETA
    GAMMA[i,,]<-otp$GAMMA
    
    PHI[i,]<-otp$PHI
    TAU[i,,]<-otp$TAU
    OR[i,,]<-otp$OR
    RR[i,,]<-otp$RR
    GV[i,,]<-otp$GV
    
  }
  sim<-list(ALPHA111a=ALPHA111a,ALPHA101a=ALPHA101a,ALPHA110a=ALPHA110a,
            ALPHA111b=ALPHA111b,ALPHA101b=ALPHA101b,ALPHA110b=ALPHA110b,
            BETA=BETA,GAMMA=GAMMA,PHI=PHI,
            GV=GV,TAU=TAU,OR=OR,RR=RR,true_OR=true_OR,true_RR=true_RR,true_GV=true_GV)
  
  return(sim)
}

######################################################
## 4. Output summary and ICC estimation ##############
######################################################
## simulation result summary
sumstat<-function(sim1,afterburn=2501:10000){
  pestimate<-c(
    mean(rowMeans(sim1$ALPHA111a[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111a[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111a[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111a[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA101a[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101a[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101a[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101a[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA110a[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110a[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110a[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110a[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA111b[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111b[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111b[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA111b[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA101b[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101b[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101b[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA101b[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$ALPHA110b[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110b[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110b[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$ALPHA110b[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$BETA[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$BETA[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$BETA[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$GAMMA[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GAMMA[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GAMMA[,3,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$TAU[,1, afterburn]),na.rm=T),
    mean(rowMeans(sim1$TAU[,2, afterburn]),na.rm=T),
    mean(rowMeans(sim1$TAU[,4, afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$PHI[,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$OR[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$OR[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$OR[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$OR[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$RR[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$RR[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$RR[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$RR[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$GV[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GV[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$GV[,3,afterburn]),na.rm=T))
  
  param<-c("alpha111a","alpha111a","alpha111a","alpha111a",
           "alpha101a","alpha101a","alpha101a","alpha101a",
           "alpha110a","alpha110a","alpha110a","alpha110a",
           "alpha111b","alpha111b","alpha111b","alpha111b",
           "alpha101b","alpha101b","alpha101b","alpha101b",
           "alpha110b","alpha110b","alpha110b","alpha110b",
           "beta","beta","beta",
           "gamma","gamma","gamma",
           "tau11","tau12","tau22","phi",
           "OR1","OR2","OR3","OR4",
           "RR1","RR2","RR3","RR4",
           "G0","G1","G2")
  
  true_val<-c(  c(alpha111a),  
                c(alpha101a),
                c(alpha110a),
                c(alpha111b) ,  
                c(alpha101b),
                c(alpha110b),
                c(beta),
                c(gamma),
                0.1,0.5*sqrt(0.01),0.1,
                phi,
                c(colMeans(sim1$true_OR)),
                c(colMeans(sim1$true_RR)),
                c(colMeans(sim1$true_GV)))
  
  rbias<-(pestimate-true_val)/true_val*100
  
  coverage<-c()
  a.true1=c(alpha111a)
  
  cla1=apply(sim1$ALPHA111a[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA111a[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA111a[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA111a[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  
  a1<-mean(a.true1[1]>=cla1[1,] & a.true1[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true1[2]>=cla2[1,] & a.true1[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true1[3]>=cla3[1,] & a.true1[3]<=cla3[2,],na.rm=T)
  a4<-mean(a.true1[4]>=cla4[1,] & a.true1[4]<=cla4[2,],na.rm=T)
  coverage<-c(coverage,a1,a2,a3,a4)
  
  a.true2=c(alpha101a)
  
  cla1=apply(sim1$ALPHA101a[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA101a[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA101a[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA101a[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true2[1]>=cla1[1,] & a.true2[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true2[2]>=cla2[1,] & a.true2[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true2[3]>=cla3[1,] & a.true2[3]<=cla3[2,],na.rm=T)
  a4<-mean(a.true2[4]>=cla4[1,] & a.true2[4]<=cla4[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3,a4)
  
  a.true3=c(alpha110a)
  
  cla1=apply(sim1$ALPHA110a[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA110a[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA110a[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA110a[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  
  a1<-mean(a.true3[1]>=cla1[1,] & a.true3[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true3[2]>=cla2[1,] & a.true3[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true3[3]>=cla3[1,] & a.true3[3]<=cla3[2,],na.rm=T)
  a4<-mean(a.true3[4]>=cla4[1,] & a.true3[4]<=cla4[2,],na.rm=T)
  coverage<-c(coverage,a1,a2,a3,a4)
  
  a.true1=c(alpha111b)
  
  cla1=apply(sim1$ALPHA111b[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA111b[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA111b[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA111b[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true1[1]>=cla1[1,] & a.true1[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true1[2]>=cla2[1,] & a.true1[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true1[3]>=cla3[1,] & a.true1[3]<=cla3[2,],na.rm=T)
  a4-mean(a.true1[4]>=cla4[1,] & a.true1[4]<=cla4[2,],na.rm=T)
  coverage<-c(coverage,a1,a2,a3,a4)
  
  a.true2=c(alpha101b)
  
  cla1=apply(sim1$ALPHA101b[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA101b[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA101b[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA101b[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true2[1]>=cla1[1,] & a.true2[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true2[2]>=cla2[1,] & a.true2[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true2[3]>=cla3[1,] & a.true2[3]<=cla3[2,],na.rm=T)
  a4<-mean(a.true2[4]>=cla4[1,] & a.true2[4]<=cla4[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3,a4)
  
  a.true3=c(alpha110b)
  
  cla1=apply(sim1$ALPHA110b[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$ALPHA110b[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$ALPHA110b[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$ALPHA110b[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  a1<-mean(a.true3[1]>=cla1[1,] & a.true3[1]<=cla1[2,],na.rm=T)
  a2<-mean(a.true3[2]>=cla2[1,] & a.true3[2]<=cla2[2,],na.rm=T)
  a3<-mean(a.true3[3]>=cla3[1,] & a.true3[3]<=cla3[2,],na.rm=T)
  a4<-mean(a.true3[4]>=cla4[1,] & a.true3[4]<=cla4[2,],na.rm=T)
  
  coverage<-c(coverage,a1,a2,a3,a4)
  
  #beta
  b.true=c(beta)
  cla1=apply(sim1$BETA[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$BETA[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$BETA[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  b1<-mean(b.true[1]>=cla1[1,] & b.true[1]<=cla1[2,],na.rm=T)
  b2<-mean(b.true[2]>=cla2[1,] & b.true[2]<=cla2[2,],na.rm=T)
  b3<-mean(b.true[3]>=cla3[1,] & b.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,b1,b2,b3)
  
  #gamma
  g.true=c(gamma)
  cla1=apply(sim1$GAMMA[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$GAMMA[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$GAMMA[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(g.true[1]>=cla1[1,] & g.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(g.true[2]>=cla2[1,] & g.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(g.true[3]>=cla3[1,] & g.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3)
  
  tau.true=c(0.1,0.5*sqrt(0.01),0.1)
  
  cla1=apply(sim1$TAU[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$TAU[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$TAU[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(tau.true[1]>=cla1[1,] & tau.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(tau.true[2]>=cla2[1,] & tau.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(tau.true[3]>=cla3[1,] & tau.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3)
  
  phi.true=phi
  cla1=apply(sim1$PHI[,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  g1<-mean(phi.true[1]>=cla1[1,] & phi.true[1]<=cla1[2,],na.rm=T)
  coverage<-c(coverage,g1)
  
  te.true=sim1$true_OR
  
  cla1=apply(sim1$OR[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$OR[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$OR[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$OR[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(te.true[,1]>=cla1[1,] & te.true[,1]<=cla1[2,],na.rm=T)
  g2<-mean(te.true[,2]>=cla2[1,] & te.true[,2]<=cla2[2,],na.rm=T)
  g3<-mean(te.true[,3]>=cla3[1,] & te.true[,3]<=cla1[2,],na.rm=T)
  g4<-mean(te.true[,4]>=cla4[1,] & te.true[,4]<=cla2[2,],na.rm=T)
  coverage<-c(coverage,g1,g2,g3,g4) 
  
  te.true=sim1$true_RR
  
  cla1=apply(sim1$RR[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$RR[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$RR[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$RR[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(te.true[,1]>=cla1[1,] & te.true[,1]<=cla1[2,],na.rm=T)
  g2<-mean(te.true[,2]>=cla2[1,] & te.true[,2]<=cla2[2,],na.rm=T)
  g3<-mean(te.true[,3]>=cla3[1,] & te.true[,3]<=cla1[2,],na.rm=T)
  g4<-mean(te.true[,4]>=cla4[1,] & te.true[,4]<=cla2[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3,g4) 
  
  g.true=sim1$true_GV
  
  cla1=apply(sim1$GV[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$GV[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$GV[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(g.true[,1]>=cla1[1,] & g.true[,1]<=cla1[2,],na.rm=T)
  g2<-mean(g.true[,2]>=cla2[1,] & g.true[,2]<=cla2[2,],na.rm=T)
  g3<-mean(g.true[,3]>=cla3[1,] & g.true[,3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3) 
  
  return(data.frame(param,true_val,pestimate,rbias,coverage))
}

icc_calculation<-function(sim1,afterburn=2501:1000){
  
  RHO1<-sim1$TAU[,1, afterburn]/(sim1$TAU[,1, afterburn]+1)
  RHO2<-sim1$TAU[,4, afterburn]/(sim1$TAU[,4, afterburn]+1)
  RHO112<-sim1$TAU[,2, afterburn]/(sqrt(sim1$TAU[,1, afterburn]+1)*
                                     sqrt(sim1$TAU[,4, afterburn]+1))
  
  RHO212<-(sim1$TAU[,2, afterburn]+0)/(sqrt(sim1$TAU[,1, afterburn]+1)*
                                         sqrt(sim1$TAU[,4, afterburn]+1))
  
  #mean ICC
  icc<-c(rho1<-mean(rowMeans(RHO1),na.rm=T),
         rho2<-mean(rowMeans(RHO2),na.rm=T),
         rho112<-mean(rowMeans(RHO112),na.rm=T),
         rho212<-mean(rowMeans(RHO212),na.rm=T))
  
  #coverage
  val.true=c(0.1/1.1,0.1/1.1,0.05/(sqrt(0.1+1)*sqrt(0.1+1)),0.05/(sqrt(0.1+1)*
                                                                    sqrt(0.1+1)))
  rbias<-(icc-val.true)/val.true*100
  
  cla1=apply(RHO1,1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(RHO2,1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(RHO112,1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(RHO212,1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(val.true[1]>=cla1[1,] & val.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(val.true[2]>=cla2[1,] & val.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(val.true[3]>=cla3[1,] & val.true[3]<=cla3[2,],na.rm=T)
  g4<-mean(val.true[4]>=cla4[1,] & val.true[4]<=cla4[2,],na.rm=T)
  
  coverage<-c(g1,g2,g3,g4)
  
  return(cbind( icc, val.true, rbias, coverage))
  
}

####################################################################
## 5. Example simulation execution and result summary ##############
####################################################################
output<-sim0(S=5000,nsim=50, N=1500, s=60)
sumstat(output,afterburn=2501:5000)
icc_calculation(output,afterburn=2501:5000)












