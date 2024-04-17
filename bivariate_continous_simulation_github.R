##################################################################################################################################################################################################################################################
##### Simulation code for the joint model with principal stratification model and outcome regression model with bivariate continuous outcomes and two study covariates in cluster randomized trials subject to missingness#####
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

gendata<-function(N=1500,s=60,cv=0.3){
  
  m=N/s  #mean cluster size
  tau2<-matrix(c(1,0.5*sqrt(2),0.5*sqrt(2),2),ncol=2,nrow=2) #cluster variance-covariance matrix
  sigma2<-matrix(c(5,0.5*sqrt(5*10),0.5*sqrt(5*10),10),ncol=2,nrow=2) #outcome variance-covariance matrix
  
  ##cluster ids and two covariates
  a0=1/cv^2
  b0=1/m/cv^2
  cl_size=round(rgamma(n=s,shape=a0,rate=b0))
  cl_size[which(cl_size<1)]<-1 #set the lower bound to 1
  
  cl<-rep(1:s,cl_size) #cluster index
  size=rep(cl_size,cl_size)

  N=length(cl)
  x1<-rnorm(N,0,10) #normal continuous 
  x2<-runif(N,-10,10)  #uniform continuous 
  X<-as.matrix(data.frame(rep(1,N),x1,x2,size))   #covariate matrix
  
  ##strata model parameters
  beta=c(-5.5,0.5,-0.7) #first Probit layer effect (strata proportion with 0.2 never survivors, 0.2 protected and 0.6 always survivors)
  gamma=c(-5.8,-0.6,0.4) #second Probit layer effect(strata proportion with 0.2 never survivors, 0.2 protected and 0.6 always survivors)

  phi2<-1 #cluster variance

  #Cluster effect 
  chi<-rnorm(s,0,sqrt(phi2)) 
  chin<-rep(chi,cl_size) 
  
  #principal strata membership (0: always non-survivor,1: protected, 2: always survivor)
  pp<-pnorm(X[,1:3]%*%beta+chin)
  G00 <- rbinom(N, 1, pp)
  G01<-rep(0,N)
  tgt<-which(G00==0)
  p1<-pnorm(X[tgt,1:3]%*%gamma+chin[tgt])
  G01[tgt]<-rbinom(length(tgt), 1, p1)
  G11<-rep(1,N)-G01-G00
  G=1*G01+2*G11
  
  #treatment assignment
  d<-rep(sample(rep(c(0,1),each=s/2),replace=F),cl_size)
  
  #outcome model parameters
  alpha111a<-c(-13,-0.5,-0.2,0.3)
  alpha101a<-c(2,2,2,1)
  alpha110a<-c(14,0.5,-0.4,-0.4)
  
  alpha111b<-c(-11,-0.4,0.3,0.3)
  alpha101b<-c(-9,-1,0.8,-1)
  alpha110b<-c(12,0.4,0.4,-0.3)
  resid<-rmvnorm(N,c(0,0),sigma2) #outcome variance-covariance matrix
  eta<-rmvnorm(s,c(0,0),tau2)  #cluster effect
  etan<-cbind(rep(eta[,1],cl_size),rep(eta[,2],cl_size)) #cluster effect at individual level
  
  #continuous outcomes
  # truncation by death
  Y1<-rep(NA,N) #outcome (NA are death)
  Y1[intersect(which(G11==1),which(d==1))]<-X[intersect(which(G11==1),which(d==1)),]%*%alpha111a
  Y1[intersect(which(G01==1),which(d==1))]<-X[intersect(which(G01==1),which(d==1)),]%*%alpha101a
  Y1[intersect(which(G11==1),which(d==0))]<-X[intersect(which(G11==1),which(d==0)),]%*%alpha110a
  Y1<-Y1+etan[,1]+resid[,1]
  
  Y2<-rep(NA,N) #outcome (NA are death)
  Y2[intersect(which(G11==1),which(d==1))]<-X[intersect(which(G11==1),which(d==1)),]%*%alpha111b
  Y2[intersect(which(G01==1),which(d==1))]<-X[intersect(which(G01==1),which(d==1)),]%*%alpha101b
  Y2[intersect(which(G11==1),which(d==0))]<-X[intersect(which(G11==1),which(d==0)),]%*%alpha110b
  Y2<-Y2+etan[,2]+resid[,2]
  

  s<-ifelse(is.na(Y1),0,1) # survival status

  ##survivor average causal effect
  a<-X[which(G==2),]%*%alpha111a-X[which(G==2),]%*%alpha110a
  b<-X[which(G==2),]%*%alpha111b-X[which(G==2),]%*%alpha110b
  TE1_i<-colMeans(a) #outcome 1 SIACE
  TE2_i<-colMeans(b) #outcome 2 SIACE
  
  X1<-data.frame(cbind(a,b,cl[which(G==2)]))
  TE1_c<-mean(by(X1$X1,X1$X3, mean)) #outcome 1 SCACE
  TE2_c<-mean(by(X1$X2,X1$X3, mean)) #outcome 2 SCACE
  
  #missingness model
  z1=c(10.5,2,-0.1,0.3)
  z3=c(-0.25,0.5,-0.5,0.2)
  ##missing both survival status and outcomes
  pr1=pnorm(X[,1:4]%*%z1)
  k<-rbinom(length((s)),1,pr1)
  s<-ifelse(k==1,s,NA)
  Y1[which(is.na(s))]<-NA
  Y2[which(is.na(s))]<-NA
  
  # only missing outcome
  mb<-which(is.na(Y1)&is.na(s)) 
  pr3=pnorm(X[-mb,1:4]%*%z3)
  k2<-rbinom((length(s)-length(mb)),1,pr3)
  Y1[-mb]<-ifelse(k2==1,Y1[-mb],NA)
  Y2[-mb]<-ifelse(k2==1,Y2[-mb],NA)
  
  #For G, the true status of 0, 1 or 2 will not be known, and death in treatment group and survival in control group have known G status   
  return(data.frame(cbind(Y1,Y2,x1,x2,d,cl,s,size,G,TE1_i,TE2_i,TE1_c,TE2_c)))
}

###################################
## 2. MCMC model estimation  ######
###################################
#df: the simulated data 
#S: length of chain

mvsampler<-function(df,S=5000){
  # extract data information
  surv=df$s # survival status
  N=nrow(df) # number of individuals
  Y1=df$Y1 #outcome 1
  Y2=df$Y2 #outcome 2
  cl=df$cl # cluster ID
  
  s<-length(unique(cl)) # number of clusters
  cs<-as.numeric(table(cl)) # size of cluster
  
  D<-df$d #treatment assignment
  
  ##starting values
  #outcome model
  
  alpha111a<-c(-13,-0.5,-0.2,0.3)
  alpha101a<-c(2,2,2,1)
  alpha110a<-c(14,0.5,-0.4,-0.4)
  
  alpha111b<-c(-11,-0.4,0.3,0.3)
  alpha101b<-c(-9,-1,0.8,-1)
  alpha110b<-c(12,0.4,0.4,-0.3)

  beta=c(-5.5,0.5,-0.7) #first Probit layer effect (strata proportion with 0.2 never survivors, 0.2 protected and 0.6 always survivors)
  gamma=c(-5.8,-0.6,0.4) #second Probit layer effect(strata proportion with 0.2 never survivors, 0.2 protected and 0.6 always survivors)
  #beta=c(-8.5,0.5,-0.7)  #first Probit layer effect(strata proportion with 0.1 never survivors, 0.1 protected and 0.8 always survivors)
  #gamma=c(-8.8,-0.6,0.4) #second Probit layer effect(strata proportion with 0.1 never survivors, 0.1 protected and 0.8 always survivors)
  
  Tau<-matrix(c(1,0.5*sqrt(2),0.5*sqrt(2),2),ncol=2,nrow=2)
  Sigma<-matrix(c(5,0.5*sqrt(5*10),0.5*sqrt(5*10),10),ncol=2,nrow=2)
  # random effect at individual level
  eta<-rmvnorm(s,c(0,0),Tau)
  etan1<-rep(eta[,1],cs)
  etan2<-rep(eta[,2],cs)
  
  ##weakly informative priors
  #outcome model regression coefficient
  Sigma111<-Sigma101<-Sigma110<-diag(2*(p+1))*1000
  a=0.0001  #prior variance terms
  b=0.0001
  cc=0.0001
  dd=0.0001
  
  #initial value for G  
  G=df$G
  
  p=3 #number of covariates(stra model)
  X<-as.matrix(cbind(rep(1,N),df$x1,df$x2,df$size)) #covariate matrix
  
  #strata model
  #regression coefficients
  beta0=gamma0=rep(0,p) #normal prior with mean zero
  Lambda=Gamma=diag(p)*10000
  # random effects
  phi=1
  chi1<-rnorm(s,0,sqrt(phi))
  chin1<-rep(chi1,cs)
  
  g=0.0001  #prior variance terms
  h=0.0001
  
  
  ##store chains for parameters/estimates with length S
  ALPHA111a<-ALPHA101a<-ALPHA110a<-matrix(NA,p+1,S)
  ALPHA111b<-ALPHA101b<-ALPHA110b<-matrix(NA,p+1,S)
  SIGMA<-TAU<-matrix(NA,4,S)
  
  TE<-matrix(NA,4,S) #sample TE by estimating the potential outcomes of G=2.
  
  BETA<-GAMMA<-matrix(NA,p,S)
  PHI<-rep(0,S)
  GV<-matrix(NA,3,S)
  
  tgts1=which(is.na(Y1)&surv==1&D==1) 
  tgts2=which(is.na(Y1)&surv==1&D==0)
  tgtss=which(is.na(Y1)&is.na(surv)) 
  miss=c(which(is.na(Y1)&is.na(surv)),which(is.na(Y1)&surv==1))
  nonmiss=c(1:N)[-miss] 
  Nst<-length(intersect(which(!is.na(Y1)),nonmiss))
  
  #missing value imputation for individuals who survive with missing outcomes
  residm1<-rmvnorm(N,c(0,0),Sigma)
  Y1[intersect(tgts1,which(G==2))]<-X[intersect(tgts1,which(G==2)),]%*%alpha111a+etan1[intersect(tgts1,which(G==2))]+residm1[intersect(tgts1,which(G==2)),1]
  Y1[intersect(tgts2,which(G==2))]<-X[intersect(tgts2,which(G==2)),]%*%alpha110a+etan1[intersect(tgts2,which(G==2))]+residm1[intersect(tgts2,which(G==2)),1]
  Y1[intersect(tgts1,which(G==1))]<-X[intersect(tgts1,which(G==1)),]%*%alpha101a+etan1[intersect(tgts1,which(G==1))]+residm1[intersect(tgts1,which(G==1)),1]
  
  Y2[intersect(tgts1,which(G==2))]<-X[intersect(tgts1,which(G==2)),]%*%alpha111b+etan2[intersect(tgts1,which(G==2))]+residm1[intersect(tgts1,which(G==2)),2]
  Y2[intersect(tgts2,which(G==2))]<-X[intersect(tgts2,which(G==2)),]%*%alpha110b+etan2[intersect(tgts2,which(G==2))]+residm1[intersect(tgts2,which(G==2)),2]
  Y2[intersect(tgts1,which(G==1))]<-X[intersect(tgts1,which(G==1)),]%*%alpha101b+etan2[intersect(tgts1,which(G==1))]+residm1[intersect(tgts1,which(G==1)),2]
  
  ## impute S and Y for both missing
  surv[tgtss]<-ifelse(G[tgtss]==2,1,ifelse(G[tgtss]==1&D[tgtss]==1,1,0))
  tgtss1=which(surv[tgtss]==1&D[tgtss]==1)
  tgtss2=which(surv[tgtss]==1&D[tgtss]==0)
  
  Y1[intersect(tgtss[tgtss1],which(G==2))]<-X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111a+etan1[intersect(tgtss[tgtss1],which(G==2))]+residm1[intersect(tgtss[tgtss1],which(G==2)),1]
  Y1[intersect(tgtss[tgtss2],which(G==2))]<-X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110a+etan1[intersect(tgtss[tgtss2],which(G==2))]+residm1[intersect(tgtss[tgtss2],which(G==2)),1]
  Y1[intersect(tgtss[tgtss1],which(G==1))]<-X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101a+etan1[intersect(tgtss[tgtss1],which(G==1))]+residm1[intersect(tgtss[tgtss1],which(G==1)),1]
  
  Y2[intersect(tgtss[tgtss1],which(G==2))]<-X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111b+etan2[intersect(tgtss[tgtss1],which(G==2))]+residm1[intersect(tgtss[tgtss1],which(G==2)),2]
  Y2[intersect(tgtss[tgtss2],which(G==2))]<-X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110b+etan2[intersect(tgtss[tgtss2],which(G==2))]+residm1[intersect(tgtss[tgtss2],which(G==2)),2]
  Y2[intersect(tgtss[tgtss1],which(G==1))]<-X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101b+etan2[intersect(tgtss[tgtss1],which(G==1))]+residm1[intersect(tgtss[tgtss1],which(G==1)),2]
  
  # latent parameter from truncated normal distribution 
  Z<-rep(NA,N)
  W<-rep(NA,N)
  Z[G==0]<-rtruncnorm(length(which(G==0)),a=0,mean=(X[G==0,1:3]%*%c(beta)+chin1[G==0]),sd=1)
  Z[G>0]<-rtruncnorm(length(which(G>0)),b=0,mean=(X[G>0,1:3]%*%c(beta)+chin1[G>0]),sd=1)
  W[G==1]<-rtruncnorm(length(which(G==1)),a=0,mean=(X[G==1,1:3]%*%c(gamma)+chin1[G==1]),sd=1)
  W[G==2]<-rtruncnorm(length(which(G==2)),b=0,mean=(X[G==2,1:3]%*%c(gamma)+chin1[G==2]),sd=1)
  W[G==0]<-NA #set back to 0 is very crucial
  
  
  for(i in 1:S){
    #update alpha111
    alpha111<-c(alpha111a,alpha111b)
    tgt<-intersect(which(D==1),which(G==2)) 
    
    InvSigma<-solve(Sigma)
    a=InvSigma[1,1]
    b=c=InvSigma[1,2]
    d=InvSigma[2,2]
    
    res111<-cbind(Y1[tgt],Y2[tgt])-cbind(etan1[tgt],etan2[tgt])
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
    
    res101<-cbind(Y1[tgt],Y2[tgt])-cbind(etan1[tgt],etan2[tgt])
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
    
    res110<-cbind(Y1[tgt],Y2[tgt])-cbind(etan1[tgt],etan2[tgt])
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
    
    #update eta1
    for(j in 1:s){
      if (D[which(cl==j)[1]]==0){
        tgt3<-intersect(which(G==2),intersect(which(D==0),which(cl==j)))
        
        V_eta=solve(length(tgt3)*InvSigma + solve(Tau))
        res_eta=matrix(cbind(Y1[tgt3],Y2[tgt3])-cbind(X[tgt3,]%*%alpha110a, X[tgt3,]%*%alpha110b),length(tgt3),2)
        M_eta=V_eta%*%(InvSigma %*% colSums(res_eta))
        
        eta[j,]<-as.numeric(rmvnorm(1,M_eta,V_eta))
        
      }else{
        tgt1<-intersect(which(G==2),intersect(which(D==1),which(cl==j)))
        tgt2<-intersect(which(G==1),intersect(which(D==1),which(cl==j)))
        
        V_eta=solve((length(tgt1)+length(tgt2))*InvSigma + solve(Tau))
        
        res_eta1=matrix(cbind(Y1[tgt1],Y2[tgt1])-cbind(X[tgt1,]%*%alpha111a, X[tgt1,]%*%alpha111b),length(tgt1),2)
        res_eta2=matrix(cbind(Y1[tgt2],Y2[tgt2])-cbind(X[tgt2,]%*%alpha101a, X[tgt2,]%*%alpha101b),length(tgt2),2)
        M_eta=V_eta%*%(InvSigma %*% (colSums(res_eta1) + colSums(res_eta2)))
        
        eta[j,]<-as.numeric(rmvnorm(1,M_eta,V_eta))
      }
    }
    etan1<-rep(eta[,1],cs)
    etan2<-rep(eta[,2],cs)
    
    #update Tau
    Tau<-solve(rwish(1,2+s,solve(matrix(c(1,0,0,1),2,2)+crossprod(eta))))
    TAU[,i]<-as.numeric(Tau)
    
    #update Sigma
    tgt1<-intersect(which(G==2),which(D==1))
    tgt2<-intersect(which(G==1),which(D==1))
    tgt3<-intersect(which(G==2),which(D==0))
    
    resid1<-c((Y1[tgt1]-X[tgt1,]%*%alpha111a-etan1[tgt1]),
              (Y1[tgt2]-X[tgt2,]%*%alpha101a-etan1[tgt2]),
              (Y1[tgt3]-X[tgt3,]%*%alpha110a-etan1[tgt3]))
    resid2<-c((Y2[tgt1]-X[tgt1,]%*%alpha111b-etan2[tgt1]),
              (Y2[tgt2]-X[tgt2,]%*%alpha101b-etan2[tgt2]),
              (Y2[tgt3]-X[tgt3,]%*%alpha110b-etan2[tgt3]))
    
    resid<-cbind(resid1,resid2)
    Nst<-length(which(!is.na(Y1)))

    Sigma<-solve(rwish(1,2+Nst,solve(matrix(c(1,0,0,1),2,2)+crossprod(resid))))
    SIGMA[,i]<-as.numeric(Sigma)
    
    #potential outcomes
    
    Y11<-X[G==2,]%*%alpha111a
    Y10<-X[G==2,]%*%alpha110a
    te1<-mean(Y11-Y10) 
    TE[1,i]<-te1 
    
    Y21<-X[G==2,]%*%alpha111b
    Y20<-X[G==2,]%*%alpha110b
    te2<-mean(Y21-Y20) 
    TE[2,i]<-te2 
    
    X1<-data.frame(cbind(Y11-Y10,Y21-Y20,cl[which(G==2)]))
    te3<-mean(by(X1$X1,X1$X3, mean))
    te4<-mean(by(X1$X2,X1$X3, mean))
    TE[3,i]<-te3
    TE[4,i]<-te4
    
    #update beta/gamma
    v=solve(crossprod(X[,1:3])+solve(Lambda))
    m=v%*%(crossprod(X[,1:3],(Z-chin1))+solve(Lambda)%*%beta0)
    beta=rmvnorm(1,m,v)
    BETA[,i]<-beta
    
    tgt1<-which(G>0)
    v=solve(crossprod(X[tgt1,1:3])+solve(Gamma))
    m=v%*%(crossprod(X[tgt1,1:3],(W[tgt1]-chin1[tgt1]))+solve(Gamma)%*%gamma0)
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
    tgt2<-intersect(which(is.na(Y1)),which(D==0))
    
    A=pnorm(X[tgt2,1:3]%*%t(beta)+chin1[tgt2]) #probability always death/(always death+protected+always survival)
    B=pnorm(X[tgt2,1:3]%*%t(gamma)+chin1[tgt2]) #probability protected/(protected+always survival)
    PA<-(1-A)*B/(A+(1-A)*B)
    GP<-rbinom(length(tgt2),1,PA)
    G[tgt2]<-GP
    
    tgt1<-intersect(which(!is.na(Y1)),which(D==1))
    P111<-rep(NA,length(c(tgt1)))
    P111<-pbvnorm(x=(Y1[tgt1]-X[tgt1,]%*%alpha111a-etan1[tgt1])/sqrt(Sigma[1,1]),y=(Y2[tgt1]-X[tgt1,]%*%alpha111b-etan2[tgt1])/sqrt(Sigma[2,2]),rho=rep(Sigma[1,2]/(sqrt(Sigma[1,1]*Sigma[2,2])),length(c(tgt1))))
    P101<-rep(NA,length(c(tgt1)))
    P101<-pbvnorm(x=(Y1[tgt1]-X[tgt1,]%*%alpha101a-etan1[tgt1])/sqrt(Sigma[1,1]),y=(Y2[tgt1]-X[tgt1,]%*%alpha101b-etan2[tgt1])/sqrt(Sigma[2,2]),rho=rep(Sigma[1,2]/(sqrt(Sigma[1,1]*Sigma[2,2])),length(c(tgt1))))
    A=(1-pnorm(X[c(tgt1),1:3]%*%t(beta)+chin1[tgt1]))*pnorm(X[c(tgt1),1:3]%*%t(gamma)+chin1[tgt1])* P101 #likelihood of protected
    B=(1-pnorm(X[c(tgt1),1:3]%*%t(beta)+chin1[tgt1]))*(1-pnorm(X[c(tgt1),1:3]%*%t(gamma)+chin1[tgt1]))*P111#likelihood of always survivor
    
    r=rbinom(length(c(tgt1)),1,(A/(A+B)))
    GP<-rep(2,length(c(tgt1)))-r
    G[c(tgt1)]<-GP

    #impute for G for missing survival status and missing outcome
    pp<-pnorm(X[tgtss,1:3]%*%t(beta)+chin1[tgtss])
    G00 <- rbinom(length(tgtss), 1, pp)
    G01<-rep(0,length(tgtss))
    tgt<-which(G00==0)
    p1<-pnorm(X[tgtss,][tgt,1:3]%*%t(gamma)+chin1[tgtss][tgt])
    G01[tgt]<-rbinom(length(tgt), 1, p1)
    G11<-rep(1,length(tgtss))-G01-G00
    G[tgtss]=1*G01+2*G11
    
    #store G percentage
    GV[,i]<-c(length(which(G==0)),length(which(G==1)),length(which(G==2)))/N
    if (i%%5000==0) print(i)
    
    #update missingness in outcomes
    residm1<-rmvnorm(N,c(0,0),Sigma)
    
    #missing value imputation for individuals who survive with missing outcomes
    Y1[intersect(tgts1,which(G==2))]<-X[intersect(tgts1,which(G==2)),]%*%alpha111a+etan1[intersect(tgts1,which(G==2))]+residm1[intersect(tgts1,which(G==2)),1]
    Y1[intersect(tgts2,which(G==2))]<-X[intersect(tgts2,which(G==2)),]%*%alpha110a+etan1[intersect(tgts2,which(G==2))]+residm1[intersect(tgts2,which(G==2)),1]
    Y1[intersect(tgts1,which(G==1))]<-X[intersect(tgts1,which(G==1)),]%*%alpha101a+etan1[intersect(tgts1,which(G==1))]+residm1[intersect(tgts1,which(G==1)),1]
    
    Y2[intersect(tgts1,which(G==2))]<-X[intersect(tgts1,which(G==2)),]%*%alpha111b+etan2[intersect(tgts1,which(G==2))]+residm1[intersect(tgts1,which(G==2)),2]
    Y2[intersect(tgts2,which(G==2))]<-X[intersect(tgts2,which(G==2)),]%*%alpha110b+etan2[intersect(tgts2,which(G==2))]+residm1[intersect(tgts2,which(G==2)),2]
    Y2[intersect(tgts1,which(G==1))]<-X[intersect(tgts1,which(G==1)),]%*%alpha101b+etan2[intersect(tgts1,which(G==1))]+residm1[intersect(tgts1,which(G==1)),2]
    
    
    #missing value imputation for individuals with both survival status and outcomes missing
    surv[tgtss]<-ifelse(G[tgtss]==2,1,ifelse(G[tgtss]==1&D[tgtss]==1,1,0))
    tgtss1=which(surv[tgtss]==1&D[tgtss]==1)
    tgtss2=which(surv[tgtss]==1&D[tgtss]==0)
    
    Y1[intersect(tgtss[tgtss1],which(G==2))]<-X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111a+etan1[intersect(tgtss[tgtss1],which(G==2))]+residm1[intersect(tgtss[tgtss1],which(G==2)),1]
    Y1[intersect(tgtss[tgtss2],which(G==2))]<-X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110a+etan1[intersect(tgtss[tgtss2],which(G==2))]+residm1[intersect(tgtss[tgtss2],which(G==2)),1]
    Y1[intersect(tgtss[tgtss1],which(G==1))]<-X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101a+etan1[intersect(tgtss[tgtss1],which(G==1))]+residm1[intersect(tgtss[tgtss1],which(G==1)),1]
    
    Y2[intersect(tgtss[tgtss1],which(G==2))]<-X[intersect(tgtss[tgtss1],which(G==2)),]%*%alpha111b+etan2[intersect(tgtss[tgtss1],which(G==2))]+residm1[intersect(tgtss[tgtss1],which(G==2)),2]
    Y2[intersect(tgtss[tgtss2],which(G==2))]<-X[intersect(tgtss[tgtss2],which(G==2)),]%*%alpha110b+etan2[intersect(tgtss[tgtss2],which(G==2))]+residm1[intersect(tgtss[tgtss2],which(G==2)),2]
    Y2[intersect(tgtss[tgtss1],which(G==1))]<-X[intersect(tgtss[tgtss1],which(G==1)),]%*%alpha101b+etan2[intersect(tgtss[tgtss1],which(G==1))]+residm1[intersect(tgtss[tgtss1],which(G==1)),2]
    
    Y1[intersect(miss,which(G==0))]<-NA
    Y1[intersect(miss,which(G==1&D==0))]<-NA
    Y2[intersect(miss,which(G==0))]<-NA
    Y2[intersect(miss,which(G==1&D==0))]<-NA
    
    #latent variable updates
    #set Z/W using the known strata
    Z[G==0]<-rtruncnorm(length(which(G==0)),a=0,mean=(X[G==0,1:3]%*%c(beta)+chin1[G==0]),sd=1)
    Z[G>0]<-rtruncnorm(length(which(G>0)),b=0,mean=(X[G>0,1:3]%*%c(beta)+chin1[G>0]),sd=1)
    W[G==1]<-rtruncnorm(length(which(G==1)),a=0,mean=(X[G==1,1:3]%*%c(gamma)+chin1[G==1]),sd=1)
    W[G==2]<-rtruncnorm(length(which(G==2)),b=0,mean=(X[G==2,1:3]%*%c(gamma)+chin1[G==2]),sd=1)
    W[G==0]<-NA #set back to 0 is very crucial
    
  }
  
  sim_result<-list(ALPHA111a=ALPHA111a,ALPHA101a=ALPHA101a,ALPHA110a=ALPHA110a,ALPHA111b=ALPHA111b,ALPHA101b=ALPHA101b,ALPHA110b=ALPHA110b,
                   BETA=BETA,GAMMA=GAMMA,TAU=TAU,SIGMA=SIGMA,PHI=PHI,GV=GV,TE=TE)
  return(sim_result)
}

###################################
## 3. Simulation Studies ##########
###################################
sim<-function(S=10000,nsim=50,N=6000, s=60,seed=1234567){
  p=3
  GV<-matrix(NA,3,S)

  ALPHA111a<-ALPHA101a<-ALPHA110a<-
    ALPHA111b<-ALPHA101b<-ALPHA110b<-array(NA,c(nsim,4,S))
  
  BETA<-GAMMA<-GV<-array(NA,c(nsim,3,S))
  TAU<-SIGMA<-array(NA,c(nsim,4,S))
  TE<-array(NA,c(nsim,4,S))
  true_TE<-array(NA,c(nsim,4))
  true_GV<-array(NA,c(nsim,3))
  PHI<-matrix(NA,nrow=nsim,ncol=S)
  
  for(i in 1:nsim){
    print(paste("replication",i))
    df<-gendata(N=N,s=s)
    N=nrow(df)
    true_TE[i,]<-c(df$TE1_i[1],df$TE2_i[1],df$TE1_c[1],df$TE2_c[1])
    true_GV[i,]<-c(sum(df$G==0)/N,
                   sum(df$G==1)/N,
                   sum(df$G==2)/N)
    otp<-try(mvsampler(df=df,S=S,dau1=0.005)) #avoid potential cases with extreme data/no intial value; not observed in simulations.
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
    SIGMA[i,,]<-otp$SIGMA
    TE[i,,]<-otp$TE
    GV[i,,]<-otp$GV
  }
  
  sim<-list(ALPHA111a=ALPHA111a,ALPHA101a=ALPHA101a,ALPHA110a=ALPHA110a,
            ALPHA111b=ALPHA111b,ALPHA101b=ALPHA101b,ALPHA110b=ALPHA110b,
            BETA=BETA,GAMMA=GAMMA,PHI=PHI,
            GV=GV,TAU=TAU,SIGMA=SIGMA,TE=TE,true_TE=true_TE,true_GV=true_GV)
  
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
    
    mean(rowMeans(sim1$SIGMA[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$SIGMA[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$SIGMA[,4,afterburn]),na.rm=T),
    
    mean(rowMeans(sim1$TE[,1,afterburn]),na.rm=T),
    mean(rowMeans(sim1$TE[,2,afterburn]),na.rm=T),
    mean(rowMeans(sim1$TE[,3,afterburn]),na.rm=T),
    mean(rowMeans(sim1$TE[,4,afterburn]),na.rm=T),
    
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
           "sigma11","sigma12","sigma22",
           "te1","te2","te3","te4",
           "G0","G1","G2")
  
  true_val<-c(  c(alpha111a),  
                c(alpha101a),
                c(alpha110a),
                c(alpha111b) ,  
                c(alpha101b),
                c(alpha110b),
                c(beta),
                c(gamma),
                1,0.5*sqrt(2),2,
                phi,
                5,0.5*sqrt(50),10,
                c(colMeans(sim1$true_TE)),
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
  
  tau.true=c(1,0.5*sqrt(2),2)
  
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
  
  sigma.true=c(5,0.5*sqrt(50),10)
  
  cla1=apply(sim1$SIGMA[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$SIGMA[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$SIGMA[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
  g1<-mean(sigma.true[1]>=cla1[1,] & sigma.true[1]<=cla1[2,],na.rm=T)
  g2<-mean(sigma.true[2]>=cla2[1,] & sigma.true[2]<=cla2[2,],na.rm=T)
  g3<-mean(sigma.true[3]>=cla3[1,] & sigma.true[3]<=cla3[2,],na.rm=T)
  
  coverage<-c(coverage,g1,g2,g3)   
  
  te.true=sim1$true_TE
  
  cla1=apply(sim1$TE[,1,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla2=apply(sim1$TE[,2,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla3=apply(sim1$TE[,3,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  cla4=apply(sim1$TE[,4,afterburn],1,function(x)quantile(x,c(0.025,0.975),na.rm=T))
  
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
  
  
  
  pestimate_var<-c(
    mean(rowVars(sim1$ALPHA111a[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111a[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111a[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111a[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$ALPHA101a[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101a[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101a[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101a[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$ALPHA110a[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110a[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110a[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110a[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$ALPHA111b[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111b[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111b[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA111b[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$ALPHA101b[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101b[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101b[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA101b[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$ALPHA110b[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110b[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110b[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$ALPHA110b[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$BETA[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$BETA[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$BETA[,3,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$GAMMA[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$GAMMA[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$GAMMA[,3,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$TAU[,1, afterburn]),na.rm=T),
    mean(rowVars(sim1$TAU[,2, afterburn]),na.rm=T),
    mean(rowVars(sim1$TAU[,4, afterburn]),na.rm=T),
    
    mean(rowVars(sim1$PHI[,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$SIGMA[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$SIGMA[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$SIGMA[,4,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$GV[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$GV[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$GV[,3,afterburn]),na.rm=T),
    
    mean(rowVars(sim1$TE[,1,afterburn]),na.rm=T),
    mean(rowVars(sim1$TE[,2,afterburn]),na.rm=T),
    mean(rowVars(sim1$TE[,3,afterburn]),na.rm=T),
    mean(rowVars(sim1$TE[,4,afterburn]),na.rm=T)
  )
  
  return(data.frame(param,true_val,pestimate,rbias,coverage,pestimate_var))
}

icc_calculation<-function(sim1,afterburn=2501:1000){

  RHO1<-sim1$TAU[,1, afterburn]/(sim1$TAU[,1, afterburn]+sim1$SIGMA[,1, afterburn])
  RHO2<-sim1$TAU[,4, afterburn]/(sim1$TAU[,4, afterburn]+sim1$SIGMA[,4, afterburn])
  RHO112<-sim1$TAU[,2, afterburn]/(sqrt(sim1$TAU[,1, afterburn]+sim1$SIGMA[,1, afterburn])*
                                     sqrt(sim1$TAU[,4, afterburn]+sim1$SIGMA[,4, afterburn]))
  
  RHO212<-(sim1$TAU[,2, afterburn]+sim1$SIGMA[,2, afterburn])/(sqrt(sim1$TAU[,1, afterburn]+sim1$SIGMA[,1, afterburn])*
                                                                 sqrt(sim1$TAU[,4, afterburn]+sim1$SIGMA[,4, afterburn]))
  
  #mean ICC
  icc<-c(rho1<-mean(rowMeans(RHO1),na.rm=T),
         rho2<-mean(rowMeans(RHO2),na.rm=T),
         rho112<-mean(rowMeans(RHO112),na.rm=T),
         rho212<-mean(rowMeans(RHO212),na.rm=T))
  
  #coverage
  val.true=c(1/6,2/12,0.08333333,0.5)
  
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
  
  pestimate_var<-c(mean(rowVars(RHO1),na.rm=T),
                   mean(rowVars(RHO2),na.rm=T),
                   mean(rowVars(RHO112),na.rm=T),
                   mean(rowVars(RHO212),na.rm=T))
  
  return(cbind( icc, val.true, rbias, coverage,pestimate_var))
  
}

####################################################################
## 5. Example simulation execution and result summary ##############
####################################################################
output<-sim0(S=5000,nsim=50, N=1500, s=60)
sumstat(output,afterburn=2501:5000)
icc_calculation(output,afterburn=2501:5000)








