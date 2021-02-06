
# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
setwd("/home/antonio/MEGA/Lavoro_sync/Working projects/AC_fuzzy-rating/simulation_study/")
set.seed(200105)


# Set design --------------------------------------------------------------
Ix = c(50,150,500) #low, medium, large sample sizes
Jx = c(5,15) 
Mx = c(3,5)
betax = c(-1.5,-4.5,-8.5) #low values indicate that the RT model is too noisy
design = expand.grid(Ix,Jx,betax,Mx); colnames(design)=c("I","J","beta","M")


# Set sim parameters ------------------------------------------------------
Eta=Omega=list()
for(k in 1:length(Ix)){
  Eta[[k]]=rnorm(n=Ix[k])
  Omega[[k]]=rnorm(n=Ix[k])
}
names(Eta)=names(Omega)=Ix

Alpha=Gamma=list("3"=list(),"5"=list())
for(k in 1:length(Jx)){
  Alpha[[1]][[k]]=pracma::repmat(rnorm(Jx[k],0,1),Mx[1],1)
  Alpha[[2]][[k]]=pracma::repmat(rnorm(Jx[k],0,1),Mx[2],1)
  Gamma[[1]][[k]]=runif(Jx[k],1,9)
  Gamma[[2]][[k]]=runif(Jx[k],1,9)
}
names(Alpha[[1]])=names(Alpha[[2]])=Jx
names(Gamma[[1]])=names(Gamma[[2]])=Jx


# Simulate data -----------------------------------------------------------
B=1000 #number of samples

for(k in 1:NROW(design)){
  print(paste0("Processing design no.: ",k))
  
  # get current design 
  I=design$I[k]; J=design$J[k]; M=design$M[k]
  betax=design$beta[k]
  
  # generate data according to the GRM-RT model
  Probs=Ydata=RTdata=DIFFdata=list()
  for(b in 1:B){
    P0 = array(NA,dim=c(I,J,M))
    DIFF = Y = RT = matrix(NA,I,J)
    for(i in 1:I){
      for(j in 1:J){
        pd = sum(mapply(function(m)exp(sum(Eta[[as.character(I)]][i]-Alpha[[as.character(M)]][[as.character(J)]][1:m,j])),1:M)) #normaliz factor for computing response probabilities
        P0[i,j,] = mapply(function(m)exp(sum(Eta[[as.character(I)]][i]-Alpha[[as.character(M)]][[as.character(J)]][1:m,j]))/pd,1:M) #probability for each response category
        DIFF[i,j] = sum(P0[i,j,]^2) #Diff index is closely related to the difficulty of responding
        Y[i,j] = which(rmultinom(n=1,size=1,prob=P0[i,j,])==1) #Sample observed responses
        RT[i,j] = 1*Gamma[[as.character(M)]][[as.character(J)]][j] + DIFF[i,j]*betax + Omega[[as.character(I)]][i] + rnorm(1,0,0.5) #Sample observed response times
      }
    }
    Probs[[b]]=P0; DIFFdata[[b]]=DIFF; Ydata[[b]]=Y; RTdata[[b]]=RT
  }
  gendata = list(n=I,J=J,M=M,beta=b,eta=Eta[[as.character(I)]],Alpha=Alpha[[as.character(M)]][[as.character(J)]],omega=Omega[[as.character(I)]],Gamma=Gamma[[as.character(M)]][[as.character(J)]],
                 Probs=Probs,DIFF=DIFFdata,Y=Ydata,RT=RTdata,design=design)  
  save(gendata,file=paste0("gendata/gendata_",k,".rda"))
}



