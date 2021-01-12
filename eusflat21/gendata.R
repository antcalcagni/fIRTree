
# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
setwd("/home/eusflat21/")
source("utils.R")
set.seed(210108)


# Set design --------------------------------------------------------------
Ix = c(50,150,300) #low, medium, large sample sizes
Jx = c(10,20)
pix = c(0,0.25,0.5,0.75)
design = expand.grid(Ix,Jx,pix); colnames(design)=c("I","J","pix")
save(design,file = "design_all.rda")

M = 5; N = M-1
Tm = matrix(c(0,NA,NA,NA,
              1,0,NA,NA,
              1,1,0,NA,
              1,1,1,0,
              1,1,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:NCOL(Tm)); print(Tm)


# Set sim parameters ------------------------------------------------------
Eta=list()
for(k in 1:length(Ix)){
  Eta[[k]]=pracma::repmat(matrix(rnorm(n=Ix[k],0,1)),1,N)
}
names(Eta)=Ix

Alpha=list()
for(k in 1:length(Jx)){
  Alpha[[k]]=t(pracma::repmat(rnorm(Jx[k],-1.75,0.25),N,1))
}
names(Alpha)=Jx


# Simulate data -----------------------------------------------------------
B=1000 #number of samples

for(k in 1:NROW(design)){
  print(paste0("Processing design no.: ",k))
  
  I=design$I[k]; J=design$J[k];pix=design$pix[k]
  
  gendata = mapply(function(b){
    datax=genIRTree(Eta[[as.character(I)]],Alpha[[as.character(J)]],Tm)
    RepMat = sgr::replacement.matrix(Q=M,p=c(pix,0.0),fake.model="slight")
    Yf = sgr::rdatarepl(datax$Y,RepMat,printfp = FALSE)$Fx  
    return(list(Y=datax$Y,Yf=Yf,P0=datax$Py,fuzzydata=list(C=datax$C,S=datax$S,L=datax$L,R=datax$R,W=datax$W,Kauff=datax$Kauff)))
    },1:B,SIMPLIFY = FALSE)
  
  genout = list(n=I,J=J,M=M,pix=pix,eta=Eta[[as.character(I)]],Alpha=Alpha[[as.character(J)]],design=design,gendata=gendata)  
  save(genout,file=paste0("gendata/gendata_",k,".rda"))
}



