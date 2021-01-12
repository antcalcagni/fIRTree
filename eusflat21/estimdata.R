
# Set environment ---------------------------------------------------------
rm(list=ls()); graphics.off()
setwd("/home/eusflat21/")
source("utils.R")
require(doParallel)

ncores=4
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)


# Set functions and parameters --------------------------------------------
M = 5; N = M-1
Tm = matrix(c(0,NA,NA,NA,
              1,0,NA,NA,
              1,1,0,NA,
              1,1,1,0,
              1,1,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:NCOL(Tm)); print(Tm)

computeQts = function(Y,Tm,J,I,N){
  
  Ydata = irtrees::dendrify(as.matrix(Y),Tm)
  mod1_glm = glmmTMB::glmmTMB(formula = value ~ 0+item+(1|person), family = binomial, data = Ydata)
  
  Alpha_est = matrix(summary(mod1_glm)$coefficients$cond[,1],J,N)
  Eta_est = pracma::repmat(matrix(apply(lme4::ranef(mod1_glm)$cond$person,1,as.numeric)),1,N)
  
  fuzzyx = getFuzzyNumbers(Eta_est,Alpha_est,Tm)
  C_est = matrix(fuzzyx$Yfuzzy$C,I,J,byrow=TRUE)
  S_est = matrix(fuzzyx$Yfuzzy$S,I,J,byrow=TRUE)
  L_est = matrix(fuzzyx$Yfuzzy$L,I,J,byrow=TRUE)
  R_est = matrix(fuzzyx$Yfuzzy$R,I,J,byrow=TRUE)
  W_est = matrix(fuzzyx$Yfuzzy$W,I,J,byrow=TRUE)
  Kauff_est = matrix(fuzzyx$Yfuzzy[,"Kauff"],I,J,byrow=TRUE)
  
  Res = cbind(C_est,S_est,L_est,R_est,W_est,Kauff_est)
  return(Res)
}



# Compute MC quantities ---------------------------------------------------
B = 1000
chunk.size = floor(B/ncores)
n_resid = rep(0,ncores)
if((ncores*chunk.size)<B){n_resid[ncores] = B - ncores*chunk.size}
combAbind = abind::abind

funToExport=NULL
dataToExport=NULL
pkgsToExport=c("glmmTMB","pracma","irtrees","lme4")

for(k in 1:24){ #no. of scenarios
  cat(paste("\n Design cell n.: ",k),sep="")
  load(paste0("gendata/gendata_",k,".rda"))

  Yout <- foreach(h=1:ncores, .combine = "combAbind", .export = c(funToExport,dataToExport),.packages = pkgsToExport) %dopar% { 
    res = array(NA, dim = c(genout$n,genout$J*6,chunk.size+n_resid[h]),dimnames = list(NULL,paste0(rep(c("C","S","L","R","W","Kauff"),each=genout$J),1:genout$J),NULL))
    iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
    for(u in iid_core){
      res[,,u-(h-1)*chunk.size] = computeQts(Y = genout$gendata[[u]]$Yf,Tm = Tm,J = genout$J,I = genout$n,N = genout$M-1)
    }
    res #it is needed to populate Yout along the third-dimension which contains the B number of samples
  }
  cat("\n > Saving current data")
  save(Yout,file=paste("estimdata_",k,".rda",sep=""))
}

doParallel::stopImplicitCluster(); parallel::stopCluster(cl)
cat("\n Done.")

