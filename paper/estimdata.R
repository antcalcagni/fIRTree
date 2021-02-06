
# Set environment ---------------------------------------------------------
rm(list=ls())
setwd("/home/antonio/MEGA/Lavoro_sync/Working projects/AC_fuzzy-rating/simulation_study/")
require(doParallel)
source("../utils.R")

ncores=4
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl,cores=ncores)


# Set functions and parameters --------------------------------------------
Tms=list()
Tx = matrix(c(0,NA,1,0,1,1),3,2,byrow=TRUE); rownames(Tx) = paste0(c(1:3));colnames(Tx) = paste0("node",1:NCOL(Tx)); Tms[[1]]=Tx
Tx = matrix(c(0,NA,NA,NA,1,0,NA,NA,1,1,0,NA,1,1,1,0,1,1,1,1),5,4,byrow=TRUE); rownames(Tx) = paste0(c(1:5));colnames(Tx) = paste0("node",1:NCOL(Tx));Tms[[2]]=Tx
names(Tms)=c("3","5"); rm(Tx)

computeQts = function(Y,RT,gendata,Tm){
  
  Ydata = irtrees::dendrify(Y,Tm)
  mod1_glm = glmmTMB::glmmTMB(formula = value ~ 0+item+(1|person), family = binomial, data = Ydata)
  
  #z_predicted = predict(mod1_glm,type = "response"); z_predicted[z_predicted<0.5]=0;z_predicted[z_predicted>=0.5]=1;
  #Ptab=table(Ydata$value,z_predicted); mod1_auc=sum(diag(Ptab))/sum(Ptab)
  
  Alpha_est = matrix(summary(mod1_glm)$coefficients$cond[,1],gendata$J,gendata$M-1)
  Eta_est = pracma::repmat(matrix(apply(lme4::ranef(mod1_glm)$cond$person,1,as.numeric)),1,gendata$M-1)
  
  fuzzyx = getFuzzyNumbers(Eta_est,Alpha_est,Tms[[as.character(gendata$M)]])
  C_est = matrix(fuzzyx$Yfuzzy[,3],gendata$n,gendata$J,byrow=TRUE)
  S_est = matrix(fuzzyx$Yfuzzy[,4],gendata$n,gendata$J,byrow=TRUE)
  Card_est = matrix(fuzzyx$Yfuzzy[,5],gendata$n,gendata$J,byrow=TRUE)
  
  Res=matrix(NA,5,gendata$J)
  px = matrix(0,gendata$n,1)
  for(j in 1:gendata$J){
    px[RT[,j]>quantile(RT[,j],c(0.5))]=1
    mod_px = glm(px~Card_est[,j],family=binomial)
    px_est=predict(mod_px,type="response");px_est[px_est<0.5]=0;px_est[px_est>=0.5]=1; 
    iid=which(is.na(Card_est[,j])); if(length(iid)>0){px=px[-iid]}; Ptab = table(px,px_est)
    Res[1,j] = sum(diag(Ptab))/sum(Ptab)
    
    px[RT[,j]>quantile(RT[,j],c(0.5))]=1
    mod_px = glm(px~sqrt(S_est[,j]),family=binomial)
    px_est=predict(mod_px,type="response");px_est[px_est<0.5]=0;px_est[px_est>=0.5]=1; 
    iid=which(is.na(S_est[,j])); if(length(iid)>0){px=px[-iid]}; Ptab = table(px,px_est)
    Res[2,j] = sum(diag(Ptab))/sum(Ptab)
    
    a0=AIC(lm(RT[,j]~C_est[,j])); a1=AIC(lm(RT[,j]~sqrt(S_est[,j])))
    Res[3,j] = which.min(c(a0,a1))-1
    
    Res[4,j] = summary(lm(RT[,j]~sqrt(S_est[,j])))$r.squared
    Res[5,j] = cor(RT[,j],sqrt(S_est[,j]))
  }
  return(Res)
}



# Compute MC quantities ---------------------------------------------------
B = length(gendata$Y)
chunk.size = floor(B/ncores)
n_resid = rep(0,ncores)
if((ncores*chunk.size)<B){n_resid[ncores] = B - ncores*chunk.size}
combAbind = abind::abind

funToExport=NULL
dataToExport=NULL
pkgsToExport=c("glmmTMB","pracma","irtrees","lme4")

for(k in 1:54){ #no. of scenarios
  cat(paste("\n Design cell n.: ",k),sep="")
  load(paste0("gendata/gendata_",k,".rda"))

  Yout <- foreach(h=1:ncores, .combine = "combAbind", .export = c(funToExport,dataToExport),.packages = pkgsToExport) %dopar% { 
    res = array(NA, dim = c(5,gendata$J,chunk.size+n_resid[h]),dimnames = list(c("auc_card","auc_s","aic_cs","r2_s","cor_rts"),NULL,NULL))
    iid_core = ((h-1)*chunk.size+1):((h*chunk.size)+n_resid[h])
    for(u in iid_core){
      res[,,u-(h-1)*chunk.size] = computeQts(Y=gendata$Y[[u]],RT=gendata$RT[[u]],gendata=gendata[1:3],Tm=Tms[[as.character(gendata$M)]])
    }
    res #it is needed to populate Yout along the third-dimension which contains the B number of samples
  }
  cat("\n > Saving current data")
  save(Yout,file=paste("estimdata_",k,".rda",sep=""))
}

doParallel::stopImplicitCluster(); parallel::stopCluster(cl)
cat("\n Done.")

