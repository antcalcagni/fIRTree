# Set environment ---------------------------------------------------------
rm(list=ls());graphics.off()
setwd("/home/paper/")
library(tikzDevice); library(xtable); library(kableExtra)
source("utils.R")


# Section 3 ---------------------------------------------------------------

## Figure 2
N=2;M=3
Tm = matrix(c(1,1,0,NA,1,0),M,N,byrow=TRUE)
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N)
PY = matrix(c(3.0,0.5,0.1,
              0.1,1.0,3.0,
              3.0,2.5,2.3,
              0.4,3,0.01),4,3,byrow=TRUE); PY=t(mapply(function(i)PY[i,]/sum(PY[i,]),1:NROW(PY)))
y = c(1,3,2,2)
cy = mapply(function(i)sum(PY[i,]*1:M),1:NROW(PY))
sy = 1/mapply(function(i)sum(PY[i,]*(1:M-cy[i])^2),1:NROW(PY))

#tikz('fig2b.tex',width=4.5,height=5)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(1:M/M,PY[i,],bty="n",xlim=c(0,1.25),ylim=c(-0.05,1),col=1,cex=2,xlab="",ylab="P(Y=m)"); segments(x0 = (1:M)/M,x1 = (1:M)/M,y0 = 0,y1=PY[i,],lty = 2,col=1)
  points(y[i]/M,0,col=2,cex=2); points(which.max(PY[i,])/M,0,col=4,cex=1.25)
  points(seq(0,1,length.out = 100),beta_fn(seq(0,1,length.out = 100),cy[i]/M,sy[i]),type="l",lty=2.5,col=1,lwd=1.5)
}
#dev.off()



# Section 4: Simulation Study ---------------------------------------------
load("simulation_data.rda")

## Table 1
X=rbind(cbind(dataout[dataout$M==3 & dataout$beta==-10.5,c("auc_s_mean","auc_s_sd")][1:3,],dataout[dataout$M==3 & dataout$beta==-10.5,c("auc_s_mean","auc_s_sd")][4:6,],
              dataout[dataout$M==3 & dataout$beta==-20.5,c("auc_s_mean","auc_s_sd")][1:3,],dataout[dataout$M==3 & dataout$beta==-20.5,c("auc_s_mean","auc_s_sd")][4:6,]),
        cbind(dataout[dataout$M==5 & dataout$beta==-10.5,c("auc_s_mean","auc_s_sd")][1:3,],dataout[dataout$M==5 & dataout$beta==-10.5,c("auc_s_mean","auc_s_sd")][4:6,],
              dataout[dataout$M==5 & dataout$beta==-20.5,c("auc_s_mean","auc_s_sd")][1:3,],dataout[dataout$M==5 & dataout$beta==-20.5,c("auc_s_mean","auc_s_sd")][4:6,]))

Xtab = cbind(matrix(paste0(round(X[,1],3)," (",round(X[,2],3),")"),ncol = 1),
             matrix(paste0(round(X[,3],3)," (",round(X[,4],3),")"),ncol = 1),
             matrix(paste0(round(X[,5],3)," (",round(X[,6],3),")"),ncol = 1),
             matrix(paste0(round(X[,7],3)," (",round(X[,8],3),")"),ncol = 1))
rownames(Xtab) = rep(paste0("$I=",c(50,150,500),"$"),2)
colnames(Xtab) = rep(paste0("$J=",c(5,15),"$"),2)
Xtab
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab:sim_1"
attributes(Xtab_tex)$align = rep("r",length(attributes(Xtab_tex)$align))
attributes(Xtab_tex)$row.names = rep(paste0("$I=",c(50,150,500),"$"),2)
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


# Section 5: Case study 1 -------------------------------------------------
load("dataFaking.Rdata")
M = max(X_CG[,1]); J = NCOL(X_CG)

#Load data (CG: honest; FMG: faking)
X_CG_prop = apply(X_CG,2,function(x)tablex(x,M))/NROW(X_CG) #honest
X_FMG_prop = apply(X_FMG,2,function(x)tablex(x,M))/NROW(X_FMG) #faking

## Table 2
X_prop = mapply(function(i)rbind(X_CG_prop[i,],X_FMG_prop[i,]),1:NROW(X_CG_prop))
X_prop = cbind(X_prop,mapply(function(j){sum(1:M*X_prop[j,])},1:NROW(X_prop)))
colnames(X_prop) = c(paste0("$Y=",1:M,"$"),"mean response")
rownames(X_prop) = paste0("item",c("1 (H)","1 (F)","4 (H)","4 (F)","7 (H)","7 (F)","8 (H)","8 (F)"))
Xtab_tex = xtable::xtable(X_prop)
attributes(Xtab_tex)$caption = "Case study 1: Observed frequency tables a function of item number and type of group (H: honest group; F: faking group)."
attributes(Xtab_tex)$label = "tab:cs1_1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


#Analysing CG condition (honest)
X=X_CG
I=NROW(X); N=M-2
Tm = matrix(c(1,0,0,1,0,1,0,NA,NA,1,1,0,1,1,1),M,N,byrow=TRUE) #Mapping matrix: response categories (rows) x nodes (cols)
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N)
print(Tm)

Ydata = irtrees::dendrify(as.matrix(X),Tm)
mod1_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Ydata)
summary(mod1_glm)

alphas1 = matrix(summary(mod1_glm)$coefficients$cond[,1],J,N) 
sigma_eta1 = summary(mod1_glm)$varcor$cond$person 
R_eta1 = attr(x = sigma_eta1,which = 'correlation') 
sigma_eta1 = matrix(as.numeric(sigma_eta1),N,N) 
etas1 = ranef(mod1_glm)$cond$person

z_predicted = predict(mod1_glm,type = "response"); z_predicted[z_predicted<0.5]=0;z_predicted[z_predicted>=0.5]=1;
X = rbind(table(Ydata$value),table(z_predicted)); rownames(X)=c("obs","fit");print(X)
roc1 = pROC::roc(Ydata$value,z_predicted); pROC::auc(roc1)

## Table 3 (honest part)
Xtab = summary(mod1_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:4,],Xtab[5:8,],Xtab[9:12,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),3)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 1: Estimates of item parameters"
attributes(Xtab_tex)$label = "tab:cs1_2"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table 4 (honest part)
Xtab = R_eta1; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(sigma_eta1)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 1: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs1_3"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


#Analysing FMG condition (faking)
X=X_FMG
I=NROW(X); N=M-2
Tm = matrix(c(1,0,0,1,0,1,0,NA,NA,1,1,0,1,1,1),M,N,byrow=TRUE) #Mapping matrix: response categories (rows) x nodes (cols)
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N)
print(Tm)

Ydata = irtrees::dendrify(as.matrix(X),Tm)
mod2_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Ydata)
summary(mod2_glm)

alphas2 = matrix(summary(mod2_glm)$coefficients$cond[,1],J,N)
sigma_eta2 = summary(mod2_glm)$varcor$cond$person 
R_eta2 = attr(x = sigma_eta2,which = 'correlation')
sigma_eta2 = matrix(as.numeric(sigma_eta2),N,N) 
etas2 = ranef(mod2_glm)$cond$person

z_predicted = predict(mod2_glm,type = "response"); z_predicted[z_predicted<0.5]=0;z_predicted[z_predicted>=0.5]=1;
X = rbind(table(Ydata$value),table(z_predicted)); rownames(X)=c("obs","fit");print(X)
roc2 = pROC::roc(Ydata$value,z_predicted); pROC::auc(roc2)

## Table 3 (faking part)
Xtab = summary(mod2_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:4,],Xtab[5:8,],Xtab[9:12,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),3)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 1: Estimates ($\\hat\\theta$) and standard errors ($\\sigma_{\\hat{\\theta}}$) for item parameters (faking model)"
attributes(Xtab_tex)$label = "tab:cs1_4"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table 4 (faking part)
Xtab = R_eta2; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(sigma_eta2)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 1: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs1_5"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

#Get fuzzy numbers from estimated IRTrees
fuzzyx = getFuzzyNumbers(Eta = etas1,Alpha = alphas1,Tm = Tm) #honest
fuzzyf = getFuzzyNumbers(Eta = etas2,Alpha = alphas2,Tm = Tm) #faking


X=X_CG_prop; X2=X_FMG_prop
fuzzyf = list(fuzzyx,fuzzyf) #single list for both honest and faking


## Figure 3
D = fuzzyf[[1]]$Yfuzzy; D_summary = summary(D$precision)[c(1,3,5)]
iid = mapply(function(k)which.max(D$precision[D$precision<=D_summary[k]]),k=1:length(D_summary)); iid = D[iid,]
D = fuzzyf[[2]]$Yfuzzy; D_summary = summary(D$precision)[c(1,3,5)]; 
iid2 = mapply(function(k)which.max(D$precision[D$precision<=D_summary[k]]),k=1:NROW(D_summary)); iid2 = D[iid2,]
iid = rbind(iid,iid2)

tikz(file='cs2.tex',width=9.5,height=4.5)
par(mfrow=c(2,3),mar=rep(1.5,4),mai=c(0.3, 0.6, 0.3, 0.15));ylabx=""
for(u in 1:NROW(iid)){
  iid[u,]
  i=iid$sbj[u]; j=iid$itm[u]; k=as.numeric(rownames(iid[u,]))
  py=as.numeric(fuzzyf[[1]]$PY[k,1:M+2]); mi=iid$mode[u]; si=iid$precision[u]
  xtitle = paste0("(",LETTERS[u],") ","Rater: ",i,", Item: ",j,", y = ",X_CG[i,j],", m = ",round(mi,2)*M)
  xsup=seq(0,1,length.out = 100);fy=beta_fn(xsup,mi,si); plot(xsup*M,fy,bty="n",ylab=ylabx,xlab="",lwd=1.25,type="l",main=xtitle,adj = 0.01)
  points(1:M,py,col=1,cex=2); segments(x0 = (1:M),x1 = (1:M),y0 = 0,y1=py,lty = 2,col=1,lwd=1.25)
  points(X_CG[i,j],0,col=2,cex=2,lwd=2); points(mi*M,0,lty=2,col="blue",cex=1.8,lwd=2,pch=20)
}
dev.off()


## Figure 4
cols=c("peachpuff4","peru")
ylab=c("P(Y=m)","mode (m)","precision (s)","cardinality","centroid")
jit=seq(0.15,0.45,length.out = 3)
xnms=rep("",2)

#tikz(file='cs1.tex',width=8.5,height=8.5)
par(mfcol=c(5,4),mar=rep(1.5,4),mai=c(0.25, 0.6, 0.15, 0.15))
for(j in 1:J){
  if(j==1){ylabx=ylab}else{ylabx=rep("",5)}
  
  plot(1:M,X[1:M,j],bty="n",xlim=c(0,M+1),ylim=c(-0.05,1),col=cols[1],cex=2,xlab="",ylab=ylabx[1],main = paste0("item ",j)); segments(x0 = (1:M),x1 = (1:M),y0 = 0,y1=X[1:M,j],lty = 1,col=cols[1],lwd=1.5)
  points((1:M)+jit[1],X2[1:M,j],col=cols[2],cex=2); segments(x0 = (1:M)+jit[1],x1 = (1:M)+jit[1],y0 = 0,y1=X2[1:M,j],lty = 1,col=cols[2],lwd=1.5)
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$mode[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$mode[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[2],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$mode,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$mode,1:2)))))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$precision[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$precision[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[3],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$precision,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$precision,1:2)))))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$cardinality[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$cardinality[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[4],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$cardinality,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$cardinality,1:2)))))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$centroid[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$centroid[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=c("honest","faking"),col=cols,ylab=ylabx[5],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$centroid,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$centroid,1:2)))))
}
#dev.off()



# Section 5: Case study 2 -------------------------------------------------
load("moral_dilemmas.Rdata")
I=NROW(X); J=4; M=5; N=4
X_foot_prop = apply(X[,4:7],2,function(x)tablex(x,M))/I 
X_trans_prop = apply(X[,8:11],2,function(x)tablex(x,M))/I 

#Sequential IRTree for both scenarios
Tm = matrix(c(0,NA,NA,NA,
              1,0,NA,NA,
              1,1,0,NA,
              1,1,1,0,
              1,1,1,1),M,N,byrow=TRUE) #Mapping matrix: response categories (rows) x nodes (cols)
rownames(Tm) = paste0(c(1:M));colnames(Tm) = paste0("node",1:N)
print(Tm)

#Analysing Footbridge case
Ydata = irtrees::dendrify(as.matrix(X[,4:7]),Tm)
mod1_glm = glmmTMB::glmmTMB(formula = value ~0+item:node+(0+node|person), family = binomial, data = Ydata)
summary(mod1_glm)

alphas1 = matrix(summary(mod1_glm)$coefficients$cond[,1],J,N) 
sigma_eta1 = summary(mod1_glm)$varcor$cond$person
R_eta1 = attr(x = sigma_eta1,which = 'correlation')
sigma_eta1 = matrix(as.numeric(sigma_eta1),N,N)
etas1 = ranef(mod1_glm)$cond$person

z_predicted = predict(mod1_glm,type = "response"); z_predicted[z_predicted<0.5]=0;z_predicted[z_predicted>=0.5]=1;
roc1 = pROC::roc(Ydata$value,z_predicted); pROC::auc(roc1)

## Table 5 (footbridge part)
Xtab = summary(mod1_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:4,],Xtab[5:8,],Xtab[9:12,],Xtab[13:16,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),4)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 2: Estimates of item parameters"
attributes(Xtab_tex)$label = "tab:cs2_1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table 6 (footbridge part)
Xtab = R_eta1; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(sigma_eta1)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 2: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs1_3"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

#Analysing Transplant case
Ydata = irtrees::dendrify(as.matrix(X[,8:11]),Tm)
mod2_glm = glmmTMB::glmmTMB(formula = value ~0+item:node+(0+node|person), family = binomial, data = Ydata)
summary(mod2_glm)

alphas2 = matrix(summary(mod2_glm)$coefficients$cond[,1],J,N)
sigma_eta2 = summary(mod2_glm)$varcor$cond$person
R_eta2 = attr(x = sigma_eta2,which = 'correlation')
sigma_eta2 = matrix(as.numeric(sigma_eta2),N,N)
etas2 = ranef(mod2_glm)$cond$person

z_predicted = predict(mod2_glm,type = "response"); z_predicted[z_predicted<0.5]=0;z_predicted[z_predicted>=0.5]=1;
roc2 = pROC::roc(Ydata$value,z_predicted); pROC::auc(roc2)

## Table 5 (transplant part)
Xtab = summary(mod2_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:4,],Xtab[5:8,],Xtab[9:12,],Xtab[13:16,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),4)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 2: Estimates of item parameters"
attributes(Xtab_tex)$label = "tab:cs2_1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table 6 (transplant part)
Xtab = R_eta2; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(sigma_eta2)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study 2: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs1_3"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

#Get fuzzy numbers from estimated IRTrees
fuzzyx = getFuzzyNumbers(Eta = etas1,Alpha = alphas1,Tm = Tm) #footbridge
fuzzyf = getFuzzyNumbers(Eta = etas2,Alpha = alphas2,Tm = Tm) #transplant

X0=X; X=X_foot_prop; X2=X_trans_prop; fuzzyf = list(fuzzyx,fuzzyf)

## Figure 5
cols=c("peachpuff4","peru")
ylab=c("P(Y=m)","mode (m)","precision (s)","cardinality","centroid")
jit=seq(0.15,0.45,length.out = 3)
xnms=rep("",2)

#tikz(file='cs4.tex',width=8.5,height=8.5)
par(mfcol=c(5,4),mar=rep(1.5,4),mai=c(0.25, 0.6, 0.15, 0.15))
for(j in 1:J){
  if(j==1){ylabx=ylab}else{ylabx=rep("",5)}
  
  plot(1:M,X[1:M,j],bty="n",xlim=c(0,M+1),ylim=c(-0.05,0.5),col=cols[1],cex=2,xlab="",ylab=ylabx[1],main = paste0("item ",j)); segments(x0 = (1:M),x1 = (1:M),y0 = 0,y1=X[1:M,j],lty = 1,col=cols[1],lwd=1.5)
  points((1:M)+jit[1],X2[1:M,j],col=cols[2],cex=2); segments(x0 = (1:M)+jit[1],x1 = (1:M)+jit[1],y0 = 0,y1=X2[1:M,j],lty = 1,col=cols[2],lwd=1.5)
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$mode[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$mode[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[2],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$mode,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$mode,1:2)))))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$precision[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$precision[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[3],ylim=c(0,250))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$cardinality[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$cardinality[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=xnms,col=cols,ylab=ylabx[4],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$cardinality,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$cardinality,1:2)))))
  
  D = cbind(c(fuzzyf[[1]]$Yfuzzy$centroid[fuzzyx$Yfuzzy$itm==j],fuzzyf[[2]]$Yfuzzy$centroid[fuzzyx$Yfuzzy$itm==j]),rep(1:2,each=I)); D = as.data.frame(D); colnames(D)=c("xvalue","cond"); D$cond = as.factor(D$cond)
  boxplot(xvalue~cond,data=D,frame=FALSE,names=c("footbridge","transplant"),col=cols,ylab=ylabx[5],ylim=c(min(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$centroid,1:2))),max(unlist(mapply(function(k)fuzzyf[[k]]$Yfuzzy$centroid,1:2)))))
}
#dev.off()


## Figure 6
sbjs = c(9,20,23,483)
xsup=seq(0,1,length.out = 100)
#tikz(file='cs2_2.tex',width=10.5,height=6.5)
par(mfrow=c(4,4),mar=rep(1.5,4),mai=c(0.5, 0.6, 0.25, 0.15))
h=0
for(u in 1:length(sbjs)){
  i=sbjs[u]
  x1=fuzzyf[[1]]$Yfuzzy[fuzzyf[[1]]$Yfuzzy$sbj==i,3:4]
  x2=fuzzyf[[2]]$Yfuzzy[fuzzyf[[2]]$Yfuzzy$sbj==i,3:4]
  
  for(j in 1:J){
    h=h+1
    if(j==1){ylabx=paste0("Rater: ",i)}else{ylabx=""}
    xtitle = paste0("(",LETTERS[h],") ","Item: ",j,", yF = ",X0[i,j+3],", yT = ",X0[i,j+7])
    fy=beta_fn(xsup,x1[j,1],x1[j,2]);plot(xsup*M+1,fy,bty="n",ylab=ylabx,xlab="",lwd=1.75,type="l",main=xtitle,col=cols[1],adj=0.05)
    fy=beta_fn(xsup,x2[j,1],x2[j,2]);points(xsup*M+1,fy,bty="n",xlab="",lwd=1.75,type="l",col=cols[2])
    points(c(round(x1[j,1],2)*M+1,round(x2[j,1],2)*M+1),y=c(0,0),col=cols,cex=1.8,lwd=2,pch=20)
    abline(v = c(round(x1[j,1],2)*M+1,round(x2[j,1],2)*M+1),lty=2,col=cols,lwd=1.5)
  }
}  
add_legend("bottom",fill = cols,legend = c("footbridge","transplant"),border = FALSE,bty = "n",ncol = 2)
#dev.off()

