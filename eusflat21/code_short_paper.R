# Set environment ---------------------------------------------------------
rm(list=ls());graphics.off()
setwd("/home/eusflat21/")
source("utils.R")


# Simulation study --------------------------------------------------------
load("design_all.rda")

## Table 1
load("simulation_data_2.rda")
X = rbind(cbind(dataout[dataout$J==10,c(4,7)],dataout[dataout$J==10,c(5,8)],dataout[dataout$J==10,c(6,9)]),
          cbind(dataout[dataout$J==20,c(4,7)],dataout[dataout$J==20,c(5,8)],dataout[dataout$J==20,c(6,9)]))

Xtab = cbind(matrix(paste0(round(X[,1],3)," (",round(X[,2],3),")"),ncol = 1),
             matrix(paste0(round(X[,3],3)," (",round(X[,4],3),")"),ncol = 1),
             matrix(paste0(round(X[,5],3)," (",round(X[,6],3),")"),ncol = 1))
rownames(Xtab) = rep(paste0("$I=",c(50,150,500),"$"),2)
colnames(Xtab) = c("$C$","$R-L$","$W$")

Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab:sim_2"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
attributes(Xtab_tex)$row.names = rep(paste0("$I=",c(50,150,500),"$"),2)
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


## Table 2
load("simulation_data_1.rda")
X = rbind(cbind(dataout[dataout$J==10 & dataout$pix==0.0,c(5,7)],dataout[dataout$J==10 & dataout$pix==0.25,c(5,7)],dataout[dataout$J==10 & dataout$pix==0.50,c(5,7)],dataout[dataout$J==10 & dataout$pix==0.75,c(5,7)]),
          cbind(dataout[dataout$J==20 & dataout$pix==0.0,c(5,7)],dataout[dataout$J==20 & dataout$pix==0.25,c(5,7)],dataout[dataout$J==20 & dataout$pix==0.50,c(5,7)],dataout[dataout$J==20 & dataout$pix==0.75,c(5,7)]))

Xtab = cbind(matrix(paste0(round(X[,1],3)," (",round(X[,2],3),")"),ncol = 1),
             matrix(paste0(round(X[,3],3)," (",round(X[,4],3),")"),ncol = 1),
             matrix(paste0(round(X[,5],3)," (",round(X[,6],3),")"),ncol = 1),
             matrix(paste0(round(X[,7],3)," (",round(X[,8],3),")"),ncol = 1))
rownames(Xtab) = rep(paste0("$I=",c(50,150,500),"$"),2)
colnames(Xtab) = paste0("$\\pi=",unique(design$pix),"$")

Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Simulation study: "
attributes(Xtab_tex)$label = "tab:sim_1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
attributes(Xtab_tex)$row.names = rep(paste0("$I=",c(50,150,500),"$"),2)
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})




# Case study --------------------------------------------------------------
source("utils.R")
load("pns.rda")

M=6;N=M-2; I=NROW(X); J=NCOL(X)
Tm = matrix(c(1,NA,0,0,
              1,NA,0,1,
              0,0,NA,NA,
              0,1,NA,NA,
              1,NA,1,0,
              1,NA,1,1),M,N,byrow=TRUE); rownames(Tm) = paste0(c(1:M));colnames(Tm) = c("M","Aw","As","E"); print(Tm)

Ydata = irtrees::dendrify(as.matrix(X),Tm)
mod1_glm = glmmTMB::glmmTMB(formula = value ~ 0+item:node+(0+node|person), family = binomial, data = Ydata)
summary(mod1_glm)

Alpha_est = matrix(summary(mod1_glm)$coefficients$cond[,1],J,N) #easiness or difficulty?
Sigma_eta = summary(mod1_glm)$varcor$cond$person #covariance for nodes
R_eta = attr(x = Sigma_eta,which = 'correlation') #correlation for nodes
Sigma_eta = matrix(as.numeric(Sigma_eta),N,N) #removing attributes from the previous object
Eta_est = ranef(mod1_glm)$cond$person

## Table 3
Xtab = summary(mod1_glm)$coefficients$cond[,1:2]
Xtab = cbind(Xtab[1:5,],Xtab[6:10,],Xtab[11:15,],Xtab[16:20,])
rownames(Xtab) = paste0("$\\alpha_{",1:J,"}$")
colnames(Xtab) = rep(c("$\\hat{\\theta}$","$\\sigma_{\\hat{\\theta}}$"),4)
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study: Estimates of item parameters"
attributes(Xtab_tex)$label = "tab:cs1"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

## Table 4
Xtab = R_eta; Xtab[upper.tri(Xtab)]=NA
Xtab = cbind(Xtab,sqrt(diag(Sigma_eta)))
rownames(Xtab) = paste0("$\\eta_{",1:N,"}$")
colnames(Xtab) = c(paste0("$\\eta_{",1:N,"}$"),"$\\hat\\sigma_\\eta$")
Xtab_tex = xtable::xtable(Xtab)
attributes(Xtab_tex)$caption = "Case study: Estimated correlation matrix and standard deviations ($\\hat{\\sigma}_\\eta$) for the latent traits"
attributes(Xtab_tex)$label = "tab:cs2"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


fuzzyx = getFuzzyNumbers(Eta_est,Alpha_est,Tm)
C_est = matrix(fuzzyx$Yfuzzy$C,I,J,byrow = TRUE)
L_est = matrix(fuzzyx$Yfuzzy$L,I,J,byrow = TRUE)
R_est = matrix(fuzzyx$Yfuzzy$R,I,J,byrow = TRUE)
W_est = matrix(fuzzyx$Yfuzzy$W,I,J,byrow = TRUE)
Kauff_est = matrix(fuzzyx$Yfuzzy$Kauff,I,J,byrow = TRUE)

## Figure 3a
tikzDevice::tikz(file='fig3a.tex',width=4.5,height=4)
i=98; j=1
dbnx = dombi_fn(L_est[i,j],C_est[i,j],R_est[i,j],W_est[i,j])
pu = as.numeric(fuzzyx$PY[fuzzyx$PY$sbj==i & fuzzyx$PY$itm==j,3:(M+2)])
plot(1:M,pu,bty="n",xlim=c(0,(M+1)),ylim=c(-0.05,1),col=1,cex=2,xlab="",ylab=""); segments(x0 = (1:M),x1 = (1:M),y0 = 0,y1=pu,lty = 2,col=1)
points(dbnx$x*(M+1),dbnx$mu,type="l",xlim=c(0,M+1),ylim=c(0,1),bty="n",xlab="",ylab="",col=1,lwd=1)
points(X[i,j],0,col=2,pch=20,lwd=4)
dev.off()

## figure 3b
tikzDevice::tikz(file='fig3b.tex',width=4.5,height=4)
i=4; j=4
dbnx = dombi_fn(L_est[i,j],C_est[i,j],R_est[i,j],W_est[i,j])
pu = as.numeric(fuzzyx$PY[fuzzyx$PY$sbj==i & fuzzyx$PY$itm==j,3:(M+2)])
plot(1:M,pu,bty="n",xlim=c(0,(M+1)),ylim=c(-0.05,1),col=1,cex=2,xlab="",ylab=""); segments(x0 = (1:M),x1 = (1:M),y0 = 0,y1=pu,lty = 2,col=1,lwd=1)
points(dbnx$x*(M+1),dbnx$mu,type="l",xlim=c(0,M+1),ylim=c(0,1),bty="n",xlab="",ylab="",col=1,lwd=1)
points(X[i,j],0,col=2,pch=20,lwd=4)
dev.off()

## Figure 4
tikzDevice::tikz(file='fig4.tex',width=4.5,height=4)
cols=c("lightsalmon","lightskyblue4","mediumorchid3")
sbjs=c(4,28,321)
X[sbjs,j]

i=sbjs[1];j=2
dbnx = dombi_fn(L_est[i,j],C_est[i,j],R_est[i,j],W_est[i,j])
plot(dbnx$x*(M+1),dbnx$mu,type="l",xlim=c(0,M+1),ylim=c(0,1),bty="n",xlab="",ylab="",col=cols[1],lwd=2.5)

i=sbjs[2];
dbnx = dombi_fn(L_est[i,j],C_est[i,j],R_est[i,j],W_est[i,j])
points(dbnx$x*(M+1),dbnx$mu,type="l",bty="n",xlab="",ylab="",col=cols[2],lwd=2.5)

i=sbjs[3];
dbnx = dombi_fn(L_est[i,j],C_est[i,j],R_est[i,j],W_est[i,j])
points(dbnx$x*(M+1),dbnx$mu,type="l",bty="n",xlab="",ylab="",col=cols[3],lwd=2.5)
add_legend("bottom",fill = cols,legend = paste0("i=",sbjs),border = FALSE,bty = "n",ncol = 4)
dev.off()







