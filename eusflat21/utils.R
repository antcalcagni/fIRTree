
beta_fn = function(x,mi,gi){
  a = 1+gi*mi; b= 1+gi*(1-mi)
  C = (((a-1)/(a+b-2))^(a-1) * (1-((a-1)/(a+b-2)))^(b-1))
  fx = 1/C * x^(a-1) * (1-x)^(b-1)
  return(fx)
}

dombi_fn = function(l,m,r,w=1){
  x=seq(from=l,to = r,length.out = 101)
  mux = rep(0,length(x))
  mux[x<=m&x>l]= 1/(1+((m-x[x<=m&x>l])/(x[x<=m&x>l]-l))^w)
  mux[x>m&x<r]= 1/(1+((r-x[x>m&x<r])/(x[x>m&x<r]-m))^-w)
  return(list(x=x,mu=mux))
}

genIRTree = function(Eta=NULL,Alpha=NULL,Tm=NULL){
  I=NROW(Eta);J=NROW(Alpha);N=NCOL(Tm);M=NROW(Tm)
  Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}
  Y=C=S=L=R=W=Kauff=matrix(NA,I,J)
  P0=array(NA,dim=c(I,J,M))
  xsup=seq(0,1,length.out = M)
  for(i in 1:I){
    for(j in 1:J){
      Py_ij = mapply(function(m){prod((exp((Eta[i,]+Alpha[j,])*Tm[m,])/(1+exp(Eta[i,]+Alpha[j,])))^Dm[m,])},1:M)
      C[i,j] = sum(xsup*Py_ij) #betafn
      S[i,j] = 1/sum(Py_ij*(xsup-C[i,j])^2) #betafn
      xtg = triangular_fn(C[i,j],S[i,j])$par #trgfn
      L[i,j] = xtg[1] #dbnfn
      R[i,j] = xtg[2] #dbfn
      W[i,j] = sum(Py_ij^2) #dbfn
      Kauff[i,j] = kaufmann_index(dombi_fn(l = xtg[1],r = xtg[2],m = xtg[3],w = W[i,j])$mu)
      Y[i,j] = which(rmultinom(n=1,size=1,prob=Py_ij)==1)
      P0[i,j,] = Py_ij
    }
  }
  colnames(Y) = colnames(C) = colnames(S) = colnames(L) = colnames(R) = colnames(W) = colnames(Kauff) = paste0("V",1:J)
  return(list(Y=Y,C=C,S=S,L=L,R=R,W=W,Kauff=Kauff,Py=P0))
}

getFuzzyNumbers =function(Eta=NULL,Alpha=NULL,Tm=NULL){
  Alpha=as.matrix(Alpha); Eta=as.matrix(Eta)
  J=NROW(Alpha);M=NROW(Tm);I=NROW(Eta);N=NCOL(Eta)
  Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}
  PY = matrix(NA,I*J,M) #JIxM (J nested within I)
  ETA = kronecker(Eta,matrix(1,J,1)) #JIxN (rep each Eta J times)
  ALPHA = kronecker(matrix(1,I,1),Alpha) #JIxN (rep each Alpha I times)
  for(m in 1:M){
    TM = matrix(1,I*J,1)%*%Tm[m,] 
    DM = matrix(1,I*J,1)%*%Dm[m,]
    PY[,m] = apply( (exp((ETA+ALPHA)*TM)/(matrix(1,J*I,N)+exp((ETA+ALPHA))))^DM, 1, prod)
  }
  xsup=seq(0,1,length.out = M)
  YIJ_stats=matrix(NA,I*J,6); colnames(YIJ_stats) = c("C","S","L","R","W","Kauff")
  YIJ_stats[,1] = mapply(function(i)sum(PY[i,]*xsup),1:(I*J))
  YIJ_stats[,2] = mapply(function(i)1/sum(PY[i,]*(xsup-YIJ_stats[i,1])^2),1:(I*J))
  YIJ_stats[,3:4] = t(mapply(function(i){triangular_fn(mi = YIJ_stats[i,1],si = YIJ_stats[i,2])$par},1:(I*J)))[,1:2]
  YIJ_stats[,5] = mapply(function(i){sum(PY[i,]^2)},1:(I*J))
  YIJ_stats[,6] = mapply(function(i){kaufmann_index(dombi_fn(l = YIJ_stats[i,3],r = YIJ_stats[i,4],m = YIJ_stats[i,1],w = YIJ_stats[i,5])$mu)},1:(I*J))
  YIJ_stats = data.frame(rep(1:I,each=J),rep(1:J,I),YIJ_stats); colnames(YIJ_stats)[1:2]=c("sbj","itm")
  PY = data.frame(rep(1:I,each=J),rep(1:J,I),PY); colnames(PY)=c("sbj","itm",paste0("m",1:M))
  return(list(PY=PY,Yfuzzy=YIJ_stats))
}

kaufmann_index = function(mux){
  return(2*sum(abs(mux-(mux>=0.5)))/length(mux))
}

tablex = function(x,M){
  return(mapply(function(m)sum(x==m),1:M))
}

triangular_fn = function(mi,si,plotx=FALSE){
  #mi: mode 
  #si: precision
  #Source: Williams, T. M. (1992). Practical use of distributions in network analysis. Journal of the Operational Research Society, 43(3), 265-270.
  lb=0;ub=1 #standard beta function
  mux = (1+mi*si)/((1+mi*si)+(1+si*(1-mi)))
  vx = 1/si; w=3.5
  cx = sqrt(w*vx - 3*(mi-mux))
  bx =  (cx + 3*mi - 3*mux)/2
  ax = mi-bx; 
  ai=ifelse(ax<lb,lb,ax)
  bi=ax+cx; bi=ifelse(bi>ub,ub,bi)
  mi=ax+bx
  xsup=seq(lb,ub,length.out=100)
  fy=extraDistr::dtriang(xsup,ai,bi,mi); fy=fy/max(fy); 
  if(plotx){plot(xsup,fy,bty="n",type="l")}
  return(list(xsup=xsup,mux=fy,par=c(ai,bi,mi)))
}

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
