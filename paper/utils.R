
genIRTree = function(Eta=NULL,Alpha=NULL,Tm=NULL){
  I=NROW(Eta);J=NROW(Alpha);N=NCOL(Tm);M=NROW(Tm)
  Dm=matrix(1,M,N); for(n in 1:N){Dm[is.na(Tm[,n]),n]=0}
  Y=C=S=matrix(NA,I,J)
  for(i in 1:I){
    for(j in 1:J){
      Py_ij = mapply(function(m){prod((exp((Eta[i,]+Alpha[j,])*Tm[m,])/(1+exp(Eta[i,]+Alpha[j,])))^Dm[m,])},1:M)
      C[i,j] = sum(1:M*Py_ij) #betafn
      S[i,j] = 1/sum(Py_ij*(1:M-C[i,j])^2) #betafn
      
      Y[i,j] = which(rmultinom(n=1,size=1,prob=Py_ij)==1)
    }
  }
  colnames(Y) = colnames(C) = colnames(S) = paste0("V",1:J)
  return(list(Y=Y,C=C,S=S))
}


getFuzzyNumbers = function(Eta=NULL,Alpha=NULL,Tm=NULL,MC=FALSE,K=5000){
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
  if(MC==TRUE){
    YIJ = t(mapply(function(i)apply( rmultinom(K,1,PY[i,]), 2, function(x)which(x==1) ),1:(I*J))) #JIxK (each row contains K replicates for the couple {i,j})
    YIJ_stats = t(apply(YIJ,1,function(x)computeFuzzyQuantities(x,M))); colnames(YIJ_stats) = c("mode","precision","cardinality","centroid")
    YIJ = data.frame(rep(1:I,each=J),rep(1:J,I),YIJ); colnames(YIJ)=c("sbj","itm",paste0("sample",1:K))
  }else{
    xsup=seq(0,1,length.out = M)
    YIJ=NULL
    YIJ_stats=matrix(NA,I*J,4); colnames(YIJ_stats) = c("mode","precision","cardinality","centroid")
    YIJ_stats[,1] = mapply(function(i)sum(PY[i,]*xsup),1:(I*J))
    YIJ_stats[,2] = mapply(function(i)1/sum(PY[i,]*(xsup-YIJ_stats[i,1])^2),1:(I*J))
    YIJ_stats[,3] = mapply(function(i){tryCatch({fx=integrate(function(x){beta_fn(x,YIJ_stats[i,1],YIJ_stats[i,2])},0,1)$value},error = function(error_condition){fx=NA})},1:(I*J))
    YIJ_stats[,4] = mapply(function(i)defuzzify_centroid(YIJ_stats[i,1],YIJ_stats[i,2],what = "mean")*M+1,1:(I*J))
  }
  YIJ_stats = data.frame(rep(1:I,each=J),rep(1:J,I),YIJ_stats); colnames(YIJ_stats)[1:2]=c("sbj","itm")
  PY = data.frame(rep(1:I,each=J),rep(1:J,I),PY); colnames(PY)=c("sbj","itm",paste0("m",1:M))
  return(list(PY=PY,Ydata=YIJ,Yfuzzy=YIJ_stats))
}

computeFuzzyQuantities = function(x,M){
  #x is a Kx1 vector containing replicates for the couple {i,j}
  #M is the number of response categories 
  xp = tablex(x,M)/length(x) #probs
  xm = sum(1:M*xp) #expected value of multiverse (mode of the fuzzy number)
  xs = 1/sum(xp*(1:M-xm)^2) #precision of multiverse (precision of the fuzzy number)
  fc = integrate(function(x){beta_fn(x,xm/M,xs)},0,1)$value #cardinality of the fuzzy number
  fm = defuzzify_centroid(xm/M,xs,what = "mean")*M #centroid of the fuzzy number
  return(c(xm,xs,fc,fm))
}


beta_fn = function(x,mi,gi){
  a = 1+gi*mi; b= 1+gi*(1-mi)
  C = (((a-1)/(a+b-2))^(a-1) * (1-((a-1)/(a+b-2)))^(b-1))
  fx = 1/C * x^(a-1) * (1-x)^(b-1)
  return(fx)
}

triangular_fn = function(mi,si,plotx=FALSE){
  #mi: mode of fuzzy beta number
  #si: precision of fuzzy beta number
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
  return(list(mux=fy,par=c(ai,bi,mi)))
}


fuzziness = function(m=0.5,s=1){
  #See: Chakravarty, S. R., & Roy, T. (1985). Measurement of fuzziness: a general approach. Theory and decision, 19(2), 163-169.
  muy = beta_fn(seq(1e-04,1-1e-04,length.out = 100),m,s);
  fy = (1/length(muy))*(sum(muy)*log10(sum(muy))-sum(muy*log10(muy)))^0.5 * ((sum(1-muy)*log10(sum(1-muy)))-sum((1-muy)*log10(1-muy)))^0.5
  return(fy)
}

defuzzify_centroid = function(mi,si,what=c("mean","mode")){
  a1 = 1+si*mi; b1 = 1+si*(1-mi) 
  if(what=="mean"){(a1)/(a1+b1)}
  else if(what=="mode"){(a1-1)/(a1+b1-2)}
}

tablex = function(x,M){
  return(mapply(function(m)sum(x==m),1:M))
}


add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


