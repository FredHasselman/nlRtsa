opp <- function (n,degree){
  basis <- sapply(0:degree,function(i) (1:n)^i)
  Q <- qr.Q(qr(basis),complete=T)
  ML_proj = Q[,1:ncol(basis)]
  ML_perp = Q[,(ncol(basis)+1):ncol(Q)]
  P = ML_proj%*%t(ML_proj)
  P_perp = ML_perp%*%t(ML_perp)
  basis = ML_proj
  return(P_Perp)
}



binTS<-function(TSmat,binSz=16,Overlap=0.5) {
  # Slice columns of TSmat in epochs of size = Scales
  require("plyr")
  if((Overlap<0)&(Overlap>=1)){break}
  N <- dim(TSmat)
  return(llply(seq(1,round(N[1]-((1-Overlap)*binSz)),by=round((1-Overlap)*binSz)),function(i) TSmat[i:min(i+binSz-1,N[1]),1:N[2]]))
}

detRES <- function(binTS, Order=1){
  return(lm(binTS~poly(1:length(binTS), degree=Order, raw=T))$residuals)
}

getScales <- function(TSvec, Fs=500, tsec=c(tmin=(4/Fs),tmax=(length(TSvec)/Fs)/4), Ns=20){
  return(Scales<-unique(round(2^seq(log2(tsec[1]*Fs),log2(tsec[2]*Fs),length.out=Ns))))
}


FqS <- function(fluctMEAN,Q){
# For each Q, return detrended fluctuation
# Adjust Q = 0
 require("plyr")   
  FQS <- aaply(Q,1,function(q) (mean(abs(fluctMEAN)^(q/2)))^(1/q))
  if(0 %in% Q) { FQS[Q==0] <- exp( (1/(2*length(fluctMEAN)))*(sum(log(abs(fluctMEAN)))) ) }
  names(FQS) <- Q
 return(FQS)
}

mfDXA <- function(TSmatX, TSmatY, Scales=getScales(TSmatX), Q=-10:10, Order=1){
  # For each Q, and Scale, return detrended covariance-based fluctuation
  require("plyr")
  
  # This should not be a problem, just to be safe, unload package:signal because of poly()  
  lp <- search()
  lsignal <- FALSE
  if("package:signal" %in% lp){
    lsignal <- TRUE
    detach(package:signal)}

  names(Scales) <- paste("s=", Scales,", v=",sep="")
  binTSmatX     <- llply(Scales,binTS,TSmat=TSmatX)
  str           <- as.numeric(summary(binTSmatX)[,1])
  
  detTSmatX <- llply(unlist(binTSmatX, recursive=F, use.names=T), function(x){apply(x, 2, detRES, Order=Order)})
  nScale    <- factor(rep(1:length(str), times = str))
  
  # Check if 2 matrices were provided
  ifelse(!is.null(TSmatY),{
    if(any(dim(TSmatX)!=dim(TSmatY))){break}
    binTSmatY <- llply(Scales,binTS,TSmat=TSmatY)
    detTSmatY <- llply(unlist(binTSmatY, recursive=F, use.names=T), function(x){apply(x,2, detRES, Order=Order)})
    XY = TRUE},{
      detTSmatY <- detTSmatX
      XY = FALSE} )
  
  covSVmat <- mapply(cov,detTSmatX,detTSmatY,SIMPLIFY=F)
  F2s      <- ldply(unique(nScale), function(i) {cbind(Location=dimnames(TSmatX)[[2]],Reduce('+', covSVmat[nScale==i])/length(covSVmat[nScale==i]))})
  rownames(F2s)<- paste(rep(paste("F2(",unlist(Scales),")",sep=""),each=dim(TSmatX)[2]),dimnames(TSmatX)[[2]],sep=" - ")
  
  plot(log2(Scales), log2(as.numeric(F2s[which(F2s$Location==names(F2s)[2]),2])))
  Hq <- ldply(1:dim(TSmatX)[[2]], function(i) apply(F2s[,2:dim(F2s)[[2]]],2, function(x,i=i) lm(log2(as.numeric(x[which(F2s$Location==names(F2s)[i+1]),i+1]))~poly(log2(Scales), degree=1))$coefficients[2])
  
    dim(TSmatX)[[2]]
    
  tst2     <- unlist(F2s, recursive=F, use.names=T)
  <- llply(F2s,function(x){x[lower.tri(x, diag=XY)] <- NA}
    names(xx[upper.tri(F2s[[1]],diag=XY)])
    
    id<-cbind(rep(1:9,each=26),rep(1:26,9))
    F2s[id[,1]==1&id[,1]==2]
    
    cadd <- function(x) Reduce("+", x, accumulate = TRUE)
    cadd(seq_len(7))
    
    F2s[[1]][2,]
    
    fluctQdep <- llply(1:length(str), function(x) FqS(fluct$detCOV[fluct$Scale==x], Q=Q)) 
    return(MFDXA<-list(Fs=fluctBINS,FqS=fluctQdep))
    
    tst<-cov(detTSmatX[[1]][,2],detTSmatY[[1]][,1])
    cov(detTSmatX[[1]],detTSmatY[[1]])
    mean(tst[,1])
    mean(tst[1,])
    
    
    l_ply(1:9, function(x){ plot(Q,fluctQdep[[x]]); par(new=F) })
    
    if(lsignal==TRUE){
      attach(package:signal)
      lsignal==FALSE}
}