#' Multi-fractal Detrended Fluctuation Analysis
#'
#' @param signal    An input signal.
#' @param qq    A vector containing a range of values for the order of fluctuation \code{q}.
#' @param mins    Minimum scale to consider.
#' @param maxs    Maximum scale to consider.
#' @param ressc
#' @param m
#'
#' @return
#' @export
#'
#' @examples
MFDFA <- function(signal,qq,mins=6,maxs=12,ressc=30,m=1){
  require(pracma)
#   reload <- FALSE
#   if("signal" %in% .packages()){
#     warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
#     reload <- TRUE
#     detach("package:signal", unload=TRUE)
#   }
  scale     <- round(2^(seq(mins,maxs,by=((maxs-mins)/ressc))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  qRMS      <- vector("list",length(qq))
  Fq        <- vector("list",length(qq))
  qRegLine  <- vector("list",length(qq))
  Hq        <- numeric(length(qq))

  Y        <- cumsum(signal-mean(signal))
  TSm      <- as.matrix(cbind(t=1:length(Y),y=Y))
  Hglobal  <- monoH(TSm,scale)

  Hadj <- 0
  if((Hglobal>1.2)&(Hglobal<1.8)){
    Y <- diff(signal)
    Hadj=1}
  if(Hglobal>1.8){
    Y <- diff(diff(signal))
    Hadj <- 2}
  if(Hglobal<0.2){
    Y <- cumsum(signal-mean(signal))
    Hadj <- -1}
  if(Hadj!=0){TSm  <- as.matrix(cbind(t=1:length(Y),y=cumsum(Y-mean(Y))))}

  for(ns in seq_along(scale)){
    RMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
    for(nq in seq_along(qq)){
      qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)] <- RMS_scale[[ns]]$V1^qq[nq]
      Fq[[nq]][ns] <- mean(qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)])^(1/qq[nq])
      if(is.inf(log2(Fq[[nq]][ns]))){Fq[[nq]][ns]<-NA}
    }
    Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
    if(is.inf(log2(Fq[[which(qq==0)]][ns]))){Fq[[which(qq==0)]][ns]<-NA}
  }

  fmin<-1
  fmax<-which(scale==max(scale))
  #for(nq in seq_along(qq)){Hq[nq] <- lm(log2(Fq[[nq]])~log2(scale))$coefficients[2]}
  Hq <- ldply(Fq,function(Fqs){lm(log2(Fqs[fmin:fmax])~log2(scale[fmin:fmax]),na.action=na.omit)$coefficients[2]})

  tq <- (Hq[,1]*qq)-1
  hq <- diff(tq)/diff(qq)
  Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])

  monoH <- function(TSm,scale){
    dfaRMS_scale <- vector("list",length(scale))
    F2 <- numeric(length(scale))
    for(ns in seq_along(scale)){
      dfaRMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
      F2[ns] <- mean(dfaRMS_scale[[ns]]$V1^2)^(1/2)
      if(is.inf(log2(F2[ns]))){F2[ns] <- NA}
    }
    return(lm(log2(F2)~log2(scale),na.action=na.omit)$coefficients[2])
  }

  detRend <- function(TS, Order=1){
    detR <- lm(TS~stat::poly(1:length(TS), degree=Order))$residuals
    return(detR)
  }

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(q=qq,Hq=Hq,tq=tq,hq=hq,Dq=Dq,Hglobal=Hglobal,Hadj=Hadj))
}


set.seed(100)
z <- dispersion(rnorm(1024))
plot(log(z$scale),log(z$sd))
#

# trace(detRend,edit=T)
# seq(1,length(X),by=4096)
#
# z<-sliceTS(TSm,scale[1])
# z[[1]][,2]
#
# Hglobal <-
#
# segments <- laply(scale,function(s) floor(length(X)/s))
# IDv <- llply(segments,slice.index,)
# segv <- function(X,segments){
#
#
#
#   seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#
# }
# for(ns in seq_along(scale)){
#   segv[ns] <- floor(length(X)/scale[ns])
#   for(v in 1:segv[ns]){
#     Index <- seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#     Cslope<- polyfit(Index,X[Index],m)
#     fit   <- polyval(Cslope,Index)
#     RMS_scale[[ns]][v] <- sqrt(mean((X[Index]-fit)^2))
#     rm(Cslope, fit, Index)
#   }
#   for(nq in seq_along(qq)){
#     qRMS[[nq]][1:segv[ns]] <- RMS_scale[[ns]]^qq[nq]
#     Fq[[nq]][ns] <- mean(qRMS[[nq]][1:segv[ns]])^(1/qq[nq])
#   }
#   Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
# }
#
# for(nq in seq_along(qq)){
#   Cslope <- polyfit(log2(scale),log2(Fq[[nq]]),1)
#   Hq[nq] <- Cslope[1]
#   qRegLine[[nq]] <- polyval(Cslope,log2(scale))
#   rm(Cslope)
# }
#
# tq <- (Hq*qq)-1
# hq <- diff(tq)/diff(qq)
# Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])
#
# plot(hq,Dq,type="l")
#
#
# qq<-c(-10,-5,seq(-2,2,.1),5,10)
