#' FD.psd
#'
#' @description Get Fractal Dimension by log-log PSD slope estimate using a variety of spectral methods.
#'
#' @param y    A vector or matrix of time series objects.
#' @param methods    Character array with any, or all (default) of: "spec.pgram", "pwelch", "hurstSpec".
#' @param nfft    Number of frequencies to be estimated (only used with methods: "pwelch" and "hurstSpec")
#'
#' @author Fred Hasselman
#'
#' @return A list object with fields:
#'
#' @export
#'
#' @details
#'
#' @examples
fd.psd <- function(y, normalize = TRUE, detrend = TRUE){
#require(oce)
require(fractal)

n <- length(y[,1])
y <- scale(y)

# # Number of frequencies estimated cannot be set! (defaults to Nyquist)
# # Use Tukey window: cosine taper with r = 0.5
# # fast = TRUE ensures padding with zeros to optimize FFT to highly composite number, however, we just pad to nextPow2
# npad <- 1+(nextn(n,factors=2)-n)/n
# npad <- nextn(n)
# if(n==npad) npad = 0

#  psd  <- spec.pgram(ys, demean=F, detrend=T, plot=F, pad=, taper=0.5)

# First check the global slope for anti-persistent noise (GT +0.20)
# If so, fit the line starting from the highest frequency
nr     <- length(psd$freq)-1
lsfit  <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr]))
glob   <- coef(lsfit)[2]

# General guideline: fit over 25% frequencies
# If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
nr <- HDEst(NFT = length(psd$freq), as.vector(SDF(ys[,1])))

exp1 <- hurstSpec(y, sdf.method = "wosa", freq.max = 0.25 )
exp2 <- hurstSpec(y, sdf.method = "wosa", freq.max = psd$freq[nr] )

# psd <- pwelch(ys, nfft = 1024, noverlap = .5, taper = 0.5, plot = FALSE)

ifelse( (glob > 0.2),
  lsfit <- lm(log(rev(psd$spec)[2:nr]) ~ log(rev(psd$freq)[2:nr])),
  lsfit <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr])) )

return(coef(lsfit)[2])
}


require(plyr)
require(reshape2)
df.y <- ldply(seq(0.05,.95,by=.05), function(H) paste0('H.',H,'=',eval(parse(text=as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024))))))
df.y <- melt(df.y)
lmModel("dfbm", HB=H)

require(rio)
ts1 <- import('ts1.txt')
ts2 <- import('ts2.txt')
ts3 <- import('ts3.txt')



ys <- ts1

fgn<-list()
fBm<-list()
Hs <- seq(0.1, .9, by=.1)
cnt=0
for(H in Hs){
  cnt<-cnt+1
fgn[[cnt]] <- ts(as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024)),start = 0, end = 1, deltat=1/1024)
fgn[[cnt]] <- (fgn[[cnt]]-min(fgn[[cnt]]))/max(fgn[[cnt]])
fBm[[cnt]] <- ts(as.numeric(lmSimulate(lmModel("dfbm", HB=H), n.sample = 1024)), start = 0, end = 1, deltat=1/1024)
fBm[[cnt]] <- (fBm[[cnt]]-min(fBm[[cnt]]))/max(fBm[[cnt]])
}
names(fgn) <- paste0('H = ',Hs)
names(fBm) <- paste0('H = ',Hs+1)

op<-par(yaxt = 'n')
stackPlot(time(fBm[[1]]),c(fgn,fBm), rescale = T)
#stackPlot(1:1024,fBm, rescale = T)
par(op)

plot(ts(as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024)), start = 0,end = 1, deltat=1/1024))
plot(ts(as.numeric(lmSimulate(lmModel("dfbm", HB=H), n.sample = 1024)), start = 0, end = 1, deltat=1/1024))

models <- c("ppl","fdp","fgn")
lag <- 100
z <- lapply(models, function(x, models, lag)
{ lmACF(lmModel(x), lag=lag)@data},
models=models, lag=lag)
names(z) <- paste(upperCase(models), "ACF")
stackPlot(seq(0,lag), z, xlab="lag")
title("Stochastic Fractal Model ACFs")


