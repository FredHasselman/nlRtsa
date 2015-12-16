#' fd.psd
#'
#' @description Get Fractal Dimension by log-log PSD slope estimate using a variety of spectral methods.
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix with columns \code{freq.norm}, \code{freq} and \code{spec}.
#' \item Model output in field \code{lmfit}.
#' \item Estimate of scaling exponent \code{alpha} based on a fit over the lowest 25% frequencies (\code{low25}), or using the HD estimate \{HD}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' }
#'
#' @export
#'
#' @details Calls function \code{\link{SDF}} to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), \code{SDF} is called using a Tukey window (\code{raised cosine \link{taper}).
#'
#' A line is fitted on the periodogram in log-log coordinates. Two fit-ranges are used: The 25% lowest frequencies and the Hurvic-De
#'
#' @examples
fd.psd <- function(y, fs = NULL, normalize = TRUE, dtrend = TRUE, plot = FALSE){
require(pracma)
require(fractal)
require(sapa)

if(!is.ts(y)){
  if(is.null(fd)){fs <- 1}
  y <- ts(y, frequency = fs)
  cat("\n\nFD.psd:\tSample rate was set to 1.\n\n")
}

N             <- length(y)
# Simple linear detrending.
if(dtrend)    y <- ts(detrend(as.vector(y), tt = 'linear'), frequency = fs)
# Normalize using N instead of N-1.
if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

# Number of frequencies estimated cannot be set! (defaults to Nyquist)
# Use Tukey window: cosine taper with r = 0.5

# fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
# However, we just pad to nextPow2, except if length already is a power of 2.
npad <- 1+(stats::nextn(N,factors=2)-N)/N
npad <- stats::nextn(N)

# if(N==npad) npad = 0
# psd  <- stats::spec.pgram(y, fast = FALSE, demean=FALSE, detrend=FALSE, plot=FALSE, pad=npad, taper=0.5)

Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = npad)
psd   <- sapa::SDF(y, taper. = Tukey, npad = npad)

info   <- capture.output(print(psd))
df.psd <- as.matrix(psd)

powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1], freq = attr(psd, "frequency")[-1]*frequency(y), spec = as.matrix(psd)[-1])

# First check the global slope for anti-persistent noise (GT +0.20)
# If so, fit the line starting from the highest frequency
nr     <- length(psd$freq)-1
lsfit  <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr]))
glob   <- coef(lsfit)[2]

# General guideline: fit over 25% frequencies
# If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
nr <- fractal::HDEst(NFT = length(psd$freq), as.vector(SDF(y)))

exp1 <- hurstSpec(y, sdf.method = "wosa", freq.max = 0.25 )
exp2 <- hurstSpec(y, sdf.method = "wosa", freq.max = psd$freq[nr] )

# psd <- pwelch(ys, nfft = 1024, noverlap = .5, taper = 0.5, plot = FALSE)

ifelse( (glob > 0.2),
  lsfit <- lm(log(rev(psd$spec)[2:nr]) ~ log(rev(psd$freq)[2:nr])),
  lsfit <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr])) )

return(coef(lsfit)[2])
}

FD.sda <- function(y, normalize = TRUE){

  if(normalize){y <- scale(y)}

  out           <- dispersion(y)
  maxbin        <- length(out$sd)-2
  lsfit         <- lm(log(out$sd[2:maxbin]) ~ log(out$scale[2:maxbin]))
  FDdata <- c(coef(lsfit)[2],sda2fd(coef(lsfit)[2]))
  }

FD.sda <- function(y, dmethod = "poly1"){
  out           <- DFA(y[,k], detrend=dmethod, sum.order=1, scale.max=trunc(length(y)/4), scale.min=4, scale.ratio=2^(1/4), verbose=FALSE)
  FDdata <- c(attributes(out)$logfit[]$coefficients['x'],dfa2fd(attributes(out)$logfit[]$coefficients['x']))
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


