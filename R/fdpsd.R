fd <- function(y, ...) UseMethod("fd")

fd.default <- function(y, ...){

    cat("No type specified.\nReturning exponential growth power law.")

    r = 1.01
    y <- growth.ac(Y0=0.001, r=r, N=2048, type = "driving")
    tsp(y) <-c(1/500,2048/500,500)
    bulk <- log1p(hist(y,plot = F, breaks = seq(0,max(y),length.out = 129))$counts)
    size <- log1p(seq(0,2047,length.out = 128))
    id<-bulk==0

    lmfit <- lm(bulk[!id] ~ size[!id])

    old <- ifultools::splitplot(2,1,1)
    plot(y, ylab = "Y", main = paste0('Exponential growth  sap: ', round(coef(lmfit)[2],digits=2), ' | r:', r))
    ifultools::splitplot(2,1,2)
    plot(size[!id],bulk[!id], xlab="Size = log(bin(Time))", ylab = "Bulk = logbin(Y)", pch=21, bg="grey60", pty="s")
    lines(size[!id], predict(lmfit),lwd=4,col="darkred")
    #legend("bottomleft",c(paste0("Range (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
}

#' @title Power Spectral Density Slope (PSD).

#' @description Estimate Alpha, Hurst Exponent and Fractal Dimension through log-log slope.
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{alpha} based on a fit over the lowest 25\% frequencies (\code{low25}), or using the HD estimate \code{HD}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @family FD estimators
#'
#' @export
#'
#' @details Calls function \code{\link[sapa]{SDF}} to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), \code{SDF} is called using a Tukey window (\code{raised cosine \link[sapa]{taper}}).
#'
#' A line is fitted on the periodogram in log-log coordinates. Two fit-ranges are used: The 25\% lowest frequencies and the Hurvich-Deo estimate (\code{\link[fractal]{HDEst}}).
#'
#' @examples
#' fd.psd(rnorm(2048), plot = TRUE)
fd.psd <- function(y, fs = NULL, normalize = TRUE, dtrend = TRUE, plot = FALSE){
    require(pracma)
    require(fractal)
    require(sapa)
    require(ifultools)

    if(!is.ts(y)){
        if(is.null(fs)){fs <- 1}
        y <- ts(y, frequency = fs)
        cat("\n\nFD.psd:\tSample rate was set to 1.\n\n")
    }

    N             <- length(y)
    # Simple linear detrending.
    if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
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

    powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1], size = attr(psd, "frequency")[-1]*frequency(y), bulk = as.matrix(psd)[-1])

    # First check the global slope for anti-persistent noise (GT +0.20)
    # If so, fit the line starting from the highest frequency
    nr     <- length(powspec[,1])
    lsfit  <- lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
    glob   <- coef(lsfit)[2]

    # General guideline: fit over 25% frequencies
    # If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
    nr <- fractal::HDEst(NFT = length(powspec[,1]), sdf = psd)

    exp1 <- hurstSpec(y, sdf.method = "direct", freq.max = 0.25, taper. = Tukey )
    exp2 <- hurstSpec(y, sdf.method = "direct", freq.max = powspec$freq.norm[nr], taper. = Tukey)

    ifelse((glob > 0.2), {
        lmfit1 <- lm(log(rev(powspec$bulk[powspec$size<=0.25])) ~ log(rev(powspec$size[powspec$size<=0.25])))
        lmfit2 <- lm(log(rev(powspec$bulk[1:nr])) ~ log(rev(powspec$size[1:nr])))
    },{
        lmfit1 <- lm(log(powspec$bulk[powspec$size<=0.25]) ~ log(powspec$size[powspec$size<=0.25]))
        lmfit2 <- lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
    })

    if(plot){
        old<- ifultools::splitplot(2,1,1)
        plot(y,ylab = "Y", main = paste0('Lowest 25%    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(exp1,digits=2), ' | FD:',round(psd2fd(coef(lmfit1)[2]),digits=2),'\nHurvic-Deo    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(exp2,digits=2), ' | FD:',round(psd2fd(coef(lmfit2)[2]),digits=2)))
        ifultools::splitplot(2,1,2)
        plot(log(powspec$bulk) ~ log(powspec$size), xlab="log(Frequency)", ylab = "log(Power)")
        lines(log(powspec$size[powspec$size<=0.25]), predict(lmfit1),lwd=3,col="darkred")
        lines(log(powspec$size[1:nr]), predict(lmfit2),lwd=3,col="darkblue")
        legend("bottomleft",c(paste0("lowest 25% (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
        par(old)
    }

    return(list(
        PLAW  = powspec,
        low25 = list(sap = coef(lmfit1)[2], H = exp1, FD = psd2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
        HD    = list(sap = coef(lmfit2)[2], H = exp2, FD = psd2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
        info  = psd)
    )
}


# SDA -------------------------------------------------

#' fd.sda
#'
#' @title Standardised Dispersion Analysis (SDA).
#'
#' @param y
#' @param normalize
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @family FD estimators
#' @examples
fd.sda <- function(y, fs = NULL, normalize = TRUE, dtrend = FALSE, scales = dispersion(y)$scale, fitRange = c(scales[1], scales[length(scales)-2]), plot = FALSE){
    require(pracma)
    require(fractal)

    if(!is.ts(y)){
        if(is.null(fs)){fs <- 1}
        y <- ts(y, frequency = fs)
        cat("\n\nfd.sda:\tSample rate was set to 1.\n\n")
    }

    N             <- length(y)
    # Simple linear detrending.
    if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
    # Normalize using N instead of N-1.
    if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

    bins          <- which(fitRange[1]==scales):which(fitRange[2]==scales)
    out           <- dispersion(y, front = FALSE)
    lmfit1        <- lm(log(out$sd) ~ log(out$scale))
    lmfit2        <- lm(log(out$sd[bins]) ~ log(out$scale[bins]))

    if(plot){
        old<- ifultools::splitplot(2,1,1)
        plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sda2fd(coef(lmfit1)[2]),digits=2),'\nRange    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sda2fd(coef(lmfit2)[2]),digits=2)))
        ifultools::splitplot(2,1,2)
        plot(log(out$sd) ~ log(out$scale), xlab="log(Bin Size)", ylab = "log(SD)")
        lines(log(out$scale), predict(lmfit1),lwd=3,col="darkred")
        lines(log(out$scale[bins]), predict(lmfit2),lwd=3,col="darkblue")
        legend("bottomleft",c(paste0("Full (n = ",length(out$scale),")"), paste0("Range (n = ",length(bins),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
        par(old)
    }

    return(list(
        PLAW  =  cbind.data.frame(freq.norm = frequency(y)/scales, size = out$scale, bulk = out$sd),
                                  fullRange = list(sap = coef(lmfit1)[2], H = 1+coef(lmfit1)[2], FD = sda2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
                                  fitRange  = list(sap = coef(lmfit2)[2], H = 1+coef(lmfit2)[2], FD = sda2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
                                  info = out)
    )
}


# DFS ---------------------------------------------

#' fd.dfa
#'
#' @title Detrended Fluctuation Analysis (DFA)
#'
#' @param y
#' @param dmethod
#'
#'
#' @return Estimate of Hurst exponent (slope of \code{log(bin)} vs. \code{log(RMSE))} and an FD estimate based on Hasselman(2013)
#' A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family FD estimators
#' @examples
fd.dfa <- function(y, fs = NULL, dtrend = "poly1", normalize = FALSE, sum.order = 1, scale.max=trunc(length(y)/4), scale.min=4, scale.ratio=2^(1/4), overlap = 0, plot = FALSE){
    require(pracma)
    require(fractal)

    if(!is.ts(y)){
        if(is.null(fs)){fs <- 1}
        y <- ts(y, frequency = fs)
        cat("\n\nfd.dfa:\tSample rate was set to 1.\n\n")
    }

    N             <- length(y)
    # Normalize using N instead of N-1.
    if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

    out1 <- DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=trunc(length(y)/2), scale.min=2, scale.ratio=2, overlap = 0, verbose=FALSE)
    out2 <- DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=scale.max, scale.min=scale.min, scale.ratio=scale.ratio, overlap = overlap, verbose=FALSE)

    lmfit1        <- lm(log(attributes(out1)$stat) ~ log(attributes(out1)$scale))
    lmfit2        <- lm(log(attributes(out2)$stat) ~ log(attributes(out2)$scale))

    if(plot){
        plot.new()
        old <- ifultools::splitplot(2,1,1)
        plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:',
                                         round(attributes(out1)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                         round(dfa2fd(coef(lmfit1)[2]),digits=2),'\nRange    sap: ',
                                         round(coef(lmfit2)[2],digits=2), ' | H:',
                                         round( attributes(out2)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                         round(dfa2fd(coef(lmfit2)[2]),digits=2)
                                         )
             )
        ifultools::splitplot(2,1,2)
        plot(log(attributes(out1)$stat) ~ log(attributes(out1)$scale), xlab="log(Bin Size)", ylab = "log(RMSE)")
        lines(log(attributes(out1)$scale), predict(lmfit1),lwd=3,col="darkred")
        lines(log(attributes(out2)$scale), predict(lmfit2),lwd=3,col="darkblue")
        legend("topleft",c(paste0("Full (n = ",length(attributes(out1)$scale),")"), paste0("Range (n = ",length(attributes(out2)$scale),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
        par(old)
    }
    return(list(
        PLAW  =  cbind.data.frame(freq.norm = scale.R(attributes(out1)$scale*frequency(y)), size = attributes(out1)$scale, bulk = attributes(out1)$stat),
                                  fullRange = list(sap = coef(lmfit1)[2], H = attributes(out1)$logfit[]$coefficients['x'] , FD = dfa2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
                                  fitRange  = list(sap = coef(lmfit2)[2], H = coef(lmfit2)[2], FD = dfa2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
                                  info = list(out1,out2))
        )
}


#
#
# require(plyr)
# require(reshape2)
# df.y <- ldply(seq(0.05,.95,by=.05), function(H) paste0('H.',H,'=',eval(parse(text=as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024))))))
# df.y <- melt(df.y)
# lmModel("dfbm", HB=H)
#
# require(rio)
# ts1 <- import('ts1.txt')
# ts2 <- import('ts2.txt')
# ts3 <- import('ts3.txt')
#
#
#
# ys <- ts1
#
# fgn<-list()
# fBm<-list()
# Hs <- seq(0.1, .9, by=.1)
# cnt=0
# for(H in Hs){
#   cnt<-cnt+1
# fgn[[cnt]] <- ts(as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024)),start = 0, end = 1, deltat=1/1024)
# fgn[[cnt]] <- (fgn[[cnt]]-min(fgn[[cnt]]))/max(fgn[[cnt]])
# fBm[[cnt]] <- ts(as.numeric(lmSimulate(lmModel("dfbm", HB=H), n.sample = 1024)), start = 0, end = 1, deltat=1/1024)
# fBm[[cnt]] <- (fBm[[cnt]]-min(fBm[[cnt]]))/max(fBm[[cnt]])
# }
# names(fgn) <- paste0('H = ',Hs)
# names(fBm) <- paste0('H = ',Hs+1)
#
# op<-par(yaxt = 'n')
# stackPlot(time(fBm[[1]]),c(fgn,fBm), rescale = T)
# #stackPlot(1:1024,fBm, rescale = T)
# par(op)
#
# plot(ts(as.numeric(lmSimulate(lmModel("fgn", HG=H), n.sample = 1024)), start = 0,end = 1, deltat=1/1024))
# plot(ts(as.numeric(lmSimulate(lmModel("dfbm", HB=H), n.sample = 1024)), start = 0, end = 1, deltat=1/1024))
#
# models <- c("ppl","fdp","fgn")
# lag <- 100
# z <- lapply(models, function(x, models, lag)
# { lmACF(lmModel(x), lag=lag)@data},
# models=models, lag=lag)
# names(z) <- paste(upperCase(models), "ACF")
# stackPlot(seq(0,lag), z, xlab="lag")
# title("Stochastic Fractal Model ACFs")
#

