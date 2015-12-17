#' fd.psd
#'
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
#' \item A data matrix \code{PSD} with columns \code{freq.norm}, \code{freq} and \code{spec}.
#' \item Model output in field \code{lmfit}.
#' \item Estimate of scaling exponent \code{alpha} based on a fit over the lowest 25\% frequencies (\code{low25}), or using the HD estimate \code{HD}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
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
        plot(log(powspec$spec) ~ log(powspec$freq), xlab="log(Frequency)", ylab = "log(Power)")
        lines(log(powspec$freq[powspec$freq<=0.25]), predict(lmfit1),lwd=3,col="darkred")
        lines(log(powspec$freq[1:nr]), predict(lmfit2),lwd=3,col="darkblue")
        legend("bottomleft",c(paste0("lowest 25% (n = ",sum(powspec$freq<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
        par(old)
    }

    return(list(
        PLAW  = powspec,
        low25 = list(sap = coef(lmfit1)[2], H = exp1, FD = psd2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
        HD    = list(sap = coef(lmfit2)[2], H = exp2, FD = psd2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
        info  = psd)
    )
}

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
#' @return Slope and FD estimate based on Hasselman(2013)
#' @export
#'
#' @family FD estimators
#' @examples
fd.sda <- function(y, fs = NULL, normalize = TRUE, dtrend = FALSE, scales = dispersion(y)$scale, fitRange = c(scales[2], scales[length(scales)-1]), plot = FALSE){
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
        plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(psd2fd(coef(lmfit1)[2]),digits=2),'\nRange    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(psd2fd(coef(lmfit2)[2]),digits=2)))
        ifultools::splitplot(2,1,2)
        plot(log(out$sd) ~ log(out$scale), xlab="log(Bin Size)", ylab = "log(SD)")
        lines(log(out$scale), predict(lmfit1),lwd=3,col="darkred")
        lines(log(out$scale[bins]), predict(lmfit2),lwd=3,col="darkblue")
        legend("topright",c(paste0("Full (n = ",length(out$scale),")"), paste0("Range (n = ",length(bins),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
        par(old)
    }

    return(list(
        PLAW  =  cbind.data.frame(freq.norm = frequency(y)/scales, size = out$scale, bulk = out$sd),
                                  fullRange = list(sap = coef(lmfit1)[2], H = 1+coef(lmfit1)[2], FD = sda2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
                                  fitRange  = list(sap = coef(lmfit1)[2], H = 1+coef(lmfit2)[2], FD = sda2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
                                  info = out)
    )
}

#' fd.dfa
#'
#' @title Detrended Fluctuation Analysis (DFA)
#'
#' @param y
#' @param dmethod
#'
#'
#' @return Estimate of Hurst exponent (slope of \code{log(bin)} vs. \code{log(RMSE))} and an FD estimate based on Hasselman(2013)
#'
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family FD estimators
#' @examples
fd.dfa <- function(y, fs = NULL, dtrend = "poly1", normalize = FALSE, sum.order = 1, scale.max=trunc(length(y)/4), scale.min=4, scale.ratio=2^(1/4), overlap = 0.5, plot = FALSE){
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
        old<- ifultools::splitplot(2,1,1)
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
        PLAW  =  cbind.data.frame(freq.norm = scale.R(attributes(out)$scale*frequency(y)), size = attributes(out)$scale, bulk = attributes(out)$stat),
                                  fullRange = list(sap = coef(lmfit1)[2], H = attributes(out)$logfit[]$coefficients['x'] , FD = dfa2fd(coef(lmfit1)[2]), fitlm1 = lmfit1),
                                  fitRange  = list(sap = coef(lmfit2)[2], H = coef(lmfit2)[2], FD = dfa2fd(coef(lmfit2)[2]), fitlm2 = lmfit2),
                                  info = out)
        )
}


#' psd2fd
#'
#' @description Conversion formula: From periodogram based self-affinity parameter estimate (\code{sap}) to an informed estimate of the (fractal) dimension (FD).
#' @param sap Self-afinity parameter estimate based on PSD slope (e.g., \code{\link{fd.psd}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @details The spectral slope will be converted to a dimension estimate using:
#'
#' \deqn{D_{PSD}\approx\frac{3}{2}+\frac{14}{33}*\tanh\left(Slope * \ln(1+\sqrt{2})\right)}{D_{PSD} ≈ 3/2 + ((14/33)*tanh(sap*log(1+sqrt(2))))}
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#' @examples
#' # Informed FD of white noise
#' psd2fd(0)
#'
#' # Informed FD of Brownian noise
#' psd2fd(-2)
#'
#' # Informed FD of blue noise
#' psd2fd(2)
psd2fd <- function(sap){return(round(3/2 + ((14/33)*tanh(sap*log(1+sqrt(2)))), digits = 2))}


#' DFA slope (H in DFA)
#'
#' @description Conversion formula: Detrended Fluctuation Analysis (DFA) estimate of the Hurst exponent (a self-affinity parameter \code{sap}) to an informed estimate of the (fractal) dimension (FD).
#'
#' @param sap Self-afinity parameter estimate based on DFA slope (e.g., \code{\link{fd.sda}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#'
#' @export
#'
#' @details The DFA slope (H) will be converted to a dimension estimate using:
#'
#' \deqn{D_{DFA}\approx 2-(\tanh(\log(3)*sap)) }{D_{DFA} ≈ 2-(tanh(log(3)*sap)) }
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @examples
#' # Informed FD of white noise
#' dfa2fd(0.5)
#'
#' # Informed FD of Pink noise
#' dfa2fd(1)
#'
#' # Informed FD of blue noise
#' dfa2fd(0.1)
dfa2fd <- function(sap){return(round(2-(tanh(log(3)*sap)), digits = 2))}



#' sda2fd
#'
#' @description Conversion formula: Standardised Dispersion Analysis (SDA) estimate of self-affinity parameter (\code{sap}) to an informed estimate of the (fractal) dimension (FD).
#'
#' @param sap Self-afinity parameter estimate based on SDA slope (e.g., \code{\link{fd.sda}})).
#'
#' @details
#'
#'
#' Note that for some signals different PSD slope values project to a single SDA slope. That is, SDA cannot distinguish between all variaties of power-law scaling in the frequency domain.
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @examples
#' # Informed FD of white noise
#' sda2fd(-0.5)
#'
#' # Informed FD of Brownian noise
#' sda2fd(-1)
#'
#' # Informed FD of blue noise
#' sda2fd(-0.9)
sda2fd <- function(sap){return(1-sap)}




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

