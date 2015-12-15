#' Dpsd
#'
#' @description Get PSD slope by Periodogram
#'
#' @param ys
#'
#' @return
#' @export
#'
#' @examples
Dpsd <- function(ys)
{
n <- length(ys)
# Number of frequencies estimated cannot be set! (defaults to Nyquist)
# Use Tukey window: cosine taper with r = 0.5
# fast = TRUE ensures padding with zeros to optimize FFT to highly composite number, however, we just pad to nextPow2
npad <- 1+(nextn(n,factors=c(2))-n)/n
psd  <- spec.pgram(ys, demean=F, detrend=F, plot=F, pad=npad, taper=0.5)

# First check the global slope for anti-persistent noise (GT +0.20)
# If so, fit the line starting from the highest frequency
nr     <- length(psd$freq)-1
lsfit  <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr]))
glob   <- coef(lsfit)[2]

# General guideline: fit over 25% frequencies
# If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
nr <- trunc((length(psd$freq)) * 0.25)
n  <- length(psd$freq)

ifelse( (glob > 0.2),
  lsfit <- lm(log(rev(psd$spec)[2:nr]) ~ log(rev(psd$freq)[2:nr])),
  lsfit <- lm(log(psd$spec[2:nr]) ~ log(psd$freq[2:nr])) )

return(coef(lsfit)[2])
}
